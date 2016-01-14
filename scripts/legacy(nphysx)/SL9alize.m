%% A script to analyze SL9 output (dump.xml or dump.mat files)
% this is work in progress.
clear
clc
addpath ../scripts

%% Load the scene
% The dumps from all PhysX projects are in .xml files. The SDK provides
% a library with methods to serialize the state of the simulation. The
% resutling  structures are not very easy to decipher but all the
% information is there - in xml format. MATLAB has a native xmlread
% function, but it returns a DOM object which is just as hard to work with
% as the original file. The excellent xml_read from the file exchange comes
% to the rescue and returns a usable tree structure. I like to save some
% individual sub-trees into separate variables to save some typing.
% EDIT: xml_read becomes too slow on large files so i had to learn to
% navigate a DOM object anyway. Now I use readNXU. Also note that DOM
% parsers are so memory inefficient that the java heap size needs to be
% increased in matlab.
[filename pathname]=uigetfile({'*.mat';'*.xml'},[],'..\Release\');
if ~ischar(filename), return, end
filename=fullfile(pathname,filename);
if (filename)
    [~,~,xt]=fileparts(filename);
    if strcmp(xt,'.xml')
        tree=readNXU(filename); % this can take a minute...
    elseif strcmp(xt,'.mat')
        load(filename); % or maybe you have the tree in a mat file already
    else
        error('don''t be a dick');
    end
    clear xt
end
if ~exist('tree','var') || ~isfield(tree,'NxuPhysicsCollection')
    error('don''t be a dick');
end
actors=tree.NxuPhysicsCollection.NxSceneDesc.NxActorDesc;

% get x,y,z etc variables in usable form
globalPose=cell2mat({actors.globalPose}');
xVec=globalPose(:,10);
yVec=globalPose(:,11);
zVec=globalPose(:,12);
if isfield(actors,'NxSphereShapeDesc')
    dummy=[actors.NxSphereShapeDesc];
    dummy=[dummy.ATTRIBUTE];
    rVec=[dummy.radius];
elseif isfield(actors,'NxConvexShapeDesc')
    dummy=[actors.NxConvexShapeDesc];
    dummy=[dummy.ATTRIBUTE];
    meshLabelVec={dummy.meshData}'; % label of mesh used by actors
    % get convex shapes in usable form: this is quite brittle due to the
    % unpredictable return type from xml_read - it sometimes defaults to
    % char...
    shapes=cell(numel(actors),1);
    rVec=nan(size(xVec));
    vVec=nan(size(xVec));
    for k=1:numel(actors)
        mesh_id=str2double(meshLabelVec{k}(12:end));
        xmlmesh=...
            tree.NxuPhysicsCollection.NxConvexMeshDesc(mesh_id+1).points';
        shapes{k}=reshape(xmlmesh,3,[]); % but it's not rotated yet!
        pose=reshape(globalPose(k,1:9),3,3)';
        shapes{k}=pose*shapes{k}; % now it's rotated
        shapes{k}=shapes{k}'; % and make it a nice column major thing
        rVec(k)=mean(diag(shapes{k}*shapes{k}'));
        [~,vVec(k)]=convhulln(shapes{k});
    end
    rVec=sqrt(rVec);
    rPart=mean(rVec);
end
rVec=rVec(:); % make sure it's column

% look for a matching orbit file
orbfile=dir([pathname '*.orb']);
if ~isempty(orbfile)
    orbit=importdata([pathname orbfile.name]);
end

%% Find clusters
% I use eclazz to find the equivalence classes of the dataset of point
% coordinates with the equivalence relation inClust such that particle i is
% inClust of particle j if the distance between their centers is less than
% the sum of their radii plus a small "skin" correction. This is clearly
% reflexive and symetric, and eclazz assumes it is transitive so it doesn't
% have to be defined as such explicitly.
skinFactor=1.2;
if exist('vVec','var'), rVec=vVec.^(1/3); end % for clustering size of a poly is not its side
dset=num2cell([xVec,yVec,zVec,rVec],2); % eclazz expects cell array data
inClust=@(v1,v2)norm(v1(1:3)-v2(1:3))<(v1(4)+v2(4))*skinFactor; % equivalence R
clustVec=eclazz(inClust,dset); % NOTE: may take a few minutes

%% View the scene
% This section uses scatter or scatter3 to get a basic look at the scene.
% The way to visualize the scene depeneds on the shapes involved. The
% simple branch deals with spheres.

% the most basic scatter. since the disrupted comet is expected to have
% very little extent outside the orbital plane i will use 2d scatter plot.
% the tricky part is to scale the size of plotted circles in the same
% dimensions as the space variables. the method i use here is a brute force
% construction of scatter plot from multiple plots of individual circles.
% the advantage is that it should handle zoom and pan without losing scale.
fh=figure; ah=axes;
axis equal
hold on
hVec=nan(size(xVec));
for k=1:length(xVec)
    hVec(k)=rectangle('position',...
        [(xVec(k)-rVec(k)),(yVec(k)-rVec(k)),rVec(k)*2,rVec(k)*2]/1000,...
        'curvature',1,'facecolor','k'); % a circle, in km
end

% pretty up the plot
set(ah,'xdir','reverse') % to match the PhysX view
set(ah,'box','on')
xlabel('X [km]')
ylabel('Y [km]')
xlim manual
ylim manual

% try to point the way to jupi. the string of pearls (sop) center of mass
% is approximately at (0,0). the negative of the orbit data therefore
% points all the way to jupi. i'm going to put an arrow on the screen
% pointing the way but with a convenient size and location.
if exist('orbit','var')
    toJupi=-orbit.data(end,2:3);
    sopExtent=sqrt(range(xVec)^2+range(yVec)^2);
    toJupi=toJupi/norm(toJupi);
    toJupi=[toJupi*0.4*sopExtent;toJupi*0.2*sopExtent];
    toJupi=toJupi/1000;
    hArrow=quiver(toJupi(1,1),toJupi(1,2),toJupi(2,1),toJupi(2,2));
end

%% Mark clusters on the plot
% let's mark clusters with a big red circle, but only if they contain more
% than nPart/100 particles
clusts=unique(clustVec);
hClust=nan(size(clusts));
clumpMassVec=nan(size(clusts));
for k=1:numel(clusts)
    elmX=xVec(clustVec==clusts(k));
    elmY=yVec(clustVec==clusts(k));
    elmZ=zVec(clustVec==clusts(k));
    elmR=rVec(clustVec==clusts(k));
    elmV=elmR.^3;
    Vtot=sum(elmV);
    if (numel(elmX))>(numel(xVec)/100)
        cmX=dot(elmX,elmV)/Vtot;
        cmY=dot(elmY,elmV)/Vtot;
        cmZ=dot(elmZ,elmV)/Vtot;
        cmR=max([range(elmX) range(elmY) range(elmZ)])/2;
        hClust(k)=rectangle('position',...
            [(cmX-cmR)/1000,(cmY-cmR)/1000,cmR/500,cmR/500],'curvature',1);
    end
    clumpMassVec(k)=numel(elmX);
end
hClust(isnan(hClust))=[];
set(hClust,'linewidth',2,'edgecolor','r')

%% Annotate plot
% since rectangles don't show up in legends (they really don't, see
% MathWorks solution ID 1-5F9YKM) we create invisible line objects with the
% same appearance as the grains and clusters just for the legend.
hFakeGrain=line(nan,nan,'linewidth',2,'linestyle','none',...
    'marker','.','color','k','displayname',['Grains (' num2str(numel(xVec)) ')']);
hFakeClust=line(nan,nan,'linewidth',2,'linestyle','none',...
    'marker','o','color','r','displayname',['Clumps (' num2str(numel(hClust)) ')']);
if exist('hArrow','var')
    set(hArrow,'displayName','to Jupiter')
else
    hArrow=line(nan,nan,'displayname',' ');
end
legend([hFakeGrain hFakeClust hArrow],'location','northwest')
textbp(['largest fractional clump = ' num2str(max(clumpMassVec)/numel(xVec))])

% try to divine run_base_name from dump directory
s=dir([fileparts(filename) filesep '*.mat']);
if ~isempty(s), title(s(1).name(1:end-4),'interpreter','none'), end
