function pData=kreutzalize(filename)

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

if ~ischar(filename), return, end

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
linearVelocity=cell2mat({actors.linearVelocity}');
mass=cell2mat({actors.mass}');
X=globalPose(:,10);
Y=globalPose(:,11);
Z=globalPose(:,12);
VX=linearVelocity(:,1);
VY=linearVelocity(:,2);
VZ=linearVelocity(:,3);

% return a single array
pData=[X Y Z VX VY VZ mass];

end

function tree=readNXU(filename)
%READNXU make a tree structure from the xml saved by PhysX

%% create the dom object.
[~,~,xt]=fileparts(filename);
if ~strcmp(xt,'.xml'), error('don''t be a dick'), end
dom=xmlread(filename);

%% create the actors struct array fields (not the complete list)
rawActors=dom.getElementsByTagName('NxActorDesc');
nActors=rawActors.getLength;
tree.NxuPhysicsCollection.NxSceneDesc.NxActorDesc=struct(...
    'globalPose',cell(1,nActors),...
    'linearVelocity',cell(1,nActors),...
    'mass',cell(1,nActors),...
    'NxBodyDesc',cell(1,nActors),...
    'density',cell(1,nActors),...
    'NxActorFlag',cell(1,nActors),...
    'group',cell(1,nActors),...
    'ATTRIBUTE',cell(1,nActors));

%% fill up the important fields: globalPose, linearVelocity, mass
for k=1:nActors % careful, java objects zero-based
    % globalPose
    sgp=char(rawActors.item(k-1).getElementsByTagName('globalPose').item(0).getTextContent());
    gPose=textscan(sgp,'%f',12);
    tree.NxuPhysicsCollection.NxSceneDesc.NxActorDesc(k).globalPose=...
        gPose{1}';
    
    % linear velocity
    slv=char(rawActors.item(k-1).getElementsByTagName('linearVelocity').item(0).getTextContent());
    lVel=textscan(slv,'%f',3);
    tree.NxuPhysicsCollection.NxSceneDesc.NxActorDesc(k).linearVelocity=...
        lVel{1}';
    
    % mass
    sms=char(rawActors.item(k-1).getElementsByTagName('mass').item(0).getTextContent());
    mass=str2double(sms);
    tree.NxuPhysicsCollection.NxSceneDesc.NxActorDesc(k).mass=...
        mass;
end

%% sniff for shape information
for k=1:nActors % careful
    thisRawActor=rawActors.item(k-1);
    thisRawSphere=thisRawActor.getElementsByTagName('NxSphereShapeDesc');
    if thisRawSphere.getLength > 0
        rad=str2double(char(thisRawSphere.item(0).getAttribute('radius')));
        tree.NxuPhysicsCollection.NxSceneDesc.NxActorDesc(k).NxSphereShapeDesc.ATTRIBUTE.radius=rad;
    end
    thisRawPoly=thisRawActor.getElementsByTagName('NxConvexShapeDesc');
    if thisRawPoly.getLength > 0
        meshLabel=char(thisRawPoly.item(0).getAttribute('meshData'));
        tree.NxuPhysicsCollection.NxSceneDesc.NxActorDesc(k).NxConvexShapeDesc.ATTRIBUTE.meshData=...
            meshLabel;
    end
end

%% the convex mesh data is outside of the actor
rawMeshes=dom.getElementsByTagName('NxConvexMeshDesc');
nMeshes=rawMeshes.getLength;
for k=1:nMeshes
    spts=char(rawMeshes.item(k-1).getElementsByTagName('points').item(0).getTextContent);
    pts=textscan(spts,'%f');
    tree.NxuPhysicsCollection.NxConvexMeshDesc(k).points=pts{1}';
end
end