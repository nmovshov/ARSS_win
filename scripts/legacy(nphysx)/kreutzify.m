%% kreutzify
clear all
close all
clc

%% Get the directory list
dirs=dir('kreutz*');
dirs=dirs([dirs.isdir]);
sortedDirNames=sort({dirs.name})';
sortedDirNames(end+1:end+8)=sortedDirNames(1:8);
sortedDirNames(1:8)=[];

%% Create a matrix of figures


%% Call kreutzalize for each run and add to the matrix

for k=1:56
filename=dir([sortedDirNames{k} '/dump*.xml']);
pData=kreutzalize([sortedDirNames{k} '/' filename.name]);
xVec=pData(:,1);
yVec=pData(:,2);
zVec=pData(:,3);
vxVec=pData(:,4);
vyVec=pData(:,5);
vzVec=pData(:,6);
mVec=pData(:,7);
rVec=ones(size(xVec))*70;

%% Find clusters
skinFactor=1.2;
dset=num2cell([xVec,yVec,zVec,rVec],2); % eclazz expects cell array data
inClust=@(v1,v2)norm(v1(1:3)-v2(1:3))<(v1(4)+v2(4))*skinFactor; % equivalence R
clustVec=eclazz(inClust,dset); % NOTE: may take a few minutes
clusts=unique(clustVec);
hClust=nan(length(clusts),8);
for j=1:numel(clusts)
    elmX=xVec(clustVec==clusts(j));
    elmY=yVec(clustVec==clusts(j));
    elmZ=zVec(clustVec==clusts(j));
    elmR=rVec(clustVec==clusts(j));
    elmM=mVec(clustVec==clusts(j));
    elmVX=vxVec(clustVec==clusts(j));
    elmVY=vyVec(clustVec==clusts(j));
    elmVZ=vzVec(clustVec==clusts(j));
    elmV=elmR.^3;
    Vtot=sum(elmV);
    Mtot=sum(elmM);
    if (numel(elmX))>0
        cmX=dot(elmX,elmM)/Mtot;
        cmY=dot(elmY,elmM)/Mtot;
        cmZ=dot(elmZ,elmM)/Mtot;
        cmVX=dot(elmVX,elmM)/Mtot;
        cmVY=dot(elmVY,elmM)/Mtot;
        cmVZ=dot(elmVZ,elmM)/Mtot;
        cmR=max([range(elmX) range(elmY) range(elmZ)])/2;
        hClust(j,:)=[cmX cmY cmZ cmVX cmVY cmVZ Mtot numel(elmX)];
    end
end
hClust(any(isnan(hClust),2),:)=[];
%%
col=1+mod((k-1),8);
row=7-floor((k-1)/8);
h=subplot(7,8,sub2ind([8 7],col,row));
%lh=plot(xVec/1000,yVec/1000,'k.','markersize',1);
lh=plot(hClust(:,1)/1000,hClust(:,2)/1000,'k.','markersize',6);
set(h,'xtick',[],'ytick',[],'xlim',[-100 100],'ylim',[-100 100])
title(num2str(size(hClust,1)));

% % save a record of all particles
% save(['forPaul/' sortedDirNames{k} '.dynamics'],'pData','-ascii')
% save(['forPaul/' sortedDirNames{k} '.clumps'],'clumpMassVec','-ascii')

end