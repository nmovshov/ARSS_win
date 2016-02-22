%% ECAFY
close all
clc

%% Get the directory list
dirs=dir('ECA*');
dirs=dirs([dirs.isdir]);
sortedDirNames=sort({dirs.name})';
%sortedDirNames(end+1:end+8)=sortedDirNames(1:8);
%sortedDirNames(1:8)=[];

%% Call kreutzalize for each run and add to the matrix
%progressbar
for k=1:30
filename=dir([sortedDirNames{k} '/dump*.xml']);
pData=kreutzalize([sortedDirNames{k} '/' filename.name]);
xVec=pData(:,1);
yVec=pData(:,2);
zVec=pData(:,3);
mVec=pData(:,7);
rVec=ones(size(xVec))*70;

% % Find clusters
% skinFactor=1.2;
% dset=num2cell([xVec,yVec,zVec,rVec],2); % eclazz expects cell array data
% inClust=@(v1,v2)norm(v1(1:3)-v2(1:3))<(v1(4)+v2(4))*skinFactor; % equivalence R
% clustVec=eclazz(inClust,dset); % NOTE: may take a few minutes
% clusts=unique(clustVec);
% clumpMassVec=nan(size(clusts));
% for j=1:numel(clusts)
%     mClust=mVec(clustVec==clusts(j));
%     clumpMassVec(j)=sum(mClust);
% end
 
 col=1+floor((k-1)/5);
 row=5-mod(k-1,5);

 h= subplot(5,6,sub2ind([6 5],col,row)); lh=plot(xVec/1000,yVec/1000,'k.','markersize',6);
 set(h,'xtick',[],'ytick',[],'xlim',[-200 200],'ylim',[-200 200])

% % save a record of all particles
% save(['forPaul/' sortedDirNames{k} '.dynamics'],'pData','-ascii')
% save(['forPaul/' sortedDirNames{k} '.clumps'],'clumpMassVec','-ascii')

%progressbar(k/16)
end