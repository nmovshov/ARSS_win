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
