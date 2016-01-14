%% A script to set up a Rubble-Pile SL9-Jupiter encounter
% The nPhysX based SL9 program needs 3-4 input files. The .ini file is the
% general setup file used in all nPhysX simulations. The .ic file specifies
% the dynamical initial conditions - a list of positions to generate actors
% in and a matching list of starting velocities. The positions are used to
% create a rubble-pile with specific shape and compactness. The velocities
% are used to give the pile initial spin. A pre-calculated orbit is defined
% in the .orb file. An optional .sh file specifies individual particle
% shapes. If this file exists the number of shapes must match the number of
% actors derived from .ic, and the shapes are created in the corresponding
% position from .ic. If this file is missing generic shapes will be created
% at the .ic positions, based on the rubble_type and grain_size parameters
% in .ini. Unfortunately, there is some overlap in the information
% contained in these files. For example, the integration time step is
% specified in .ini, but implied in .orb. The number of actors is specified
% in .ini but implied in .ic. The tiding mass is specified in .ini but
% implied in .ic and .orb. The grain size is specified in .ini but implied
% in .sh. To make sure these parameters match in a given simulation, it is
% best to create all 3 files in a single script - this one.
%
% To complicate matters more, there are two distinct setup scenarios. One
% is to setup a scene from scratch and the other is to use a saved scene
% from a .xml file (or the equivalent tree structure in .mat). There are
% many overlapping calculations but many contradictory ones, so I decided
% to make two separate branches of this script for these two options.

clear
clc
close all
physunits off % on for debugging off for performance
si=setUnits;
sol=solsys;
addpath ../scripts

%% Parameters
% run
runBaseName='test';
freshSetup=false;
backupRate=0;
dumpRate=0;
captureRate=0;
gravSoft=0;
% rubble
expandby=1;
tPart='poly';
restitution=0.83;
friction=0.5;
rhoBulk=0.3*si.g/si.cm^3; % bulk density (approx!!!)
% spin
saxis=[0 0 -1]; % spin axis (remember orbit plane is xy)
T=inf*si.hr; % spin period
% orbit
skipOrbit=false;
bigM=sol.Jupiter.mass;
bigG=si.gravity;
q=1.31*sol.Jupiter.diameter/2; % perijove
e=0.997; % eccentricity (<1; for hyperbolic orbits see teic.m)
Roche=1.51*(bigM/rhoBulk)^(1/3);
r_ini=4*Roche; % assumed pre-peri-
r_end=7.5*sol.Jupiter.diameter; % assumed post-peri-

%% ***Fresh Branch***
if freshSetup

%% more rubble parameters
nPart=1024; %#ok<UNRCH> % desired number of "particles"
tPart='spheres'; % type of particle
a1=1.5*si.km; a2=1.5*si.km; a3=1.5*si.km; % semi- axes lengths (x,y,z)
rPart=(0.5*a1*a2*a3/nPart)^(1/3); % approx size of one particle
rhoPart=rhoBulk/0.5; % material density
mPart=rhoPart*4/3*pi*rPart^3;

%% Time scale
% Based on the number and size of desired particles or on the shape of the
% desired ellipsoid, a dynamical time scale is calculated.
if ~exist('dt','var')
dt=1e-3*sqrt(a1^3/(2*bigG*nPart*mPart));
end

%% Positions
% The goal here is to create a rubble pile with a specified numebr of
% elements of a given size, that results in something that looks like a
% triaxial ellipsoid. The shape of the elements is not specified, so there
% is some uncertainty and the "packing" algorithm needs to allow for this.
% In practice this means that the generated pile is initially quite fluffy
% and may need to be relaxed in the actual simulation.

% ok, start by calculating the necessary spacing between particles
spc=((4/3*pi*a1*a2*a3)/nPart)^(1/3);
if spc<2*rPart, error('ouch - too tight!'); end

% make an array to hold the position vectors and start the filling process
% at the right tip
pos=zeros(nPart,3)*si.m;
kp=1;
x=-a1;

% now keep filling between the lines
while x<=a1
    ymax=sqrt(round((1-x^2/a1^2)*a2^2));
    y=-ymax;
    while y<=ymax
        zmax=sqrt(round((1-x^2/a1^2-y^2/a2^2)*a3^2));
        z=-zmax;
        while z<=zmax
            pos(kp,:)=[x y z];
            kp=kp+1;
            if kp>nPart, break, end
            z=z+spc;
        end
        if kp>nPart, break, end
        y=y+spc;
    end
    if kp>nPart, break, end
    x=x+spc;
end

%% ***end of Fresh Branch***
else
%% ***saved .xml or .mat branch***
% load from file,
if strcmp(tPart,'spheres') && exist('stdSpheres.mat','file')
    load stdSpheres
elseif strcmp(tPart,'poly') && exist('stdPoly.mat','file')
    load stdPoly
else
[filename pathname]=uigetfile({'*.mat';'*.xml'});
filename=fullfile(pathname,filename);
if (filename)
    [~,~,xt]=fileparts(filename);
    if strcmp(xt,'.xml')
        tree=readNXU(filename); % custom xml reader, may break if nxustream changes
    elseif strcmp(xt,'.mat')
        load(filename); % or maybe you have the tree in a mat file already
    else
        error('don''t be a dick');
    end
end
end
% and verify DOM structure exists
if ~exist('tree','var') || ~isfield(tree,'NxuPhysicsCollection')
    error('don''t be a dick');
end
actors=tree.NxuPhysicsCollection.NxSceneDesc.NxActorDesc;

%% Get x,y,z etc variables in usable form
globalPose=cell2mat({actors.globalPose}');
xVec=globalPose(:,10)*si.m*expandby;
yVec=globalPose(:,11)*si.m*expandby;
zVec=globalPose(:,12)*si.m*expandby;

% make an array to hold the position vectors
pos=[xVec yVec zVec];
nPart=size(pos,1);
 
%% Sniff for shape information
if isfield(actors,'NxSphereShapeDesc') % they are most likely spheres
    if ~strcmp(tPart,'spheres')
    warning('sphere shapes detected - tPart reset to spheres')
    tPart='spheres';
    end
    dummy=[actors.NxSphereShapeDesc];
    dummy=[dummy.ATTRIBUTE];
    rVec=[dummy.radius]*expandby;
    vVec=4*pi/3*rVec.^3;
    if numel(unique(rVec))==1
        fprintf('uniform spheres detected - using single radius primitives\n')
        rPart=rVec(1)*si.m;
    else
        rPart=mean(rVec);
    end
elseif isfield(actors,'NxConvexShapeDesc') % they are most likely polys
    if ~strcmp(tPart,'poly')
        warning('poly shapes detected - tPart reset to specific')
        tPart='specific';
    end
    dummy=[actors.NxConvexShapeDesc];
    dummy=[dummy.ATTRIBUTE];
    meshLabelVec={dummy.meshData}'; % label of mesh used by actors
    % get convex shapes in usable form: this is quite brittle due to the
    % unpredictable format of the nxustream2 xmls
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
        shapes{k}=expandby*shapes{k}'; % and make it a nice column major thing
        rVec(k)=mean(diag(shapes{k}*shapes{k}'));
        [~,vVec(k)]=convhulln(shapes{k});
    end
    rVec=sqrt(rVec);
    rPart=mean(rVec);
end

%% Using shape information to adjust density and time step
% I use an iterative two-pass convex hull method to find the volume
% enclosed by the rubble pile. If convergance to 1% is not achieved by 32
% iterations the last one is used and a warning is thrown.
[~,Vold]=convhulln(pos);
for npts=12:32
    Vnew=BulkVolume(xVec,yVec,zVec,rVec,npts);
    if abs(Vnew-Vold)/Vold < 0.01, break, end
    Vold=Vnew;
end
if npts==32
    warning('two-pass convhull method did not converge by 32 steps, using last found volume')
end
Vbulk=Vnew*si.m^3;
fprintf('detected a pile with radius %g m\n',double((3*Vbulk/4/pi)^(1/3)));
rhoPart=rhoBulk*Vbulk/sum(vVec);%TODO FIX THIS
fprintf('using grain density of %g kg/m^3\n',double(rhoPart));
mPart=rhoPart*4/3*pi*rPart^3;
if ~exist('dt','var')
dt=1e-3*sqrt(Vbulk/(2*bigG*nPart*mPart));
end
%% ***end saved .mat or .xml branch***
end

%% Velocities
% Each particle must be given a velocity such that the result for the
% complete pile is a the specified spin period. We will be looking at the
% motion in the center-of-mass frame so there is no orbital velocity.
w=2*pi/T*saxis/norm(saxis);
w=repmat(w,nPart,1);
vel=cross(w,pos);

%% Orbit
% I want to create a section of an elliptical orbit. I want to set t=0 at
% perihelion, negative pre- and positive post. This would give negative
% angles for pre-peri from the kepler equation, and so the orbiting body
% will be approaching from below the x-axis, usually from the 4th quadrant,
% and swing up and left into the 1st quadrant. Note that in PhysX view the
% positive x direction is to the left, so the orbit will start in the 3rd
% quadrant and swing to the 2nd.
if ~skipOrbit
% First figure out the time span (t=0 at peri-). Start by calculating the
% derived orbital elements.
a=q/(1-e); % semi-major axis
Q=a*(1+e); % ap- (just for input checking)
if r_ini<q || r_end<q || r_ini>Q || r_end>Q
    error('don''t be a dick')
end
P=2*pi*sqrt(a^3/(bigG*bigM)); % orbital period
% Then get the initial and final eccentric anomalies from the specified
% initial and final distances.
cosE_ini=(1-r_ini/a)/e;
E_ini=-acos(cosE_ini); % pre-peri- negative angle
cosE_end=(1-r_end/a)/e;
E_end=acos(cosE_end); % post-per- positive angle
if ~(isreal(E_ini)&&isreal(E_end))
    error('something''s wrong')
end
% Then get the initial and final times corresponding to these angles
% (remember that t=0 is at peri- so the initial time will be negative),
t_ini=P/(2*pi)*(E_ini-e*sin(E_ini));
t_end=P/(2*pi)*(E_end-e*sin(E_end));
% and create a time vector
if (t_end-t_ini)/dt > 1e6, error('looong orbit, please reconsider'); end
tVec=t_ini:dt:t_end;

% Now the hard part - inverting Kepler's equation.
% Shunning premature optimization i will try the naive approach first.
% Numerically converging on each time point independently. To not be
% comepletely obtuse, I will at least use the previously found time point
% as a starting guess...
kepler=@(tau,psi)tau-P/(2*pi)*(psi-e*sin(psi)); % Kepler's equation as fun
EVec=ones(size(tVec));
EVec(1)=double(E_ini); % we know the first point.
progressbar(0,5);
for k=2:numel(tVec)
    fun=@(psi)double(kepler(tVec(k),psi));
    guess=EVec(k-1);
    EVec(k)=fzero(fun,guess);
    progressbar(k/numel(tVec));
end

% Finally, translate eccentric anomalies to meaningful space coordinates.
fVec=acos((cos(EVec)-e)./(1-e*cos(EVec))); % orbit phase angle
fVec=fVec.*sign(EVec); % make pre-peri angles negative
aVec=a*(1-e^2)./(1+e*cos(fVec)); % orbit equation
xVec=aVec.*cos(fVec);
yVec=aVec.*sin(fVec);

else
    fprintf('skipping orbit calculation\n')
end


%% Part IV - Writing the files
%% .ic file
% Initial conditions files for nPhysX programs can have a header that
% must end in ---END HEADER--- and is followed by a rectangular 6-column
% matrix of arbitrary size.
fprintf('writing .ic file\n')
% write the header
fid=fopen([runBaseName '.ic'],'wt');
fprintf(fid,'%s\n%s\n%s\n%s\n%s\n',...
    '---BEGIN HEADER---',...
    'Initial conditions file for sgpile experiments.',...
    'Columns are (for center of mass, mks is implied):',...
    '[x] [y] [z] [v_x] [v_y] [v_z]',...
    '---END HEADER---');
fclose(fid);

% write the data
ic=[double(pos) double(vel)];
save([runBaseName '.ic'],'ic','-ascii','-double','-append');

%% .sh file
% shape file for nPhysX programs can have a header that must end in ---END
% HEADER--- followed by shape information for each actor.
if exist('shapes','var') || numel(unique(rVec))>1
fprintf('writing .sh file\n')
% write the header
fid=fopen([runBaseName '.sh'],'wt');
fprintf(fid,'%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n',...
    '---BEGIN HEADER---',...
    'This is the shapes file for sgpile. The shape specification for *all* actors are in',...
    'this one file, so order matters. The sections separated by *SHAPE* are the vertex',...
    'positions belonging to the actor whose center of mass is on the corresponding line',...
    'of the initial conditions file. A section that starts with *SPHERE* must have just',...
    'one number, to be interpreted as a radius.',...
    '---END HEADER---');

% write the data
if strcmp(tPart,'specific') || strcmp(tPart,'poly')
    for k=1:numel(shapes)
        fprintf(fid,'SHAPE\n');
        fprintf(fid,'%f %f %f\n',shapes{k}');
    end
elseif strcmp(tPart,'spheres')
    for k=1:numel(rVec)
        fprintf(fid,'SPHERE\n');
        fprintf(fid,'%f\n',rVec(k));
    end
end
fclose(fid);
end

%% .ini file
% The general ini file for the program. most parameters are static but it's
% a pain to modify the few dynamic parameters by hand.
fprintf('writing .ini file\n')
fid=fopen([runBaseName '.ini'],'wt');
fprintf(fid,'%s\n','; This is the configuration file for program SL9. It''s fairly robust but');
fprintf(fid,'%s\n','; don''t push it too much. Adding new keys and sections in any order is fine.');
fprintf(fid,'%s\n','; Try to avoid spurious white spaces in keys and values. The weakest point');
fprintf(fid,'%s\n','; is that numerical values are read as strings and converted using atof(),');
fprintf(fid,'%s\n','; which is hard to verify as invalid arguments produce implementation ');
fprintf(fid,'%s\n','; dependent results.');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','[program]');
fprintf(fid,'%s\n',['base_name=' runBaseName]);
fprintf(fid,'%s\n',['dump_output_frequency=' num2str(dumpRate)]);
fprintf(fid,'%s\n',['dump_backup_frequency=' num2str(backupRate)]);
fprintf(fid,'%s\n','hud_refresh_rate=1');
fprintf(fid,'%s\n','start_at_frame=0');
fprintf(fid,'%s\n','quit_when_finished=true');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','[PhysX]');
fprintf(fid,'%s\n',['skin_width=' num2str(0.001*rPart)]);
fprintf(fid,'%s\n','adaptive_force=true');
fprintf(fid,'%s\n','cone_friction=true');
fprintf(fid,'%s\n','linear_damping=0.00');
fprintf(fid,'%s\n','angular_damping=0.00');
fprintf(fid,'%s\n',['bounce_eps=' num2str(1/3*sqrt(8*bigG*rhoPart*rPart^2))]);
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','[simulation]');
fprintf(fid,'%s\n',['integrator_dt=' num2str(dt)]);
fprintf(fid,'%s\n','real_time=false');
fprintf(fid,'%s\n','timing_method=variable');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','[physical]');
fprintf(fid,'%s\n',['big_G=' num2str(bigG)]);
fprintf(fid,'%s\n','gravity_type=all_pairs');
fprintf(fid,'%s\n',['gravity_softening_factor=' num2str(gravSoft)]);
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','[material]');
fprintf(fid,'%s\n',['restitution=' num2str(restitution)]);
fprintf(fid,'%s\n',['dynamicFriction=' num2str(friction)]);
fprintf(fid,'%s\n',['staticFriction=' num2str(friction)]);
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','[experiment]');
fprintf(fid,'%s\n','ts_method=const');
fprintf(fid,'%s\n','experiment_type=tidal_encounter');
fprintf(fid,'%s\n',['tiding_mass=' num2str(bigM)]);
fprintf(fid,'%s\n',['rubble_type=' tPart]);
fprintf(fid,'%s\n','uniform_rubble=true');
fprintf(fid,'%s\n',['grain_size=' num2str(rPart*2)]);
fprintf(fid,'%s\n',['grain_density=' num2str(rhoPart)]);
fprintf(fid,'%s\n','control_mode=auto');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','[DevIL]');
fprintf(fid,'%s\n',['frame_capture_frequency=' num2str(captureRate)]);
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','[GLUT]');
fprintf(fid,'%s\n','tracking_camera=true');
fprintf(fid,'%s\n',['CamZbuf.near=' num2str(0.1*rPart)]);
fprintf(fid,'%s\n',['CamZbuf.far=' num2str(800*rPart)]);
fprintf(fid,'%s\n','CamSpeed=1');
fclose(fid);

%% .orb file
if ~skipOrbit
fprintf('writing .orb file\n')
fid=fopen([runBaseName '.orb'],'wt');
fprintf(fid,'%s\n%s\n%s\n%s\n',...
    '---BEGIN HEADER---',...
    'Orbital data in cartesian coordinates (most likely mks).',...
    'Columns are:',...
    '[t] [x] [y]',...
    '---END HEADER---');
fclose(fid);
orbit=double([tVec' xVec' yVec']);
save([runBaseName '.orb'],'orbit','-ascii','-double','-append')
end

%% .mat file (a record of parameters)
fprintf('writing .mat file (record of parameters)\n')
save(runBaseName);
