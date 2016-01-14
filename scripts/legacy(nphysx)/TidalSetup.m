%% A script to set up a Rubble-Pile tidal encounter
% The nPhysX based cuTidal program needs 3-4 input files. The .ini file is the
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
% This setup script will read a saved scene from an xml file, that was
% likely saved from a run of a freely gravitating pile.

clear
clc
close all
physunits off % on for debugging off for performance
si=setUnits;
sol=solsys;
addpath ../scripts

%% PART I - Parameters
% run
runBaseName='test';
backupRate=0;
dumpRate=0;
captureRate=100;
gravSoft=0;
dt=2*si.s;
% rubble
expandby=0.9277;
tPart='poly';
restitution=0.83;
friction=0.5;
rhoBulk=2.0*si.g/si.cm^3; % bulk density (approx!!!)
% spin
saxis=[0 0 -1]; % spin axis (remember orbit plane is xy)
T=inf*si.hr; % spin period
% orbit
skipOrbit=false;
orbitType='hyperbolic';
bigM=sol.Earth.mass;
bigG=si.gravity;
q=1.01*sol.Earth.diameter/2;
vinf=9.0*si.km/si.s; % v at "infinity", for elliptical orbits use e)
e=0.999915; % eccentricity (<1; for hyperbolic orbits use vinf)
Roche=1.51*(bigM/rhoBulk)^(1/3);
r_ini=12*sol.Earth.diameter/2; % assumed pre-peri-
r_end=38*sol.Earth.diameter/2;


%% PART II - Read saved .xml or .mat
% load from file,
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
    end
    tPart='specific';
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
rhoPart=rhoBulk*Vbulk/sum(vVec);
fprintf('using grain density of %g kg/m^3\n',double(rhoPart));
mPart=rhoPart*4/3*pi*rPart^3;
if ~exist('dt','var')
dt=1e-3*sqrt(Vbulk/(2*bigG*nPart*mPart));
end

%% PART III - Velocities
% Each particle must be given a velocity such that the result for the
% complete pile is a the specified spin period. We will be looking at the
% motion in the center-of-mass frame so there is no orbital velocity.
w=2*pi/T*saxis/norm(saxis);
w=repmat(w,nPart,1);
vel=cross(w,pos);

%%  PART IV - Orbit
if ~skipOrbit
% I want to create a section of an elliptical or hyperbolic orbit. I want
% to set t=0 at perihelion, negative pre- and positive post. This would
% give negative angles for pre-peri from the kepler equation, and so the
% orbiting body will be approaching from below the x-axis, usually from the
% 4th quadrant, and swing up and left into the 1st quadrant. Note that in
% PhysX view the positive x direction is to the left, so the orbit will
% start in the 3rd quadrant and swing to the 2nd.

%% For elliptical orbits
if (strcmp(orbitType,'elliptical'))

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

%% For hypoerbolic orbits.
elseif strcmp(orbitType,'hyperbolic')

% Start by deriving orbital elements
if r_ini<q || r_end<q
    error('don''t be a dick')
end
a=bigG*bigM/vinf^2;
e=1+q/a;
% Then get the initial and final eccentric anomalies from the specified
% initial and final distances.
coshF_ini=(1+r_ini/a)/e;
F_ini=-acosh(double(coshF_ini)); % pre-peri- negative angle
coshF_end=(1+r_end/a)/e;
F_end=acosh(double(coshF_end)); % post-per- positive angle
if ~(isreal(F_ini)&&isreal(F_end))
    error('something''s wrong')
end
% Then get the initial and final times corresponding to these angles
% (remember that t=0 is at peri- so the initial time will be negative),
t_ini=sqrt(a^3/(bigG*bigM))*(e*sinh(F_ini)-F_ini);
t_end=sqrt(a^3/(bigG*bigM))*(e*sinh(F_end)-F_end);
% and create a time vector
if (t_end-t_ini)/dt > 1e6, error('looong orbit, please reconsider'); end
tVec=t_ini:dt:t_end;

% Now the hard part - inverting Kepler's equation.
% Shunning premature optimization i will try the naive approach first.
% Numerically converging on each time point independently. To not be
% comepletely obtuse, I will at least use the previously found time point
% as a starting guess...
tauVec=tVec*sqrt(bigG*bigM/a^3); % dimensionless time vector
kepler=@(tau,ef)tau-(e*sinh(ef)-ef);
FVec=ones(size(tauVec));
FVec(1)=double(F_ini); % we know the first point.
progressbar(0,5);
for k=2:numel(tauVec)
    fun=@(ef)double(kepler(tauVec(k),ef));
    guess=FVec(k-1);
    FVec(k)=fzero(fun,guess);
    progressbar(k/numel(tVec));
end

% Finally, translate eccentric anomalies to meaningful space coordinates.
el=q*sqrt(bigG*bigM*(1/a+2/q));
aVec=a*(e*cosh(FVec)-1);
cosfVec=(el^2./(aVec*bigG*bigM)-1)/e;
fVec=acos(cosfVec);
fVec=fVec.*sign(FVec); % make pre-peri angles negative
xVec=aVec.*cos(fVec);
yVec=aVec.*sin(fVec);


else
    error('unknown orbit type');
end
else
    fprintf('skipping orbit calculation\n')
end


%% Part V - Writing the files
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
