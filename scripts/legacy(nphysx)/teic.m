%% A script to generate initial conditions for a hyperbolic tidal encounter
clear
clc
close all
physunits off % on for debugging off for performance
si=setUnits;
sol=solsys;

%% Parameters
% run
name='test';
% shape
nPart=147; % desired number of "particles"
rPart=245*si.m; % size of one particle, side of a bounding box.
rhoPart=3.6*si.g/si.cm^3;
mPart=rhoPart*4/3*pi*rPart^3;
a1=2.8*si.km; a2=1.7*si.km; a3=1.5*si.km; % semi- axes lengths (x,y,z)
% spin
saxis=[0 0 -1]; % spin axis (remember orbit plane is xy)
T=6*si.hr; % spin period
% orbit
Roche=3.4*sol.Earth.diameter/2;
r_ini=1*Roche;
q=1.1*sol.Earth.diameter/2; % periapse distance
vinf=0*si.km/si.s; % encounter velocity, zero for parabolic
% tider
G=si.gravity;
bigM=sol.Earth.mass;

%% Positions
% The goal here is to create a rubble pile with a specified numebr of
% elements of a given size, that results in something that looks like a
% triaxial ellipsoid. The shapes of the elements is not specified, so there
% is some uncertainty and the "packing" algorithm needs to allow for this.
% In practice this means that the generated pile is initially quite fluffy
% and needs to be relaxed in the actual simulation.

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

%% Velocities
% Each particle must be given a velocity. The desired result for the
% complete pile is an orbit of the center of mass with some specified
% periapse and specified velocity at infinity.

% start by calculating the CMvel and CMpos vectors that correspond to the
% desired orbital parameters. Given a vinf magnitude, the orbital energy
% (per unit mass) is
E=1/2*vinf^2;
% Now given the desired perihelion q we can solve for l, the required angular
% momentum (per unit mass).
l=sqrt(2*E+2*G*bigM/q)*q;
% We can also solve for the eccentricity in terms of vinf (look at (e-1)/q):
e=1+(q*vinf^2)/(G*bigM);
% With given energy and angular momentum (per unit mass) we can find the
% components of the velocity and position vectors for any radial distance r:
costeta=(l^2/(G*bigM*r_ini)-1)/e;
sinteta=sqrt(1-costeta^2);
CMpos=[r_ini*costeta r_ini*sinteta];
Rot=[costeta -sinteta;sinteta costeta];
CMvel=[-sqrt(2*E+2*G*bigM/r_ini-l^2/r_ini^2);-l/r_ini]; % the signs matter!
CMvel=(Rot*CMvel)';
CMvel=[CMvel 0*si.m/si.s];
vel=repmat(CMvel,nPart,1);

% Now, besides the orbital velocity, we add spin velocity:
w=2*pi/T*saxis/norm(saxis);
w=repmat(w,nPart,1);
vel=vel+cross(w,pos);

%% Time step
% I choose the minimum of the dynamical time scale and orbital time scale.
dtDyn=1e-3*sqrt(a1^3/(2*G*nPart*mPart));
dtOrb=1e-4*r_ini/sqrt(CMvel*CMvel');
dt=min([dtDyn dtOrb]);

%% Scaling
% We need to choose appropriate scale units. Remember, these values will be
% used in a single precision code.
% LS=1*si.km; % length unit
% MS=1e10*si.kg; % mass unit
% TS=10*si.s; % time unit
% pos=pos/LS;
% vel=vel/(LS/TS);
% G=G/(LS^3*MS^-1*TS^-2) % make sure to use in sgpile.ini
% r=rPart/LS
% rho=3500*si.kg/si.m^3/(MS/LS^3) % make sure to use in sgpile.ini
% M=bigM/MS

%% Output (ic file)
% Initial conditions files for the sgpile program can have a header that
% must end in ---END HEADER--- and is followed by a rectangular 6-column
% matrix of arbitrary size.

% write the header
fid=fopen(['../sgpile/' name '/' name '.ic'],'wt');
fprintf(fid,'%s\n%s\n%s\n%s\n%s\n',...
    '---BEGIN HEADER---',...
    'Initial conditions file for sgpile experiments.',...
    'Columns are (for center of mass, mks is implied):',...
    '[x] [y] [z] [v_x] [v_y] [v_z]',...
    '---END HEADER---');
fclose(fid);

% write the data
ic=[double(pos) double(vel)];
save(['../sgpile/' name '/' name '.ic'],'ic','-ascii','-append');

%% Output (ini file)
% the general ini file for the program. most parameters are static but it's
% a pain to modify the few dynamic parameters by hand.
fid=fopen(['../sgpile/' name '/' name '.ini'],'wt');
fprintf(fid,'%s\n','; This is the configuration file for program SGPILE. It''s fairly robust but');
fprintf(fid,'%s\n','; don''t push it too much. Adding new keys and sections in any order is fine.');
fprintf(fid,'%s\n','; Try to avoid spurious white spaces in keys and values. The weakest point');
fprintf(fid,'%s\n','; is that numerical values are read as strings and converted using atof(),');
fprintf(fid,'%s\n','; which is hard to verify as invalid arguments produce implementation ');
fprintf(fid,'%s\n','; dependent results.');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','[program]');
fprintf(fid,'%s\n',['base_name=' name]);
fprintf(fid,'%s\n','dump_output_frequency=0');
fprintf(fid,'%s\n','dump_backup_frequency=0');
fprintf(fid,'%s\n','hud_refresh_rate=1');
fprintf(fid,'%s\n','start_at_frame=0');
fprintf(fid,'%s\n','quit_when_finished=false');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','[PhysX]');
fprintf(fid,'%s\n',['skin_width=' num2str(0.001*rPart)]);
fprintf(fid,'%s\n','adaptive_force=true');
fprintf(fid,'%s\n','cone_friction=true');
fprintf(fid,'%s\n','linear_damping=0.00');
fprintf(fid,'%s\n','angular_damping=0.00');
fprintf(fid,'%s\n',['bounce_eps=' num2str(1/3*sqrt(8*G*rhoPart*rPart^2))]);
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','[simulation]');
fprintf(fid,'%s\n',['integrator_dt=' num2str(dt)]);
fprintf(fid,'%s\n','real_time=false');
fprintf(fid,'%s\n','timing_method=variable');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','[physical]');
fprintf(fid,'%s\n',['big_G=' num2str(G)]);
fprintf(fid,'%s\n','gravity_type=all_pairs');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','[material]');
fprintf(fid,'%s\n','restitution=0.8');
fprintf(fid,'%s\n','dynamicFriction=0.0');
fprintf(fid,'%s\n','staticFriction=0.0');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','[experiment]');
fprintf(fid,'%s\n','ts_method=const');
fprintf(fid,'%s\n','experiment_type=tidal_encounter');
fprintf(fid,'%s\n',['tiding_mass=' num2str(bigM)]);
fprintf(fid,'%s\n',['tiding_center.x=' num2str(-CMpos(1))]);
fprintf(fid,'%s\n',['tiding_center.y=' num2str(-CMpos(2))]);
fprintf(fid,'%s\n','rubble_type=spheres');
fprintf(fid,'%s\n','uniform_rubble=true');
fprintf(fid,'%s\n',['grain_size=' num2str(rPart*2)]);
fprintf(fid,'%s\n',['grain_density=' num2str(rhoPart)]);
fprintf(fid,'%s\n','control_mode=auto');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','[DevIL]');
fprintf(fid,'%s\n','frame_capture_frequency=0');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','[GLUT]');
fprintf(fid,'%s\n','tracking_camera=false');
fprintf(fid,'%s\n',['CamZbuf.near=' num2str(0.1*rPart)]);
fprintf(fid,'%s\n',['CamZbuf.far=' num2str(800*rPart)]);
fprintf(fid,'%s\n','CamSpeed=1');
fclose(fid);