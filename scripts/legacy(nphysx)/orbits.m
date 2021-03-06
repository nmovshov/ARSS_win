%% A script to generate elliptical orbits with specified q e and dt
% I want to create a section of an elliptical orbit. I want to set t=0 at
% perihelion, negative pre- and positive post. This would give negative
% angles for pre-peri from the kepler equation, and so the orbiting body
% will be approaching from below the x-axis, usually from the 4th quadrant,
% and swing up and left into the 1st quadrant. Note that in PhysX view the
% positive x direction is to the left, so the orbit will start in the 3rd
% quadrant and swing to the 2nd.
clear
clc
close all
physunits on % on for debugging off for performance
si=setUnits;
sol=solsys;

%% Parameters
name='test';
q=0.98*si.astronomical_unit;
e=0.0167;
dt=1*si.day;
r_ini=1.013*si.astronomical_unit; % assumed pre-peri-
r_end=1.013*si.astronomical_unit;
bigM=sol.Sun.mass;
bigG=si.gravity;

%% The time span (t=0 at peri-)
% start by calculating the derived orbital elements
a=q/(1-e); % semi-major axis
Q=a*(1+e); % ap- (just for input checking)
if r_ini<q || r_end<q || r_ini>Q || r_end>Q
    error('don''t be a dick')
end
P=2*pi*sqrt(a^3/(bigG*bigM)); % orbital period
% then get the initial and final eccentric anomalies from the specified
% initial and final distances
cosE_ini=(1-r_ini/a)/e;
E_ini=-acos(cosE_ini); % pre-peri- negative angle
cosE_end=(1-r_end/a)/e;
E_end=acos(cosE_end); % post-per- positive angle
if ~(isreal(E_ini)&&isreal(E_end))
    error('something''s wrong')
end
% then get the initial and final times corresponding to these angles
% (remember that t=0 is at peri- so the initial time will be negative)
t_ini=P/(2*pi)*(E_ini-e*sin(E_ini));
t_end=P/(2*pi)*(E_end-e*sin(E_end));
% and create a time vector
tVec=t_ini:dt:t_end;

%% The hard part - inverting Kepler's equation
% Shunning premature optimization i will try the naive approach first.
% Numerically converging on each time point independently. To not be
% comepletely obtuse, I will at least use the previously found time point
% as a starting guess...
kepler=@(tau,psi)tau-P/(2*pi)*(psi-e*sin(psi)); % Kepler's equation as fun
EVec=ones(size(tVec));
EVec(1)=double(E_ini); % we know the first point.
for k=2:numel(tVec)
    fun=@(psi)double(kepler(tVec(k),psi));
    guess=EVec(k-1);
    EVec(k)=fzero(fun,guess);
end

%% Translate eccentric anomalies to meaningful space coordinates
fVec=acos((cos(EVec)-e)./(1-e*cos(EVec))); % orbit phase angle
fVec=fVec.*sign(EVec); % make pre-peri angles negative
rVec=a*(1-e^2)./(1+e*cos(fVec)); % orbit equation
xVec=rVec.*cos(fVec);
yVec=rVec.*sin(fVec);

%% Save orbit data to file
% a .orb file with header will be created (and overrun) in the current
% directory
fid=fopen([name '.orb'],'wt');
fprintf(fid,'%s\n%s\n%s\n%s\n',...
    '---BEGIN HEADER---',...
    'Orbital data in cartesian coordinates (most likely mks).',...
    'Columns are:',...
    '[t] [x] [y]',...
    '---END HEADER---');
fclose(fid);
orbit=double([tVec' xVec' yVec']);
save([name '.orb'],'orbit','-ascii','-append')