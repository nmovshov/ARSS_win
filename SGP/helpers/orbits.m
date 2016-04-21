%% A script to generate elliptical or hyperbolic orbits
% I want to create a (section of an) orbit. I want to set t=0 at perihelion,
% negative pre- and positive post. This would give negative angles for pre-peri
% from the kepler equation, and so the orbiting body will be approaching from
% below the x-axis, usually from the 4th quadrant, and swing up and left into
% the 1st quadrant.
clear
clc
close all

%% Parameters in code units
orbit_name = 'test';
orbit_type = 'bound'; % ['bound'|'hyperbolic']
bigM = 2e6;
bigG = 6.674e-5;
q = 20;
e = 0.4; % bound only
vinf = 0.1; % hyperbolic only
r_ini = 1; % for bound this is a fraction of Q (for full orbit use 1) and ...
r_end = 0.9; % ... for hyperbolic this is a multiple of q
dt = 0.02;

switch orbit_type
    case 'bound'
        % Start by calculating the derived orbital elements
        a = q/(1 - e); % semi-major axis
        Q = a*(1 + e); % ap- (just for input checking)
        r_ini = r_ini*Q;
        r_end = r_end*Q;
        if (r_ini < q) || (r_end < q) || (r_ini > Q) || (r_end > Q)
            error('don''t be a dick')
        end
        P = 2*pi*sqrt(a^3/(bigG*bigM)); % orbital period
        
        % Calculate initial and final eccentric anomalies
        cosE_ini = (1 - r_ini/a)/e;
        E_ini = -acos(cosE_ini); % pre-peri- negative angle
        if e == 0, E_ini = -pi*r_ini/Q; end
        cosE_end=(1-r_end/a)/e;
        E_end=acos(cosE_end); % post-per- positive angle
        if e == 0, E_end = pi*r_ini/Q; end
        if ~(isreal(E_ini) && isreal(E_end))
            error('something''s wrong')
        end
        
        % Calculate initial and final times corresponding to these "angles"
        t_ini = P/(2*pi)*(E_ini - e*sin(E_ini));
        t_end = P/(2*pi)*(E_end - e*sin(E_end));
        
        % And create a time vector
        tVec = t_ini:dt:t_end;
        
        % Invert Kepler's equation
        % Shunning premature optimization I will try the naive approach first,
        % numerically converging on each time point independently. To not be
        % comepletely dense I will at least use the previously found time point
        % as a starting guess.
        kepler = @(tau,psi)tau-P/(2*pi)*(psi-e*sin(psi)); % Kepler's equation
        EVec = ones(size(tVec));
        EVec(1) = double(E_ini); % We know the first point
        for k=2:numel(tVec)
            fun = @(psi)double(kepler(tVec(k),psi));
            guess = EVec(k-1);
            EVec(k)=fzero(fun,guess);
        end
        
        % Translate eccentric anomalies back to meaningful space coordinates
        fVec = acos((cos(EVec) - e)./(1 - e*cos(EVec))); % orbit phase angle
        fVec = fVec.*sign(EVec); % make pre-peri angles negative
        rVec = a*(1 - e^2)./(1 + e*cos(fVec)); % orbit equation
        xVec = rVec.*cos(fVec);
        yVec = rVec.*sin(fVec);
    case 'hyperbolic'
        disp 'coming soon'
    case default
        error('Unknown orbit type')
end

%% Save orbit to file
% File with header will be created (and overrun) in the current directory
fid = fopen([orbit_name '.orb'],'wt');
fprintf(fid,'# Orbital (Cartesian) coordinates in code units.\n');
fprintf(fid,'# Columns are:\n');
fprintf(fid,'# [t]  [x]  [y]\n\n');
fclose(fid);
orbit = double([tVec' xVec' yVec']);
save([orbit_name '.orb'],'orbit','-ascii','-append')
