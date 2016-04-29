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
orbit_name = '../sgp';
orbit_type = 'bound'; % ['bound'|'hyperbolic']
bigM = 1.9e21;
bigG = 6.674e-5;
q = 9e5;
e = 0.997; % bound only
vinf = 0.1; % hyperbolic only
r_ini = 0.01; % for bound this is a fraction of Q (for full orbit use 1) and ...
r_end = 0.01; % ... for hyperbolic this is a multiple of q
dt = 0.02; % does NOT have to match ARSS dt

%% Orbit calculation
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
            warning('throwing out impaginary noise')
            E_ini = real(E_ini);
            E_end = real(E_end);
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
        progressbar(0);
        for k=2:numel(tVec)
            fun = @(psi)double(kepler(tVec(k),psi));
            guess = EVec(k-1);
            EVec(k)=fzero(fun,guess);
            progressbar(k/numel(tVec));
        end
        
        % Translate eccentric anomalies back to meaningful space coordinates
        fVec = acos((cos(EVec) - e)./(1 - e*cos(EVec))); % orbit phase angle
        fVec = fVec.*sign(EVec); % make pre-peri angles negative
        rVec = a*(1 - e^2)./(1 + e*cos(fVec)); % orbit equation
        xVec = rVec.*cos(fVec);
        yVec = rVec.*sin(fVec);
    case 'hyperbolic'
        % Start by calculating derived orbital elements
        r_ini = r_ini*q;
        r_end = r_end*q;
        if r_ini < q || r_end < q
            error('don''t be a dick')
        end
        a = bigG*bigM/vinf^2;
        e = 1 + q/a;
        
        % Calculate initial and final eccentric anomalies
        coshF_ini=(1 + r_ini/a)/e;
        F_ini = -acosh(double(coshF_ini)); % pre-peri- negative angle
        coshF_end = (1 + r_end/a)/e;
        F_end = acosh(double(coshF_end)); % post-per- positive angle
        if ~(isreal(F_ini)&&isreal(F_end))
            error('something''s wrong')
        end
        
        % Calculate initial and final times corresponding to these "angles"
        t_ini = sqrt(a^3/(bigG*bigM))*(e*sinh(F_ini) - F_ini);
        t_end = sqrt(a^3/(bigG*bigM))*(e*sinh(F_end) - F_end);
        
        % Create a time vector
        if (t_end - t_ini)/dt > 1e5, error('Long orbit: please reconsider'); end
        tVec = t_ini:dt:t_end;
        
        % Invert Kepler's equation
        % Shunning premature optimization i will try the naive approach first,
        % numerically converging on each time point independently. To not be
        % comepletely dense, I will at least use the previously found time point
        % as a starting guess.
        tauVec = tVec*sqrt(bigG*bigM/a^3); % dimensionless time vector
        kepler = @(tau, ef)tau - (e*sinh(ef) - ef); % Kepler's equation
        FVec = ones(size(tauVec));
        FVec(1) = double(F_ini); % we know the first point.
        progressbar(0);
        for k=2:numel(tauVec)
            fun = @(ef)double(kepler(tauVec(k), ef));
            guess = FVec(k-1);
            FVec(k) = fzero(fun, guess);
            progressbar(k/numel(tVec));
        end
        
        % Translate eccentric anomalies back to meaningful space coordinates
        el = q*sqrt(bigG*bigM*(1/a + 2/q));
        aVec = a*(e*cosh(FVec) - 1);
        cosfVec = (el^2./(aVec*bigG*bigM) - 1)/e;
        fVec = acos(cosfVec);
        fVec = fVec.*sign(FVec); % make pre-peri angles negative
        xVec = aVec.*cos(fVec);
        yVec = aVec.*sin(fVec);
        
    case default
        error('Unknown orbit type')
end

%% Save orbit to file
% File with header will be created (and overrun) in the current directory
fid = fopen([orbit_name '.orb'],'wt');
fprintf(fid,'################################################\n');
fprintf(fid,'# Orbit coordinates (in code units)\n');
fprintf(fid,'# Columns are:\n');
fprintf(fid,'# [t]  [x]  [y]\n');
fprintf(fid,'################################################\n');
fclose(fid);
orbit = double([tVec' xVec' yVec']);
save([orbit_name '.orb'],'orbit','-ascii','-append')
