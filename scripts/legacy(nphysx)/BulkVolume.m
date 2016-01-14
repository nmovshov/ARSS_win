function V=BulkVolume(xVec,yVec,zVec,rVec,nPoints)
%BULKVOLUME     Estimated bulk volume of a rubble pile
% BULKVOLUME uses a convex hull to estimate the volume of a rubble pile.
% The volume returned is that of the convex hull containing the rubble
% pile, including internal cavities. The rubble pile is defined by center
% locations of the individual particles - xVec, yVec, and zVec, along with
% size information - rVec. If the particles are not spherical then rVec is
% simply a size scale, or, to ensure complete coverage, rVec can be the
% radii of the circumscribing spheres around the individual particles.
%
% Algorithm. Given a set of points in 3d space convhulln can find the
% convex hull of the set and its volume. However, since the particles
% residing at these positions are not point masses, the convex hull based
% on center positions does not actually contain the rubble pile. To get a
% convex hull that contains the rubble pile we take the convex hull of the
% points making the surfaces of those particles who are on the "surface" of
% the pile. So we run convhull once on the center coordinates to find which
% particles make the outside of the rubble pile. Then we generate
% coordinates on the surface of these particles and use convhull with these
% points. The only tricky part is to decide how many points can adequately
% represent the surface of a given particle. NOTE: in this implementation a
% contant numebr of points is used, since I anticipate using this procedure
% mostly on similar sized grains.

%% minimal input filter
error(nargchk(4,5,nargin))
if ~(isnumeric(xVec)&&isnumeric(yVec)&&isnumeric(zVec)&&isnumeric(rVec))
    error('don''t be a dick')
end
if ~(isreal(xVec)&&isreal(yVec)&&isreal(zVec)&&isreal(rVec))
    error('don''t be a dick')
end

%% minimal input massaging
xVec=xVec(:);
yVec=yVec(:);
zVec=zVec(:);
if isscalar(rVec)
    rVec=rVec*ones(size(xVec));
else
    rVec=rVec(:);
end
if nargin==4, nPoints=10; end

%% first pass with convhull
[K,V]=convhulln([xVec yVec zVec]);
fprintf('first pass volume %g\n',V)

%% constrcut shells of outer particles
nPoints=nPoints^2; % make it a sqrt-able number
K=unique(K(:));
Points=[];
phi=linspace(0,pi,fix(sqrt(nPoints)));
theta=linspace(0,2*pi,fix(sqrt(nPoints)));
for i=1:numel(K)
    shell=nan(nPoints,3); % pre-initialize for speed
    for j=1:numel(phi)
        for k=1:numel(theta)
            idx=numel(phi)*(j-1)+k;
            shell(idx,1)=xVec(K(i))+rVec(K(i))*cos(phi(j))*cos(theta(k));
            shell(idx,2)=yVec(K(i))+rVec(K(i))*cos(phi(j))*sin(theta(k));
            shell(idx,3)=zVec(K(i))+rVec(K(i))*sin(phi(j));
        end
    end
    Points=[Points;shell]; %#ok<AGROW>
end

%% second pass with convhull
[~,V]=convhulln(Points);
fprintf('second pass volume %g (%d points per shell)\n',V,nPoints)