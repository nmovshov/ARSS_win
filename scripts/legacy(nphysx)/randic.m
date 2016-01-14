function x=randic(n,L,r)
%RANDIC    randomized non-intersecting positions
% randic(n,L,r) returns a list of positions (n-by-3 matrix) inside a 3D
% space of extent ±L in each dimension, such that each point is at least 2r
% units away from each other point. No serious attempt at optimization was
% made, so assume complexity is O(n^2).

if nargin<3, error('don''t be a dick'); end
if fix(n)~=n, error('don''t be a dick'); end
if (2*L)^3<(n*4*pi/3*r^3), error('that''s impossible, don''t be a dick'); end
if (2*L)^3<2*(n*4*pi/3*r^3)
    warning('that''s highly unlikely, consider aborting') %#ok<WNTAG>
end
x=zeros(n,3);
good=0;
while good<n
    x0=-L+2*L*rand(1,3);
    diffmat=x(1:good,:)-repmat(x0,good,1);
    if any(sqrt(dot(diffmat,diffmat,2))<(2*sqrt(2)*r)), continue, end
    good=good+1;
    x(good,:)=x0;
end