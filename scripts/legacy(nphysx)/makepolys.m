function polys=makepolys(n,m,gam)
%MAKEPOLY  Make randomized convex ployhedron
% MAKEPOLY(n) will create n polyhedra with m vertices each, and write them
% to a file called polys.pol. The vertices lie on the unit sphere and no
% two are less than gam degrees (default 10) apart.

if nargin<3, gam=10; end
gam=gam*pi/180;
%s=RandStream('mt19937ar','seed',sum(100*clock)); % randomize timer
polys=cell(n,1);
%% make n polys
for j=1:n
    polyverts=nan(m,3);
    k=0;
    while k<m
        %% make a point on the unit sphere
        x=rand;
        y=rand;
        if (x^2+y^2)>1, continue, end
        z=sqrt(1-x^2-y^2)*sign(s.rand-0.5);
        cand=[x y z];
        %% reject the point if the smallest angular opening for another
        %% point is less than gam
        if any((polyverts*cand')>cos(gam)), continue, end
        %% add point to the list
        k=k+1;
        polyverts(k,:)=cand;
    end
    %% add object to the list
    polys{j}=polyverts;
end
%% write polys to a file.
for j=1:n
    dlmwrite('polys.pol',polys{j},'-append','delimiter',' ','newline','pc');
end
