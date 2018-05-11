function Boxgrid = makeBox2d(xsc,Nsc,L,NB,idx)

hB = L/NB;

f1 = 'id'; % box id
f2 = 'N';  % number of sources in box
f3 = 'index_x'; % index of each source in xsc, fsc
f4 = 'pt';

% pre-allocate
Boxgrid(NB^2) = struct(f1,[],f2,[],f3,[],f4,1);

for j = 1:NB
for i = 1:NB
    Boxgrid(idx([i,j])) = struct(f1,[i-1,j-1],f2,0,f3,[],f4,1);
end
end

% loop through the sources and determine how many sources in each box
for n = 1 : Nsc
    xn  = xsc(n,:);
    BoxID = floor(mod(xn,L)/hB);
    j = BoxID + 1;
    
    Boxgrid(idx(j)).N = Boxgrid(idx(j)).N + 1;
    
end

% loop through boxes to preallocate space
for k = 1 : NB^2
   nsc = Boxgrid(k).N;
   Boxgrid(k).index_x = zeros(1,nsc);
end

% loop though sources to determine exactly which source belongs to which
% box
for n = 1 : Nsc

    xn  = xsc(n,:);
    BoxID = floor(mod(xn,L)/hB);
    j = BoxID + 1;
    k = Boxgrid(idx(j)).pt;
    
    Boxgrid(idx(j)).index_x(k) = n;
    
    Boxgrid(idx(j)).pt = k + 1;
    
end