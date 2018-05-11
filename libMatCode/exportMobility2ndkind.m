function exportMobility2ndkind( geom, geom_dir, TOL,NB, output_dir )


fname = fopen(geom_dir);
tline=fgetl(fname); 
td = str2num(tline);
Nbody = td(2);
tline=fgetl(fname);
td = str2num(tline);
npts = td(1); 
% r = td(2); x1 = td(3);
% tline = fgetl(fname); td = str2num(tline);
% x2 = td(3);
% dist = x2-x1-2*r;


N = 3*Nbody;
I = eye(N);
N2nd = zeros(N,N);

% zero slip
for n = 1:3*Nbody
   
    F = I(n,:);
    F = reshape(F,[3,Nbody])';
    
    N2nd(:,n) = mobilityDL_manybody(geom,geom_dir,TOL,NB,F);
    
    
end

save(sprintf('%s/%s_mob2nd_B%02dN%03d.mat',output_dir,geom,Nbody,npts),'N2nd')
