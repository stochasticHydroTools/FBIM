clear;
addpath('../FirstKindMob/')
addpath('../periodicStokes');
addpath('./testSqrtMsing');


geom = 'nstarfish';
fid = fopen(sprintf('%s_input',geom));
tline = fgetl(fid);
td = str2num(tline);
L = td(1);
Nbody = td(2);

if strcmp(geom,'ndisk') 

    npts = zeros(1,Nbody);
    rs = zeros(1,Nbody);
    xc = zeros(Nbody,2);
    phi = zeros(1,Nbody);
    ht = zeros(1,Nbody);

    for nbod = 1:Nbody
        tline = fgetl(fid);
        td = str2num(tline);

        npts(nbod) = td(1);
        rs(nbod) = td(2);
        xc(nbod,1) = td(3);
        xc(nbod,2) = td(4);
        phi(nbod) = td(5);
        ht(nbod) = 2*pi/npts(nbod);
    end

    Ntot = sum(npts);
    x = zeros(Ntot,2);
    tau = zeros(Ntot,2);
    nv = zeros(Ntot,2);
    kappa = zeros(Ntot,1);
    wgts = zeros(Ntot,1);


    nstart=0;
    for nbod = 1:Nbody
        
        i0 = nstart + (1:npts(nbod));
        t = (0:npts(nbod)-1)' * ht(nbod);

        xx = rs(nbod)*cos(t+phi(nbod)); 
        yy = rs(nbod)*sin(t+phi(nbod));

        x(i0,:) = [xx+xc(nbod,1), yy+xc(nbod,2)];
        xp = -yy;
        yp =  xx;
        nx = yp; 
        ny = -xp;
        dsdt = sqrt(nx.^2+ny.^2);
        tau(i0,:) = [xp./dsdt, yp./dsdt];
        nv(i0,:) = [nx./dsdt, ny./dsdt];

        xpp = -xx;
        ypp = -yy;
        kappa(i0,:) = (xp.*ypp-yp.*xpp) ./ dsdt.^3;
        wgts(i0,:) = dsdt*ht(nbod);

        nstart = nstart+ npts(nbod);

    end;

elseif strcmp(geom,'nstarfish')
    
    npts = zeros(1,Nbody);
    rs = zeros(1,Nbody);
    as = zeros(1,Nbody);
    omegas = zeros(1,Nbody);
    xc = zeros(Nbody,2);
    phi= zeros(1,Nbody);
    ht = zeros(1,Nbody);
    
    for nbod = 1:Nbody
        tline = fgetl(fid);
        td = str2num(tline);

        npts(nbod) = td(1);
        rs(nbod) = td(2);
        as(nbod) = td(3);
        omegas(nbod) = td(4);
        xc(nbod,1) = td(5);
        xc(nbod,2) = td(6);
        phi(nbod) = td(7);
        ht(nbod) = 2*pi/npts(nbod);
    end 
    
    Ntot = sum(npts);
    x = zeros(Ntot,2);
    tau = zeros(Ntot,2);
    nv = zeros(Ntot,2);
    kappa = zeros(Ntot,1);
    wgts = zeros(Ntot,1);
    
        
    nstart=0;
    for nbod = 1:Nbody
        
        i0 = nstart + (1:npts(nbod));
        t = (0:npts(nbod)-1)' * ht(nbod);
        omega=omegas(nbod);
        a=as(nbod);
        r1=rs(nbod);
        
        rho  = 1+a*cos(omega*t);
        rhop = -a*omega*sin(omega*t);
        rhopp= -a*omega^2*cos(omega*t);

        xx = r1*cos(t+phi(nbod)); yy = r1*sin(t+phi(nbod));

        x(i0,:) = [xc(nbod,1) + rho.*xx, xc(nbod,2)+rho.*yy];
        xp = rhop.*xx + rho.*(-yy);
        yp = rhop.*yy + rho.*xx;
        nx = yp;
        ny = -xp;
        dsdt = sqrt(nx.^2+ny.^2);
        tau(i0,:) = [xp./dsdt,yp./dsdt];
        nv(i0,:) = [nx./dsdt,ny./dsdt];
        xpp = rhopp.*xx - 2*rhop.*yy - rho.*xx;
        ypp = rhopp.*yy + 2*rhop.*xx - rho.*yy;

        kappa(i0,:) = (xp.*ypp-yp.*xpp) ./ dsdt.^3;
        wgts(i0,:) = dsdt*ht(nbod);
        
        nstart = nstart + npts(nbod);
        
    end

end

fclose(fid);

figure(1); clf;
plot(x(:,1),x(:,2),'bo'); hold on;
quiver( x(:,1), x(:,2), nv(:,1), nv(:,2) )
axis equal
axis([0,L,0,L]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TOL = 1E-16;
NB = 2;


epsilon = TOL/10;
cR = 1/10;
cF = 1;

xi = (NB/L) * sqrt(log(1/epsilon)+ log(cR));
Mmax = (xi*L/pi) * sqrt(log(1/epsilon)+ log(cF));
Mmax = round(Mmax);
M  = 2*Mmax;

Wk = generateWk(M);

P = 33;
m = sqrt(P*pi);

h = L/M;

eta = (P*h*xi/m)^2;
alpha = 2*xi^2/eta;
spread_ampl = 2*xi^2/(pi*eta);


wnum = fftshift(-M/2 : M/2-1) * (2*pi/L);
[k2,k1] = meshgrid(wnum, wnum);
kn = sqrt(k1.^2 + k2.^2);


scale_ampl1 =  sqrt(1+kn.^2/(4*xi^2)) ./ (kn.^2);
scale_ampl1(1,1) = 0; % zero mode

scale_ampl2 = exp(-(1-eta)*kn.^2/(8*xi^2));
scale_ampl3 = exp(-eta*kn.^2/(8*xi^2));

Hk = zeros(M,M,2);

Hk(:,:,1) = -scale_ampl1 .* scale_ampl2 .* (-k2) .* Wk;
Hk(:,:,2) = -scale_ampl1 .* scale_ampl2 .* ( k1) .* Wk; 


H(:,:,1) = real(ifft2(Hk(:,:,1)));
H(:,:,2) = real(ifft2(Hk(:,:,2)));


uF = zeros(Ntot,2);
for l = 1 : Ntot
    
    xm = x(l,:);
    
    a = fastgridding2d(xm, alpha, P, M, h);
    a = spread_ampl * a;

    % no need to multiply by h^2
    uF(l,1) = sum(sum(H(:,:,1) .* a));
    uF(l,2) = sum(sum(H(:,:,2) .* a));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% 
% uF_dir = zeros(Ntot,2);
% 
% wnum = (-M : M) * (2*pi/L);
% [k2,k1] = meshgrid(wnum, wnum);
% kn = sqrt(k1.^2 + k2.^2);
% 
% % Hasimoto
% scale_ampl2 = sqrt(1+kn.^2/(4*xi^2)) ./ (kn.^2);
% 
% scale_ampl2(M+1,M+1) = 0; % zero mode
% 
% V = L^2;
% 
% % precompute
% e = exp(-kn.^2/(8*xi^2));
% Hk1 = scale_ampl2 .* (-k2) .* e;
% Hk2 = scale_ampl2 .* ( k1) .* e;
% 
% 
% for m = 1 : Ntot
%     
%    xm = x(m,:);
%    
%    e_xm = exp(-1i*(k1*xm(1) + k2*xm(2)) );
%    
%    uF_dir(m,1) = sum(sum(Hk1 .* e_xm));
%    uF_dir(m,2) = sum(sum(Hk2 .* e_xm));
%    
%     
% end
% 
% 
% imag(uF_dir)



