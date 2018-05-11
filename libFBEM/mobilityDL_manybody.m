
function [UW,relresv, t_gmres, t_exportStressletRsum] = mobilityDL_manybody(geom,geom_dir,TOL, NB, FT,uslip, ...
                                                       StressletRsum_sparse, boxes)

% load configuration
% [x,xc,rs,npts,L,nv,tau,wgts,kappa,cU,cW,Nbody] = load_config(geom, geom_dir);
gc = load_config(geom, geom_dir);
Ntot = sum(gc.np);

x = gc.x_pos;
L = gc.L;
xc = gc.q;
Nbody = gc.nb;
npts = gc.np;
tau = gc.tv;
nv = gc.nv;
kappa = gc.kappa;
cU=gc.cU;
cW=gc.cW;
wgts = gc.wgts;

Pspread = 21;
paramEwald_xc  = fastEwaldParameters(xc, L, NB , TOL, TOL, 'Stokeslet', Pspread);
paramEwald_x  = fastEwaldParameters(x, L, NB , TOL, TOL, 'Stokeslet', Pspread);


if nargin < 7
    tstart=tic;
    StressletRsum_sparse = exportStressletRsum2d(geom, geom_dir, TOL, NB);
    t_exportStressletRsum = toc(tstart);
    fprintf('elapsed time StressletRsum = %f\n', t_exportStressletRsum);
    
%     % evaluate sum due to Stokeslets and Rotlets
%     idx = @(j) j(1)+(j(2)-1)*NB;
%     boxes = makeBox2d(x,Ntot,L,NB,idx);
    
end 




u_stokeslet = eval_pStokeslet2d(x, Ntot, xc, Nbody, FT(:,1:2), L, paramEwald_xc);
u_rotlet    = eval_pRotlet2d(x, Ntot, xc, Nbody, FT(:,3), L, TOL, NB);

rhs = -(u_stokeslet+u_rotlet);



if nargin < 6
    uslip = zeros(Ntot,2);
end

% right hand side vector
rhs = rhs - uslip;
rhs = rhs';
rhs = rhs(:);

evalDL = @(q) completedDLmobility(q, x, npts, nv, tau, kappa, wgts,...
                                  xc, TOL, L, NB,paramEwald_x.box_arr,cU,cW, StressletRsum_sparse);
                              
% GMRES 
tstart = tic;
maxit = 80; gmresTOL = 1E-12; 
[dens,fl,relres,iter,resv] = gmres(evalDL , rhs, [], gmresTOL , maxit);
t_gmres = toc(tstart);

relresv = resv/norm(rhs);

fprintf('The number of iter = %d\n', length(resv)-1)
fprintf('DL mobility GRMES elapsed time = %.4f seconds \n', t_gmres)

%%
% compute mobility from density
U = zeros(Nbody,2);
W = zeros(Nbody,1);

nstart = 0;
for nbod = 1 : Nbody
    
    nb = npts(nbod);
    
    
    i0 = nstart/2 + (1:nb);
    i1 = nstart + (1:2:2*nb);
    i2 = nstart + (2:2:2*nb);
    

    
    rperp = [-x(i0,2)+xc(nbod,2), x(i0,1)-xc(nbod,1)];
    
    
    qj = [dens(i1).*wgts(i0), dens(i2).*wgts(i0) ];
    U(nbod,:) = 1/cU(nbod) * sum(qj);
    W(nbod) = 1/cW(nbod) * sum( qj(:,1).*rperp(:,1) + qj(:,2).*rperp(:,2));
    
    nstart=nstart+2*nb;
    
end


UW = [U,W];
UW = UW'; UW = UW(:);
