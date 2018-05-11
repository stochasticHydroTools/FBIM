%--------------------------------------------------------------------------
% inputf: geometry input file
% iFT   : force, torque input file
% nb : number of points on each body
% norder: order of Alpert quadrature7
% sflg : slip velocity flag
%--------------------------------------------------------------------------


function [UW_1st, res_1st, t_GRMES, t_Rsum] = ...
         mobilitySL_manybody(geom, geom_dir, geom_ref_dir, output_dir, norder, TOL, NB, FT,uslip, ...
                             StokesletRsum2d_sparse, Msing_alpert_ref, pinvM)


%% load geometry
[x,xc,rs,npts,L,nv,tau,wgts,kappa,cU,cW,Nbody] = load_config(geom, geom_dir);
Ntot = sum(npts);
ht = 2*pi./npts;


if nargin <= 9
    %% load matrices
    [xi,r_cutoff,idx,boxes,M,P,m_shape] = fastEwaldParameters(x,L,NB,TOL,'Stokeslet');
    
    Msing_trap_ref = exportMsing_trap_ref(geom,geom_ref_dir, TOL, NB, output_dir) / ht(1);
    Mreg_ref = exportMreg(geom, geom_ref_dir, TOL, NB, output_dir) / ht(1);

    t0 = tic;
    StokesletRsum2d_sparse = exportStokesletRsum2d(geom, geom_dir, TOL, NB);
    t_Rsum = toc(t0);

    %% build preconditioner
    Msing_alpert_ref = importdata(sprintf('%sndisk_input_ref_MsingAlpert_B01N%03d_order%02d_xi%03d.dat', output_dir, npts(1), norder, floor(xi))) / ht(1);
    M1st_ref = (Msing_alpert_ref + Msing_trap_ref + Mreg_ref);

    [eV,eD] = eig(M1st_ref);
    eDinv = 1./diag(eD);
    for j = 1:length(diag(eD))
        if eD(j,j) < 0 || abs(eD(j,j)) < 1e-6 % set negative or spurious eigenvalues to zero
            eDinv(j) = 0;
            eD(j,j) = 0; 
        end
    end
    pinvM = real(eV * diag(eDinv) / eV);
end

evalF = @(q) eval_nbody1stkind(q, x, xc, npts, Nbody, Msing_alpert_ref, L, NB, TOL, StokesletRsum2d_sparse);
applyP = @(q) blkdiag_precond_1stkind(q, pinvM, x, xc, Nbody, npts );

FT2 = FT'; FT2 = FT2(:);
slip = uslip'; slip = slip(:);
qq = [slip ; -FT2];

maxit = 80; gmresTOL = 1e-9; gmresRestart = [];

t0 = tic;
[y1, flg1, relres1, iter1, resvec1] = gmres(evalF,qq, gmresRestart, gmresTOL, maxit, applyP); 
t_GRMES = toc(t0);

res_1st = resvec1 / norm(applyP(qq));
UW_1st = y1(2*Ntot+1:end);


end
% if strcmp(geom,'ndisk') 
%     
%     fid = fopen(sprintf('%s_input',geom));
%     
%     tline =  fgetl(fid);
%     td = str2num(tline);
%     L = td(1);
%     Nbody = td(2);
% 
%     npts = zeros(1,Nbody);
%     rs = zeros(1,Nbody);
%     xc = zeros(Nbody,2);
%     ht = zeros(1,Nbody);
%     
%     for nbod = 1:Nbody
%         tline = fgetl(fid);
%         td = str2num(tline);
%         
%         npts(nbod) = td(1);
%         rs(nbod) = td(2);
%         xc(nbod,1) = td(3);
%         xc(nbod,2) = td(4);
%         ht(nbod) = 2*pi/npts(nbod);
%     end
%     
%     if nargin == 2
%         fid2 = fopen(sprintf('%s_FT',geom));
%         FT = zeros(Nbody,3);
%         for nbod = 1:Nbody
%             td2 = fgetl(fid2);
%             FT(nbod,:) = str2num(td2);
%             
%         end
%         fclose(fid2);
%     end
%     
% 
%     Ntot = sum(npts);
%     x = zeros(Ntot,2);
%     tau = zeros(Ntot,2);
%     nv = zeros(Ntot,2);
%     kappa = zeros(Ntot,1);
%     wgts = zeros(Ntot,1);
% 
% 
%     nstart=0;
%     for nbod = 1:Nbody
%         
%         i0 = nstart + (1:npts(nbod));
%         t = (0:npts(nbod)-1)' * ht(nbod);
% 
%         xx = rs(nbod)*cos(t); 
%         yy = rs(nbod)*sin(t);
% 
%         x(i0,:) = [xx+xc(nbod,1), yy+xc(nbod,2)];
%         xp = -yy;
%         yp =  xx;
%         nx = yp; 
%         ny = -xp;
%         dsdt = sqrt(nx.^2+ny.^2);
%         tau(i0,:) = [xp./dsdt, yp./dsdt];
%         nv(i0,:) = [nx./dsdt, ny./dsdt];
% 
%         xpp = -xx;
%         ypp = -yy;
%         kappa(i0,:) = (xp.*ypp-yp.*xpp) ./ dsdt.^3;
%         wgts(i0,:) = dsdt*ht(nbod);
% 
%         nstart = nstart+ npts(nbod);
% 
%     end;
% 
% elseif strcmp(geom,'nstarfish')
%     
%     
%     fid = fopen(sprintf('%s_input',geom));
%     fid2 = fopen(sprintf('%s_FT',geom));
%     
%     tline =  fgetl(fid);
%     td = str2num(tline);
%     L = td(1);
%     Nbody = td(2);
%     
%     npts = zeros(1,Nbody);
%     rs = zeros(1,Nbody);
%     as = zeros(1,Nbody);
%     omegas = zeros(1,Nbody);
%     xc = zeros(Nbody,2);
%     ht = zeros(1,Nbody);
%     FT = zeros(Nbody,3);
% 
%     for nbod = 1:Nbody
%         tline = fgetl(fid);
%         td = str2num(tline);
%         
%         td2 = fgetl(fid2);
%         FT(nbod,:) = str2num(td2);
% 
%         npts(nbod) = td(1);
%         rs(nbod) = td(2);
%         as(nbod) = td(3);
%         omegas(nbod) = td(4);
%         xc(nbod,1) = td(5);
%         xc(nbod,2) = td(6);
%         ht(nbod) = 2*pi/npts(nbod);
%     end 
%     
%     Ntot = sum(npts);
%     x = zeros(Ntot,2);
%     tau = zeros(Ntot,2);
%     nv = zeros(Ntot,2);
%     kappa = zeros(Ntot,1);
%     wgts = zeros(Ntot,1);
%     
%     nstart=0;
%     for nbod = 1:Nbody
%         
%         i0 = nstart + (1:npts(nbod));
%         t = (0:npts(nbod)-1)' * ht(nbod);
%         omega=omegas(nbod);
%         a=as(nbod);
%         r1=rs(nbod);
%         
%         rho  = 1+a*cos(omega*t);
%         rhop = -a*omega*sin(omega*t);
%         rhopp= -a*omega^2*cos(omega*t);
% 
%         xx = r1*cos(t); yy = r1*sin(t);
% 
%         x(i0,:) = [xc(nbod,1) + rho.*xx, xc(nbod,2)+rho.*yy];
%         xp = rhop.*xx + rho.*(-yy);
%         yp = rhop.*yy + rho.*xx;
%         nx = yp;
%         ny = -xp;
%         dsdt = sqrt(nx.^2+ny.^2);
%         tau(i0,:) = [xp./dsdt,yp./dsdt];
%         nv(i0,:) = [nx./dsdt,ny./dsdt];
%         xpp = rhopp.*xx - 2*rhop.*yy - rho.*xx;
%         ypp = rhopp.*yy + 2*rhop.*xx - rho.*yy;
% 
%         kappa(i0,:) = (xp.*ypp-yp.*xpp) ./ dsdt.^3;
%         wgts(i0,:) = dsdt*ht(nbod);
%         
%         nstart = nstart + npts(nbod);
%         
%     end
% 
% end
% 
% fclose(fid);
% 
% % 
% Msing_trap = importdata(sprintf('%sMsingtrap_N%d.mat',geom,Ntot));
% Msing_alpert = importdata(sprintf('%sMsing_alpert_N%03d_order%02d',geom,Ntot,norder));
% Mreg  = importdata(sprintf('%sMreg_N%d.mat',geom,Ntot));
% 
% % dist = xc(2,1)-xc(1,1)-2*rs(1);
% % 
% % Msing_alpert = importdata(sprintf('./twobodytest/%sMsing_alpert_N%03d_order%02d_dist%s',geom,npts(1),norder,num2str(dist)));
% % Msing_alpert = (Msing_alpert+Msing_alpert')/2;
% % Mreg = importdata(sprintf('./twobodytest/%sMreg_N%03d_dist%s.mat',geom,npts(1),num2str(dist)));
% % Msing_trap = importdata(sprintf('./twobodytest/%sMsingtrap_N%03d_dist%s.mat',geom,npts(1),num2str(dist)));
% 
% M1stkind = (Msing_trap + Msing_alpert + Mreg)/ht(1);
% 
% K = zeros(2*Ntot,3*Nbody);
% nstart = 0;
% for nbod = 1:Nbody
%     Ksub = zeros(2*npts(nbod),3);
%     v1 = [ones(1,npts(nbod)) ; zeros(1,npts(nbod))];
%     Ksub(1:end,1) = v1(:);
%     Ksub(1:end,2) = flipud(v1(:));
%     i0 = (1:npts(nbod)) + nstart;
%     Ksub(1:2:end,3) = -(x(i0,2)-xc(nbod,2));
%     Ksub(2:2:end,3) =  (x(i0,1)-xc(nbod,1));
%     
%     K(2*nstart+(1:2*npts(nbod)),3*nbod-2:3*nbod) = Ksub;
%     nstart = nstart + npts(nbod);
% end
% 
% FT = FT'; FT=FT(:);
% 
% rhs = [zeros(2*Ntot,1) ; -FT];
% A = [M1stkind -K ; -K' zeros(3*Nbody,3*Nbody)];
% Y=A\rhs;
% 
% UW = Y(end-(3*Nbody-1):end);
