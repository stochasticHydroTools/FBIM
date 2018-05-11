function u = eval_pStokeslet2d(x_tar,N_tar, x_src, N_src, f_src, L, paramEwald, StokesletRsum_sparse)

% cR = 100;
% cF = 1;
% 
% xi = (NB/L) * sqrt(log(1/tol)+ log(cR));
% Mmax = (xi*L/pi) * sqrt(log(1/tol)+ log(cF));
% Mmax = round(Mmax);
% M  = 2*Mmax;
% 
% rc = L/NB;
% 
% %P = 33; 
% %P = 21;
% m = sqrt(P*pi);
% 
% 
% idx = @(j) j(1)+(j(2)-1)*NB;
% boxes = makeBox2d(x_src,N_src,L,NB,idx);

xi=paramEwald.xi;
NB=paramEwald.nbox;
M = paramEwald.Mfourier;
rc = paramEwald.r_cutoff;
m = paramEwald.m_shape;
idx = paramEwald.pts_idx;
boxes=paramEwald.box_arr;
P=paramEwald.Pspread;

tstartR=tic;
if nargin == 8
    ff = f_src'; ff = ff(:);
    uR = StokesletRsum_sparse * ff;
    uR = reshape(uR,[2,N_tar])';
else
    uR = stokesletRsum2d(x_tar,N_tar,x_src,f_src,boxes, xi, rc, NB, L,idx);
end
tR = toc(tstartR);

tstartF=tic;
uF = stokesletFsum2d(x_tar,N_tar,x_src,f_src,N_src,L,M,P,xi,m);
tF = toc(tstartF);

u = (uR+uF)/(4*pi);

% fprintf('Stokeslet\n')
% fprintf('xi=%.16f, M=%d, P=%d, m=%.4f\n',xi,M,P,m);
% fprintf('Real=%.4f, k-space=%.4f\n',tR,tF)

% pinf = 3;
% uR_d = stokesletRsum2d_dir(x_tar,N_tar,x_src,f_src, N_src, L, xi, pinf);
% uF_d = stokesletFsum2d_dir(x_tar,N_tar,x_src,f_src,N_src,L,M,xi);
% 
% errR = sqrt(sum((uR-uR_d).^2,2));
% errR = max(errR);
% errF = sqrt(sum((uF-uF_d).^2,2));
% errF = max(errF);
% 
% disp('Stokeslet');
% fprintf('xi=%f, M=%f, P=%f, m=%f\n', xi, M, P, m);
% fprintf('errR=%e, errF=%e\n', errR, errF);

