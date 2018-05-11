function [u,xi] = eval_pStokeslet2dwgt_dir(x_tar,N_tar, x_src, N_src, f_src, ...
                                           L,TOL,NB,P,rwgt,fwgt)



cR = 100;
cF = 1;

xi = (NB/L) * sqrt(log(1/TOL)+ log(cR));
Mmax = (xi*L/pi) * sqrt(log(1/TOL)+ log(cF));
Mmax = round(Mmax);
M  = 2*Mmax;

% %P = 33; 
% P = 21;
m = sqrt(P*pi);

uR=zeros(N_tar,2);
uF=zeros(N_tar,2);

if rwgt
    tstart=tic;
    pinf = 2; % nearby two boxes
    uR = stokesletRsum2d_dir(x_tar, N_tar, x_src, f_src, N_src, L, xi, pinf);
    tR = toc(tstart);
end

if fwgt
    tstart=tic;
    uF = stokesletFsum2d(x_tar,N_tar,x_src,f_src,N_src,L,M,P,xi,m);
    tF = toc(tstart);
end

u = (uR+uF)/(4*pi);

fprintf('Stokeslet\n')
fprintf('xi=%.16f, M=%d, P=%d, m=%.4f\n',xi,M,P,m);
if rwgt
    fprintf('Real=%.4f\n',tR)
end

if fwgt
    fprintf('k-space=%.4f\n',tF)
end
% pinf = 3;
% uR_d = stokesletRsum2d_dir(x_tar,N_tar,x_src,f_src, N_src, L, xi, pinf);
% 
% errR = sqrt(sum((uR-uR_d).^2,2));
% errR = max(errR);
% 
% 
% disp('Stokeslet');
% fprintf('xi=%f, M=%f, P=%f, m=%f\n', xi, M, P, m);
% fprintf('errR=%e\n', errR);
