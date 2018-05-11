%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eval_pRotlet2d.m
% 
% evaluate flow induced by Rotlets at source locations at target locations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function u = eval_pRotlet2d(x_tar, N_tar, x_src, N_src, t_src, L, tol, NB)

cR = 1000;
cF = 10;

xi = (NB/L) * sqrt(log(1/tol)+ log(cR));
Mmax = (xi*L/pi) * sqrt(log(1/tol)+ log(cF));
Mmax = round(Mmax);
M  = 2*Mmax;

rc = L/NB;

%P = 33; 
P = 21;
m = sqrt(P*pi);



idx = @(j) j(1)+(j(2)-1)*NB;
boxes = makeBox2d(x_src,N_src,L,NB,idx);

tic;
uR = rotletRsum2d(x_tar,N_tar,x_src,t_src,boxes, xi, rc, NB, L,idx);
tR = toc;

tic;
uF = rotletFsum2d(x_tar,N_tar,x_src,t_src,N_src,L,M,P,xi,m);
tF = toc;

u = (uR+uF)/(4*pi);

fprintf('Rotlet\n')
fprintf('xi=%.8f, M=%d, P=%d, m=%.4f\n',xi,M,P,m);
fprintf('Real=%.4f, k-space=%.4f\n',tR,tF)

% pinf = 3;
% uR_d = rotletRsum2d_dir(x_tar, N_tar, x_src, t_src, N_src, L, xi, pinf);
% uF_d = rotletFsum2d_dir(x_tar, N_tar, x_src, t_src, N_src, xi, L, M);
% 
% 
% errR = sqrt(sum((uR-uR_d).^2,2));
% errR = max(errR);
% 
% errF = sqrt(sum((uF-uF_d).^2,2));
% errF = max(errF);
% disp('Rotlet');
% fprintf('xi=%f, M=%f, P=%f, m=%f\n', xi, M, P, m);
% fprintf('errR=%e, errF=%e\n', errR, errF);

