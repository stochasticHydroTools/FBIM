function u = eval_pStresslet2d(x_tar, N_tar, x_src, N_src, q_src, nv_src, ...
                               L, tol, NB, Q, boxes, Rsum_SPmatrix)
                           
                           
V = L^2;
rc = L/NB;
cR = 100;
cF = 100;

xi = (NB/L) * sqrt(log(1/tol)+ log(cR));
Mmax = (xi*L/pi) * sqrt(log(1/tol)+ log(cF));
Mmax = round(Mmax);
M  = 2*Mmax;

% Mmax = (2*xi) * sqrt(log(sqrt(Q/V)/tol)) * (L/2/pi);
% Mmax = round(Mmax);
% M = 2*Mmax;

%P = 33; 
P = 21;
m = sqrt(pi*P);

idx = @(j) j(1)+(j(2)-1)*NB;
%boxes = makeBox2d(x_src,N_src,L,NB,idx);

t1 = tic;
if nargin == 12
    qq = q_src'; qq = qq(:);
    uR = Rsum_SPmatrix * qq;
    uR = reshape(uR, [2,N_tar])';
else
    uR = stressletRsum2d(x_tar,N_tar,x_src,q_src,nv_src,boxes, xi, rc, NB, L,idx);
end
tR = toc(t1);

t2 = tic;
uF = stressletFsum2d(x_tar,N_tar,x_src,q_src,nv_src,N_src,L,M,P,xi,m);
tF = toc(t2);

u = (uR+uF)/(4*pi);

fprintf('Stresslet\n')
fprintf('xi=%.8f, M=%d, P=%d, m=%.4f\n',xi,M,P,m);
fprintf('Real=%.4f, k-space=%.4f\n\n',tR,tF)
%   
% pinf = 3;
% uR_d = stressletRsum2d_dir(x_tar,N_tar,x_src,q_src,nv_src,N_src,L,xi,pinf);
% uF_d = stressletFsum2d_dir(x_tar,N_tar,x_src,q_src,nv_src,N_src,L,Mmax,xi);
%  
%  
% errR = sqrt(sum((uR-uR_d).^2,2));
% errR = max(errR);
% errF = sqrt(sum((uF-uF_d).^2,2)); errF = max(errF);
% disp('Stresslet');
% fprintf('xi=%f, M=%f, P=%f, m=%f\n', xi, M, P, m);
% fprintf('errR=%e, errF=%e\n', errR, errF);
