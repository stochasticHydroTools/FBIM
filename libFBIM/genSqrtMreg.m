function [ureg_half,xi,M] = genSqrtMreg(x,TOL, NB, L)


cR = 100;
cF = 1;

xi = (NB/L) * sqrt(log(1/TOL)+ log(cR));
Mmax = sqrt(2)*(xi*L/pi) * sqrt(log(1/TOL)+ log(cF));
Mmax = round(Mmax);
M  = 2*Mmax;

Wk = generateWk(M);
P = 21;
m = sqrt(P*pi);

h = L/M;

eta = (P*h*xi/m)^2;
alpha = 2*xi^2/eta;
spread_ampl = 2*xi^2/(pi*eta);


wnum = fftshift(-M/2 : M/2-1) * (2*pi/L);
[k2,k1] = meshgrid(wnum, wnum);
kn = sqrt(k1.^2 + k2.^2);


tstart = tic;

scale_ampl1 =  sqrt(1+kn.^2/(4*xi^2)) ./ (kn.^2);
scale_ampl1(1,1) = 0; % zero mode

scale_ampl2 = exp(-(1-eta)*kn.^2/(8*xi^2));

Hk = zeros(M,M,2);


Hk(:,:,1) = -scale_ampl1 .* scale_ampl2 .* (-k2) .* Wk;
Hk(:,:,2) = -scale_ampl1 .* scale_ampl2 .* ( k1) .* Wk; 

H(:,:,1) = real(ifft2(Hk(:,:,1)));
H(:,:,2) = real(ifft2(Hk(:,:,2)));


ureg_half = zeros(length(x),2);
for l = 1 : length(x)
    
    xm = x(l,:);
    
%    [a,ix,iy] = fastgridding2d(xm, alpha, P, M, h);
    
    [a,ix,iy] = fastgridding2d_mex(xm(1),xm(2), alpha, h, P);
    ix = mod(ix,M) + 1;
    iy = mod(iy,M) + 1;
    
    a = spread_ampl * a;
    

    % no need to multiply by h^2
    ureg_half(l,1) = sum(sum(H(ix,iy,1) .* a));
    ureg_half(l,2) = sum(sum(H(ix,iy,2) .* a));

end

tend = toc(tstart);
% fprintf('xi = %.16f, M = %d, P = %d, m = %.4f\n', xi, M, P, m);
% fprintf('k-space = %.4f\n',tend)


V = L^2;
ureg_half = ureg_half'; 
ureg_half = sqrt(V) * ureg_half(:);


