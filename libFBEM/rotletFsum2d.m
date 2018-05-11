function uF = rotletFsum2d(x_tar, N_tar, x_src, t_src, N_src, L, M, P, xi, m)

h = L/M;

eta = (P*h*xi/m)^2;
alpha = 2*xi^2/eta;
spread_ampl = 2*xi^2/(pi*eta);

H = zeros(M,M);
uF = zeros(N_tar,2);
Ht_hat = zeros(M,M,2);
Ht = zeros(M,M,2);


%% STEP 1: smear each source with gaussian on the grid
for n = 1 : N_src
   
    xn = x_src(n,:);
    tn = t_src(n);
    
    % [fg,ix,iy] = fastgridding2d(xn,alpha,P,M,h);
    
    % use mex
    [fg,ix,iy] = fastgridding2d_mex(xn(1),xn(2), alpha, h, P);
    ix = mod(ix,M) + 1;
    iy = mod(iy,M) + 1;
    
    
    H(ix,iy) = H(ix,iy) + tn * spread_ampl * fg;
    
end


%% STEP 2: take the fft of H
H_hat = fft2(H);

%% STEP 3: multiply H_hat by the decomposition kernel
wnum = fftshift(-M/2 : M/2-1) * (2*pi/L);
[k2,k1] = meshgrid(wnum, wnum);
kn2 = k1.^2 + k2.^2;
scale_ampl = exp(-(1-eta)*kn2/(4*xi^2));

scale_ampl2 = -2*1i*pi./kn2.*(1+kn2/(4*xi^2));
scale_ampl2(1,1) = 0;

for j = 1 : 2;
    
    if j == 1
        kj = -k2;
    else
        kj = k1;
    end
    
    Ht_hat(:,:,j) = scale_ampl .* scale_ampl2 .* kj .* H_hat;
    
end

%% STEP 4: take the iFFT to get Ht
Ht(:,:,1) = real(ifft2(Ht_hat(:,:,1)));
Ht(:,:,2) = real(ifft2(Ht_hat(:,:,2)));


%% STEP 5: evalute integral with Trapezoidal rule to get uF for each source
for l = 1 : N_tar
    
    xm = x_tar(l,:);
    
%    [a,ix,iy] = fastgridding2d(xm, alpha, P, M, h);
    
    % use mex
    [a,ix,iy] = fastgridding2d_mex(xm(1),xm(2), alpha, h, P);
    ix = mod(ix,M) + 1;
    iy = mod(iy,M) + 1;
    
    a = spread_ampl * a;

    uF(l,1) = sum(sum(Ht(ix,iy,1) .* a))*h^2;
    uF(l,2) = sum(sum(Ht(ix,iy,2) .* a))*h^2;

end
