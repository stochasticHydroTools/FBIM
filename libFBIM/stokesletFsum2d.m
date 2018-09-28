function uF = stokesletFsum2d(x_tar,N_tar,x_src,f_src,N_src,L,M,P,xi,m)

h = L/M;

eta = (P*h*xi/m)^2;
alpha = 2*xi^2/eta;
spread_ampl = 2*xi^2/(pi*eta);

H  = zeros(M,M,2);
uF = zeros(N_tar,2);

%% STEP 1: smear each source with gaussian on the grid
% for n = 1 : N_src
%    
%     xn = x_src(n,:);
%     fn = f_src(n,:);
%     
%     fg = fastgridding2d(xn,alpha,P,M,h);
%     
%     H(:,:,1) = H(:,:,1) + fn(1) * spread_ampl * fg;
%     H(:,:,2) = H(:,:,2) + fn(2) * spread_ampl * fg;
%     
% end

for n = 1 : N_src
   
    xn = x_src(n,:);
    fn = f_src(n,:);

%     [fg,ix,iy] = fastgridding2d(xn,alpha,P,M,h);
    
    % use mex
    [fg,ix,iy] = fastgridding2d_mex(xn(1),xn(2), alpha, h, P);
    ix = mod(ix,M) + 1;
    iy = mod(iy,M) + 1;
    
    H(ix,iy,1) = H(ix,iy,1) + fn(1) * spread_ampl * fg;
    H(ix,iy,2) = H(ix,iy,2) + fn(2) * spread_ampl * fg;
    
end


%% STEP 2: take the fft of H
r1 = H(:,:,1);
H_hat1 = fft2(r1);
r2 = H(:,:,2);
H_hat2 = fft2(r2);



%% STEP 3: multiply H_hat by the decomposition kernel
wnum = fftshift(-M/2 : M/2-1) * (2*pi/L);
[k2,k1] = meshgrid(wnum, wnum);
kn = sqrt(k1.^2 + k2.^2);
scale_ampl  = exp(-(1-eta)*kn.^2/(4*xi^2));

% Hasimoto
scale_ampl2 = 4*pi*(1+kn.^2/(4*xi^2)) ./ (kn.^4);
scale_ampl2(1,1) = 0; % zero mode

kH = k1.*H_hat1 + k2.*H_hat2;
Ht_hat1 = scale_ampl .* scale_ampl2 .* ( kn.^2 .* H_hat1 - k1.*kH);
Ht_hat2 = scale_ampl .* scale_ampl2 .* ( kn.^2 .* H_hat2 - k2.*kH);

%% STEP 4: take the iFFT to get Ht
Ht(:,:,1) = real(ifft2(Ht_hat1));
Ht(:,:,2) = real(ifft2(Ht_hat2));

%% STEP 5: evalute integral with Trapezoidal rule to get uF for each source
% for l = 1 : N_tar
%     
%     xm = x_tar(l,:);
%     
%     a = fastgridding2d(xm, alpha, P, M, h);
%     a = spread_ampl * a;
% 
%     uF(l,1) = sum(sum(Ht(:,:,1) .* a))*h^2;
%     uF(l,2) = sum(sum(Ht(:,:,2) .* a))*h^2;
% 
% end

for l = 1 : N_tar
    
    xm = x_tar(l,:);
    
%     [a,ix,iy] = fastgridding2d(xm, alpha, P, M, h);

    % use mex
    [a,ix,iy] = fastgridding2d_mex(xm(1),xm(2), alpha, h, P);
    ix = mod(ix,M) + 1;
    iy = mod(iy,M) + 1;
    
    a = spread_ampl * a;

    uF(l,1) = sum(sum(Ht(ix,iy,1) .* a))*h^2;
    uF(l,2) = sum(sum(Ht(ix,iy,2) .* a))*h^2;

end

