function uF = stressletFsum2d(x_tar, N_tar, x_src, q_src, nv_src, N_src, ...
                              L, M, P, xi, m)

h = L/M;
eta   = (P*h*xi/m)^2;
alpha = 2*xi^2/eta;
spread_ampl = 2*xi^2/(pi*eta);

H  = zeros(M,M,2,2);
H_hat = zeros(M,M,2,2);
uF = zeros(N_tar,2);
Ht_hat = zeros(M,M,2);
Ht = zeros(M,M,2);


%% STEP 1: smear each source with gaussian on the grid
for n = 1 : N_src
   
    xn = x_src(n,:);
    qn = q_src(n,:);
    nv = nv_src(n,:);
    
    
%    [fg,ix,iy] = fastgridding2d(xn,alpha,P,M,h); 

    % use mex
    [fg,ix,iy] = fastgridding2d_mex(xn(1),xn(2), alpha, h, P);
    ix = mod(ix,M) + 1;
    iy = mod(iy,M) + 1;

    for l = 1:2
        for m = 1 : 2;
            H(ix,iy,l,m) = H(ix,iy,l,m) + qn(l)*nv(m) *spread_ampl * fg;
        end
    end
end

%% STEP 2: take the fft of H
for l = 1:2
    for m = 1:2
        H_hat(:,:,l,m) = fft2(H(:,:,l,m));
    end
end



%% STEP 3: multiply H_hat by the decomposition kernel
wnum = fftshift(-M/2 : M/2-1) * (2*pi/L);
[k2,k1] = meshgrid(wnum, wnum);
kn2 = k1.^2 + k2.^2;
kn4 = kn2.^2;
scale_ampl  = exp(-(1-eta)*kn2/(4*xi^2));


k1sq = k1.*k1;
k2sq = k2.*k2;
k1k2 = k1.*k2;
A0   = 4*pi*1i./kn4;
A0(1,1) = 0; % zero mode
for j = 1:2
    
    if j == 1 
        kj = k1;
    else
        kj = k2;
    end
    
    A1 = k1.*H_hat(:,:,1,j) + k2.*H_hat(:,:,2,j);
    A2 = kj.*(H_hat(:,:,1,1)+H_hat(:,:,2,2));
    A3 = kj.*( k1sq.*H_hat(:,:,1,1) + k1k2.*(H_hat(:,:,1,2)+H_hat(:,:,2,1)) + k2sq.*H_hat(:,:,2,2) );
    A4 = k1.*H_hat(:,:,j,1) + k2.*H_hat(:,:,j,2);
    
    Ht_hat(:,:,j) = A0.* ( (1+kn2/(4*xi^2)).*(kn2.*(A1+A2)-2*A3) + kn2.*A4 ) .* scale_ampl;
    
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


