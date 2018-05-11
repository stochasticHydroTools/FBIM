function uF = rotletFsum2d_dir(x_tar, N_tar, x_src, t_src, N_src, xi, L, M)

uF = zeros(N_tar,2);

wnum = (-M : M) * (2*pi/L);
[k2,k1] = meshgrid(wnum, wnum);
kn2 = k1.^2 + k2.^2;

V = L^2;

ampl0 = 2*pi./kn2 .* (1+kn2/(4*xi^2)) .* exp(-kn2/(4*xi^2));
ampl0(M+1,M+1) = 0; % zero mode

ampl1 = ampl0 .* (-k2);
ampl2 = ampl0 .* k1;


for m = 1 : N_tar
    
    xm = x_tar(m,:);
    
    q = 0; 
    for n = 1 : N_src
       xn = x_src(n,:);
       dxmn = xm - xn;
       tn = t_src(n);
       en = sin( k1.*dxmn(1)+k2.*dxmn(2));
       q = q + en * tn;
    end
    
    s1 = ampl1 .* q;
    s2 = ampl2 .* q;
    
    uF(m,1) = 1/V * sum(sum(s1));
    uF(m,2) = 1/V * sum(sum(s2));
    
    
end