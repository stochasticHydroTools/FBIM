function uF = stressletFsum2d_dir(x_tar, N_tar, x_src, q_src, nv_src, ...
                                  N_src, L, M, xi)

uF = zeros(N_tar,2);

V = L^2;

wnum = (-M : M) * (2*pi/L);
[k2,k1] = meshgrid(wnum, wnum);
kn2 = k1.^2 + k2.^2;
kn4 = kn2.^2;
e = exp(-kn2/(4*xi^2));
A0   = -4*pi./kn4;
A0(M+1,M+1) = 0;

A1 = 1 + kn2/(4*xi^2);

for m = 1 : N_tar
    
    xm = x_tar(m,:);
    
    q1 = 0; q2 = 0;
    for n = 1 : N_src
        
        xn = x_src(n,:);
        dxmn = xm - xn;
        qn  = q_src(n,:);
        nvn = nv_src(n,:);
        
        kq  = k1.*qn(1) + k2.*qn(2);
        qnv = sum(qn .* nvn);
        knv = k1.*nvn(1) + k2.*nvn(2);
        en   = sin(k1*dxmn(1) + k2*dxmn(2));

        
        q1 = q1 + (A1.*( kn2.*(kq.*nvn(1) + qnv*k1) - 2*k1.*kq.*knv) + kn2*qn(1).*knv ) .* en ; 
        q2 = q2 + (A1.*( kn2.*(kq.*nvn(2) + qnv*k2) - 2*k2.*kq.*knv) + kn2*qn(2).*knv ) .* en; 
       
    end
    
    uF(m,1) = 1/V * sum(sum( A0.* q1 .* e));
    uF(m,2) = 1/V * sum(sum( A0.* q2 .* e));

   
end
