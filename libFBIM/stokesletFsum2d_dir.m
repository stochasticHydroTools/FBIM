function uF = stokesletFsum2d_dir(x_tar, N_tar, x_src, f_src, N_src, L, M, xi)

uF = zeros(N_tar,2);

wnum = (-M : M) * (2*pi/L);
[k2,k1] = meshgrid(wnum, wnum);
kn = sqrt(k1.^2 + k2.^2);

% Hasimoto
scale_ampl2 = 4*pi*(1+kn.^2/(4*xi^2)) ./ (kn.^4);

scale_ampl2(M+1,M+1) = 0; % zero mode

V = L^2;

% precompute
e = exp(-kn.^2/(4*xi^2));
Bkn = scale_ampl2 .* kn.^2;
Bk1k1 = scale_ampl2 .* k1.*k1;
Bk1k2 = scale_ampl2 .* k1.*k2;

Bk2k1 = scale_ampl2 .* k2.*k1;
Bk2k2 = scale_ampl2 .* k2.*k2;

for m = 1 : N_tar
    
   xm = x_tar(m,:);
   
   q1 = 0; q2 = 0;
   for n = 1 : N_src
      xn = x_src(n,:);
      dxmn = xm - xn;
      fn   = f_src(n,:);
      en   = exp( -1i*(k1*dxmn(1) + k2*dxmn(2)) );
      q1   = q1 + en * fn(1);
      q2   = q2 + en * fn(2);
   end
   
   Bq1 = Bkn.*q1 - (Bk1k1.*q1 + Bk1k2.*q2);
   Bq2 = Bkn.*q2 - (Bk2k1.*q1 + Bk2k2.*q2);
   
   uF(m,1) = 1/V * sum(sum(Bq1.*e));
   uF(m,2) = 1/V * sum(sum(Bq2.*e));
    
end

uF = real(uF);


