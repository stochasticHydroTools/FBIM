function y = blkdiag_precond_1stkind(q, pinvM, x, xc, Nbody, npts)

Ntot = sum(npts);
slip = q(1:2*Ntot);
F = -q(2*Ntot+1:end);

U = zeros(3*Nbody,1);
lambda = zeros(2*Ntot,1);
Ft = zeros(3,1);
KU = zeros(2*Ntot,1);

RK = zeros(2*npts(1),3);
e1 = [1;0];
e2 = [0;1];

nstart1 = 0;
nstart2 = 0;
kstart = 0;
for kbod = 1:Nbody
   
   xd = x(nstart1+1, : ) - xc(kbod,:);
   [th,rr] = cart2pol(xd(1), xd(2));
      
   % 2x2 rotation matrix
   % rM = [cos(th), -sin(th); sin(th), cos(th)];
   C = cos(th); S = sin(th);   

   % rotate slip to reference config
   for j = 1:npts(kbod)
      % j2 = nstart2 + (2*j-1:2*j);
      % lambda(j2) = rM' * slip(j2);
        j1 = nstart2 + (2*j-1);
        j2 = nstart2 + 2*j;
        lambda(j1) =  C*slip(j1) + S*slip(j2);
        lambda(j2) = -S*slip(j1) + C*slip(j2);
    end
   
   % solve M * lambda = slip
   kbod2 = nstart2 + (1:2*npts(kbod));
   lambda(kbod2) = pinvM * lambda(kbod2);
   
   % rotate back to original frame 
   for j = 1:npts(kbod)
      % j2 = nstart2 + (2*j-1:2*j);
      % lambda(j2) = rM * lambda(j2);
        j1 = nstart2 + (2*j-1);
        j2 = nstart2 + 2*j;
        tmp1 =  C*lambda(j1) - S*lambda(j2);
        tmp2 =  S*lambda(j1) + C*lambda(j2);
        lambda(j1) = tmp1; 
        lambda(j2) = tmp2;
 
   end
         
   % formulate xperp, yperp
   kbod1 = (1:npts(kbod)) + nstart1;
   xperp = -(x(kbod1,2) - xc(kbod,2));
   yperp =  (x(kbod1,1) - xc(kbod,1));
   
   
   % apply K transpose to lambda
   i1 = (1:2:2*npts(kbod));
   i2 = (2:2:2*npts(kbod));
   Ft(1) = sum(lambda(i1+nstart2));
   Ft(2) = sum(lambda(i2+nstart2));
   Ft(3) = dot(lambda(i1+nstart2), xperp) + dot(lambda(i2+nstart2), yperp);
   
   % apply rotation to K to get RK
   for j = 1:npts(kbod)
      % j1 = 2*j-1 : 2*j;
      % RK(j1,1) = rM' * e1;
      % RK(j1,2) = rM' * e2;
      % RK(j1,3) = rM' * [xperp(j) ; yperp(j)];
        j1 = 2*j-1;
        j2 = 2*j;
        RK(j1,1) =  C*e1(1) + S*e1(2);
        RK(j2,1) = -S*e1(1) + C*e1(2);
      
        RK(j1,2) =  C*e2(1) + S*e2(2);
        RK(j2,2) = -S*e2(1) + C*e2(2);
        
        RK(j1,3) =  C*xperp(j) + S*yperp(j);
        RK(j2,3) = -S*xperp(j) + C*yperp(j);

   end

   N1st = RK' * pinvM * RK;
      
   % solve for U where N * U = F - Ft
   k1 = (1:3) + kstart;
   U(k1) = N1st \ (F(k1) - Ft);
   
   % get K*U
   for j = 1:npts(kbod)
      KU(nstart2+2*j-1) = U(1+kstart) + U(3+kstart) * xperp(j);
      KU(nstart2+2*j  ) = U(2+kstart) + U(3+kstart) * yperp(j);
   end
   
   
   % rotate slip+KU to reference config
   for j = 1:npts(kbod)
        
        j1 = nstart2+ (2*j-1);
        j2 = nstart2+ 2*j;
        sKU1 = slip(j1) + KU(j1);
        sKU2 = slip(j2) + KU(j2);
        lambda(j1) =  C*sKU1 + S*sKU2;
        lambda(j2) = -S*sKU1 + C*sKU2;
%      j2 = nstart2 + (2*j-1:2*j);
%      lambda(j2) = rM' * (slip(j2)+KU(j2));
   end
   
   lambda(kbod2) = pinvM * lambda(kbod2);

   % rotate back to original frame 
   for j = 1:npts(kbod)
        
        j1 = nstart2 + (2*j-1);
        j2 = nstart2 + 2*j;
        tmp1 = C*lambda(j1) - S*lambda(j2);
        tmp2 = S*lambda(j1) + C*lambda(j2);
        lambda(j1) = tmp1;
        lambda(j2) = tmp2;
        
%      j2 = nstart2 + (2*j-1:2*j);
%      lambda(j2) = rM * lambda(j2);
   end
   
   nstart1 = nstart1 +   npts(kbod);
   nstart2 = nstart2 + 2*npts(kbod);
   kstart = kstart + 3;
end

y = [lambda ; U];
