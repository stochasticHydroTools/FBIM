function uR = stokesletRsum2d_dir(x_tar, N_tar, x_src, f_src, N_src, L, xi, pinf)

uR = zeros(N_tar,2);

for m = 1 : N_tar
    
   xm = x_tar(m,:);
   
   for n = 1 : N_src
      
       x = x_src(n,:);
       f = f_src(n,:);
       
       % periodic images
       
       for p1 = -pinf : pinf
       for p2 = -pinf : pinf
           
           p = [p1,p2]*L;
           rx = xm - (x+p);
           r  = norm(rx);
           r2 = r^2;
           rxf = sum(rx.*f);
           xir = xi * r;
           xir2= xir^2;
           
           C = 1/2*expint(xir2);
           D = exp(-xir2);
           
           if r > 0 
              uR(m,1) = uR(m,1) + C*f(1) + (rxf/r2*rx(1)-f(1))*D;
              uR(m,2) = uR(m,2) + C*f(2) + (rxf/r2*rx(2)-f(2))*D;
           end
           
       end
       end
       
   end
    
    
end