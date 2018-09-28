function A = rotletRsum2d_dir(x_tar, N_tar, x_src, t_src, N_src, L, xi, pinf)

A = zeros(N_tar,2);

for m = 1 : N_tar
    
    xm = x_tar(m,:);
    
    for n = 1 : N_src
       
        xn = x_src(n,:);
        tn = t_src(n);
        
        % periodic images
        
        for p1 = -pinf : pinf
        for p2 = -pinf : pinf
            
            p  = [p1,p2]*L;
            rx = xm - (xn+p);
            r  = norm(rx);
            r2 = r^2;
            xir2 = (xi*r)^2;
            
            
            D = exp(-xir2);
            
            
            if r > 0
               A(m,1) = A(m,1) - rx(2)*(1-xir2)*tn*D/r2;
               A(m,2) = A(m,2) + rx(1)*(1-xir2)*tn*D/r2;
            end
            
            
            
        end
        end
        
        
        
    end

end