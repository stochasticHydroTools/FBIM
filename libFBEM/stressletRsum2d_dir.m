function uR = stressletRsum2d_dir(x_tar, N_tar, x_src, q_src, nv_src, ...
                                  N_src, L, xi, pinf)

uR = zeros(N_tar,2);

for m = 1 : N_tar
    
    xm = x_tar(m,:);
    
    
    for n = 1 : N_src
        
        xn = x_src(n,:); 
        qn = q_src(n,:);
        nvn= nv_src(n,:);
        
        
        for p1 = -pinf : pinf
        for p2 = -pinf : pinf
            
            p = [p1,p2]*L;
            rx = xm - (xn + p);
            r  = norm(rx);
            r2 = r^2;
            r4 = r2^2;
            xir2 = (xi*r)^2;
            D = exp(-xir2);
            
            if r > 0
                
                qnv = sum(qn.*nvn);
                rq  = sum(rx.*qn);
                rnv = sum(rx.*nvn);

                
                uR(m,1) = uR(m,1) + ((2*xi^2)*( qnv*rx(1) + rq*nvn(1) ) -...
                          4*rx(1)/r4*rq*rnv*(1+xir2)) * D;
                     
                uR(m,2) = uR(m,2) + ((2*xi^2)*( qnv*rx(2) + rq*nvn(2) ) -...
                          4*rx(2)/r4*rq*rnv*(1+xir2)) * D;
             
            end
           
        end
        end
         
    end
    
end