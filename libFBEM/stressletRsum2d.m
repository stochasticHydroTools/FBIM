function uR = stressletRsum2d(x_tar, N_tar, x_src, q_src, nv_src, Boxgrid, ...
                              xi, rc, NB, L, idx)

hB = L/NB;
uR = zeros(N_tar,2);

for m = 1 : N_tar
   
    xm =  mod(x_tar(m,:),L); % maybe outside domain
    
    BoxID = floor( xm/hB );
    nlist = [-1,0,1];
    
    % for each neighbor
    for i1 = nlist
    for i2 = nlist
       
        ii = [i1,i2];
        % neighbor ID
        nID = BoxID + ii;
        % determine p vector shift
        a = (nID==[-1 -1]) + (nID==[NB NB]);
        p = a.*ii * L;
        
        % indices of neighbor box 
        j = mod(nID,NB) + 1;
        nj = idx(j); % convert

        nsc = Boxgrid(nj).N;
        ix  = Boxgrid(nj).index_x;
        xn  = x_src(ix,:); xn = mod(xn,L);
        qn  = q_src(ix,:);
        nvn = nv_src(ix,:);
        
        
        % for each src in a neighbor box
        for j = 1 : nsc
            rx = xm - (xn(j,:)+p);
            r  = norm(rx);
            r2 = r^2;
            r4 = r2^2;
            xir2 = xi^2*r2;
            D = exp(-xir2);
            
            if (r<=rc && r >0)
               
                qnv = sum(qn(j,:).*nvn(j,:));
                rq  = sum(rx.*qn(j,:));
                rnv = sum(rx.*nvn(j,:));
                
                uR(m,1) = uR(m,1) + ((2*xi^2)*( qnv*rx(1) + rq*nvn(j,1) ) -...
                          4*rx(1)/r4*rq*rnv*(1+xir2)) * D;
                      
                uR(m,2) = uR(m,2) + ((2*xi^2)*( qnv*rx(2) + rq*nvn(j,2) ) -...
                          4*rx(2)/r4*rq*rnv*(1+xir2)) * D;
                
            end
        end
    end
    end
    
end
