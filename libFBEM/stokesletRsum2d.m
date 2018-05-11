function uR = stokesletRsum2d(x_tar,N_tar,x_src,f_src,Boxgrid,xi,rc,NB,L,idx)

hB = L/NB;
uR = zeros(N_tar,2);

for m = 1 : N_tar
    
    xm = mod(x_tar(m,:),L); % maybe outside domain
    
    % box ID of target point
    BoxID = floor(xm/hB);
    nlist = [-1,0,1];
    
    % interaction with src in neighbour boxes
    for i1 = nlist
    for i2 = nlist
        
        ii = [i1,i2];
        % neighbor box ID, including itself
        nID = BoxID + ii;
        
        % shift vector p
        a = (nID==[-1 -1]) + (nID==[NB NB]);
        p = a.*ii * L;
        
        
        % indices of boxes 
        j = mod(nID,NB) + 1;
        nj = idx(j); % convert

        
        nsc = Boxgrid(nj).N;
        ix  = Boxgrid(nj).index_x;
        xn  = x_src(ix,:); xn = mod(xn,L);
        fn  = f_src(ix,:);

        % for each sc in an interaction box
        for j = 1:nsc
            
            % calculate periodic wrap-around distance, r = min(d1,d);
            rx = xm-(xn(j,:)+p);
            r = norm(rx);
            
            r2 = r^2;
            if (r <= rc && r > 0)
                % calculate sum
                xf = sum(rx .* fn(j,:));
                xir2 = (xi*r)^2;
                C = 1/2*expint_eone(xir2);
                D = exp(-xir2);
                
                uR(m,1) = uR(m,1) + C*fn(j,1) + (xf/r2*rx(1)-fn(j,1))*D;
                uR(m,2) = uR(m,2) + C*fn(j,2) + (xf/r2*rx(2)-fn(j,2))*D;
            end
        end
                
    end
    end
    
end