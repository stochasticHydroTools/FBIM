function uR = rotletRsum2d(x_tar, N_tar, x_src, t_src, Boxgrid, xi, rc, NB, L, idx)

hB = L/NB;
uR = zeros(N_tar,2);
nlist = [-1,0,1];

for m = 1 : N_tar
    
   xm = mod(x_tar(m,:),L);
   
   BoxID = floor(xm/hB);
   
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
       nj = idx(j); % convert index
       
       nsc = Boxgrid(nj).N;
       ix = Boxgrid(nj).index_x;
       xn = x_src(ix,:); xn = mod(xn,L);
       Tn = t_src(ix);
       
       % for each sc in a neighbor box
       for j = 1:nsc
           rx = xm - (xn(j,:)+p);
           r  = norm(rx);
           r2 = r^2;
           
           if (r<=rc && r>0)
               
               xir2 = (xi*r)^2;
               D = exp(-xir2);
               
               uR(m,1) = uR(m,1) - rx(2)*(1-xir2)*Tn(j)*D / r2;
               uR(m,2) = uR(m,2) + rx(1)*(1-xir2)*Tn(j)*D / r2;    
           end    
       end
   end
   end
   
end