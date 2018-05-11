function MStokesletRsum2d_sparse = exportStokesletRsum2d_BD(paramEwald,paramComp )


%% build Msing_trap by pairwise Ewald code
e1 = [1,0];
e2 = [0,1];


xi = paramEwald.xi;
rc = paramEwald.r_cutoff;
boxes = paramEwald.box_arr;
NB = paramEwald.nbox;
Ntot = paramComp.ntot;
x = paramComp.x_pos;
L = paramComp.L;
idx = paramEwald.pts_idx;

% preallocate 
Nalloc = round( (Ntot/NB^2) * 9 * Ntot * 4 );

x_index = zeros(1,Nalloc);
y_index = zeros(1,Nalloc);
MStokeslet_v = zeros(1,Nalloc);


nlist = [-1,0,1];
count = 0;
for nb = 1:NB^2
    
    ntar = boxes(nb).N; % number of targets in a box
    idx_list = boxes(nb).index_x;
    BoxID = boxes(nb).id;
    if ntar  > 0
       
        for m = 1:ntar
           
            % for each target, compute interaction with neighbors
            im = idx_list(m);
            xm = mod(x(im,:),L);
            
            for i1=nlist
            for i2=nlist
                ii=[i1,i2];
                % neighbor ID
                nID = BoxID + ii;
                % determine p vector shift
                a = (nID==[-1 -1]) + (nID==[NB NB]);
                p = a.*ii * L;
        
                % indices of neighbor box 
                j = mod(nID,NB) + 1;
                nj = idx(j); % convert
                
                nsc = boxes(nj).N;
                ix  = boxes(nj).index_x;
                xn  = x(ix,:); xn = mod(xn,L);
                
                for j = 1:nsc
                
                    in = ix(j);
                    rx = xm - (xn(j,:)+p);
                    r = norm(rx);
                    r2 = r^2;
                     
                    if (r<=rc && r>0)
                        
                        xir2 = xi^2 * r2;
                        C = 1/2*expint_eone(xir2);
                        D = exp(-xir2);
                        
                        x_index(count+1) = 2*im-1;
                        y_index(count+1) = 2*in-1;
                        
                        x_index(count+2) = 2*im;
                        y_index(count+2) = 2*in-1;
                        
                        x_index(count+3) = 2*im-1;
                        y_index(count+3) = 2*in;
                        
                        x_index(count+4) = 2*im;
                        y_index(count+4) = 2*in;
                        
                        xf1 = sum(rx .* e1);
                        MStokeslet_v(count+1) = C*e1(1) + (xf1/r2*rx(1)-e1(1))*D;
                        MStokeslet_v(count+2) = C*e1(2) + (xf1/r2*rx(2)-e1(2))*D;
                                                    
                        xf2 = sum(rx .* e2);
                        
                        MStokeslet_v(count+3) = C*e2(1) + (xf2/r2*rx(1)-e2(1))*D;
                        MStokeslet_v(count+4) = C*e2(2) + (xf2/r2*rx(2)-e2(2))*D;
                        
                        count = count + 4;
                    end
                    
                end
                
            end
            end
            
            
        end
       
        
        
        
    end
    
    
end

MStokesletRsum2d_sparse = sparse(x_index(1:count),y_index(1:count),MStokeslet_v(1:count));
%save(sprintf('%s/StokesletRsum2d_sparse_B%03d_xi%03d_NB%d.mat', output_dir, Nbody, floor(xi),NB),'MStokesletRsum2d_sparse');

