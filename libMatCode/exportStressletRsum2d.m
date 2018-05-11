function MStressletRsum2d_sparse = exportStressletRsum2d(geom, geom_dir, TOL, NB )


% [x,xc,rs,npts,L,nv,tau,wgts,kappa,cU,cW,Nbody] = load_config(geom, geom_dir);
gc = load_config(geom, geom_dir);
Ntot = sum(gc.np);
ht = 2*pi./gc.np;

x = gc.x_pos;
L = gc.L;
nv = gc.nv;

figure(1);clf;
plot(mod(x(:,1),L), mod(x(:,2),L),'.')



%% build Msing_trap by pairwise Ewald code
e1 = [1,0];
e2 = [0,1];

param = fastEwaldParameters(x,L,NB,TOL,TOL,'Stresslet',21);

xi=param.xi;
rc=param.r_cutoff;
idx=param.pts_idx;
boxes=param.box_arr;
M=param.Mfourier;
P=param.Pspread;
m=param.m_shape;

% preallocate 
Nalloc = round( (Ntot/NB^2) * 9 * Ntot * 4 );

x_index = zeros(1,Nalloc);
y_index = zeros(1,Nalloc);
MStresslet_v = zeros(1,Nalloc);




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
            xm = mod(x(im,:),gc.L);
            
            for i1=nlist
            for i2=nlist
                ii=[i1,i2];
                % neighbor ID
                nID = BoxID + ii;
                % determine p vector shift
                a = (nID==[-1 -1]) + (nID==[NB NB]);
                p = a.*ii * gc.L;
        
                % indices of neighbor box 
                j = mod(nID,NB) + 1;
                nj = idx(j); % convert
                
                nsc = boxes(nj).N;
                ix  = boxes(nj).index_x;
                xn  = x(ix,:); xn = mod(xn,L);
                nvn = nv(ix,:);
                
                for j = 1:nsc
                
                    in = ix(j);
                    rx = xm - (xn(j,:)+p);
                    r = norm(rx);
                    r2 = r^2;
                    r4 = r2^2;
                    xir2 = xi^2*r2;
                    D = exp(-xir2);
                    
                    
                    if (r<=rc && r>0)
                        
                        x_index(count+1) = 2*im-1;
                        y_index(count+1) = 2*in-1;
                        
                        x_index(count+2) = 2*im;
                        y_index(count+2) = 2*in-1;
                        
                        x_index(count+3) = 2*im-1;
                        y_index(count+3) = 2*in;
                        
                        x_index(count+4) = 2*im;
                        y_index(count+4) = 2*in;
                        
                        
                        
                        
                        qnv = sum(e1.*nvn(j,:));
                        rq  = sum(rx.*e1);
                        rnv = sum(rx.*nvn(j,:));
                        
                        MStresslet_v(count+1) = ((2*xi^2)*( qnv*rx(1) + rq*nvn(j,1) ) -...
                                                        4*rx(1)/r4*rq*rnv*(1+xir2)) * D;
                        MStresslet_v(count+2) = ((2*xi^2)*( qnv*rx(2) + rq*nvn(j,2) ) -...
                                                        4*rx(2)/r4*rq*rnv*(1+xir2)) * D; 
                                                    
                                                    
                        qnv = sum(e2.*nvn(j,:));
                        rq  = sum(rx.*e2);
                        rnv = sum(rx.*nvn(j,:));
                        
                        MStresslet_v(count+3) = ((2*xi^2)*( qnv*rx(1) + rq*nvn(j,1) ) -...
                                                        4*rx(1)/r4*rq*rnv*(1+xir2)) * D;
                        MStresslet_v(count+4) = ((2*xi^2)*( qnv*rx(2) + rq*nvn(j,2) ) -...
                                                        4*rx(2)/r4*rq*rnv*(1+xir2)) * D; 
                        
                        count = count + 4;
                    end
                    
                end
                
            end
            end
            
            
        end
       
        
        
        
    end
    
    
end

MStressletRsum2d_sparse = sparse(x_index(1:count),y_index(1:count),MStresslet_v(1:count));
%save(sprintf('%s/StressletRsum2d_sparse_B%03d_xi%03d_NB%d.mat', output_dir, Nbody, floor(xi),NB),'MStressletRsum2d_sparse');

