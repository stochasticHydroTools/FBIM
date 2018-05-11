function Mreg = exportMreg( geom, geom_dir, TOL, NB, P, output_dir )


gc = load_config(geom, geom_dir);
Ntot = sum(gc.np);
ht = 2*pi./gc.np;

%% build Mreg by pairwise Ewald code
e1 = [1,0];
e2 = [0,1];

Mreg = zeros(2*Ntot, 2*Ntot);

rwgt = 0;
fwgt = 1;

% build off-diag blocks
for m = 1:Ntot
    xm = gc.x_pos(m,:);
    for n = m : Ntot
        xn = gc.x_pos(n,:);
        if m ~= n
            [u1,xi] = eval_pStokeslet2dwgt_dir(xm,1,xn,1,e1,gc.L,TOL,NB,P,rwgt,fwgt);
            u2 = eval_pStokeslet2dwgt_dir(xm,1,xn,1,e2,gc.L,TOL,NB,P,rwgt,fwgt);
            Mreg(2*m-1:2*m,2*n-1:2*n) = [u1',u2']*ht(1);
        end
    end
end


Mreg = Mreg + Mreg';

% build diag blocks
for j = 1:Ntot
    xj = gc.x_pos(j,:);
    u1 = eval_pStokeslet2dwgt_dir(xj,1,xj,1,e1,gc.L,TOL,NB,P,rwgt,fwgt) ;
    u2 = eval_pStokeslet2dwgt_dir(xj,1,xj,1,e2,gc.L,TOL,NB,P,rwgt,fwgt) ;
    Mreg(2*j-1:2*j,2*j-1:2*j) = [u1',u2']*ht(1); 
    
end

save(sprintf('%s/%sMreg_B%02dN%03d_xi%.2f.mat',output_dir,geom,gc.nb,gc.np(1), xi),'Mreg');