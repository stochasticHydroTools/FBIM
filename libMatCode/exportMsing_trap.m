function Msing_trap = exportMsing_trap( geom, geom_dir, TOL, NB, output_dir )


[x,xc,rs,npts,L,nv,tau,wgts,kappa,cU,cW,Nbody] = load_config(geom, geom_dir);
Ntot = sum(npts);
ht = 2*pi./npts;


%% build Msing_trap by pairwise Ewald code
e1 = [1,0];
e2 = [0,1];

Msing_trap = zeros(2*Ntot, 2*Ntot);

rwgt = 1;
fwgt = 0;

% build off-diag blocks
for m = 1:Ntot
    xm = x(m,:);
    for n = m : Ntot
        xn = x(n,:);
        if m ~= n
            [u1,xi] = eval_pStokeslet2dwgt_dir(xm,1,xn,1,e1,L,TOL,NB,rwgt,fwgt);
            u2 = eval_pStokeslet2dwgt_dir(xm,1,xn,1,e2,L,TOL,NB,rwgt,fwgt);
            Msing_trap(2*m-1:2*m,2*n-1:2*n) = [u1',u2']*ht(1);
        end
    end
end


Msing_trap = Msing_trap + Msing_trap';

save(sprintf('%s/%sMsing_trap_B%02dN%03d_xi%.2f.mat', ...
             output_dir, geom,Nbody,npts(1), xi),'Msing_trap');
