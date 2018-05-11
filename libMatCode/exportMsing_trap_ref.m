function Msing_trap = exportMsing_trap_ref( geom, geom_dir, TOL, NB, P, output_dir )


gc = load_config(geom, geom_dir);

% [x,xc,rs,npts,L,nv,tau,wgts,kappa,cU,cW,Nbody] = load_config(geom, geom_dir);
Ntot = sum(gc.np);
ht = 2*pi./gc.np;


%% build Msing_trap by pairwise Ewald code
e1 = [1,0];
e2 = [0,1];

Msing_trap = zeros(2*Ntot, 2*Ntot);

rwgt = 1;
fwgt = 0;

% build off-diag blocks
for m = 1:Ntot
    xm = gc.x_pos(m,:);
    for n = m : Ntot
        xn = gc.x_pos(n,:);
        if m ~= n
            tic
            [u1,xi] = eval_pStokeslet2dwgt_dir(xm,1,xn,1,e1,gc.L,TOL,NB, P,rwgt,fwgt);
                 u2 = eval_pStokeslet2dwgt_dir(xm,1,xn,1,e2,gc.L,TOL,NB, P,rwgt,fwgt);
            Msing_trap(2*m-1:2*m,2*n-1:2*n) = [u1',u2']*ht(1);
            toc
        end
    end
end


Msing_trap = Msing_trap + Msing_trap';

save(sprintf('%s/%sMsing_trap_ref_B01N%03d_xi%.2f.mat',output_dir,geom,gc.np(1),xi),'Msing_trap','xi');