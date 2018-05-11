function KDmat = exportKDmatrix( geom, geom_dir, TOL, NB, output_dir )


[x,xc,rs,npts,L,nv,tau,wgts,kappa,cU,cW,Nbody] = load_config(geom, geom_dir);
Ntot = sum(npts);

idx = @(j) j(1)+(j(2)-1)*NB;
boxes = makeBox2d(x,Ntot,L,NB,idx);


%% construct DL matrix explicitly
DLmatrix = zeros(2*Ntot,2*Ntot);
for j = 1:Ntot
   
   q1 = zeros(Ntot,2);
   q2 = zeros(Ntot,2);
   q1(j,1) = 1;
   q2(j,2) = 1;
   
   u = evalpDL(q1,x,Ntot,nv,tau,kappa,wgts,TOL,L,NB,1,boxes);
   u = u';
   DLmatrix(:,2*j-1) = u(:);
   
   u = evalpDL(q2,x,Ntot,nv,tau,kappa,wgts,TOL,L,NB,1,boxes);
   u = u';
   DLmatrix(:,2*j) = u(:);
    
end

nstart=0;
Kp = zeros(2*Ntot,2*Ntot);
K = zeros(3*Nbody,2*Ntot);
for nbod = 1:Nbody
    npt = npts(nbod);
    idx = (1:npt)+nstart;
    idx2 = (1:2*npt) + 2*nstart;
    
    wt = wgts(idx);
    KT = kron(wt',eye(2,2)) / cU(nbod);
    KT = repmat(KT,[npt,1]);
    
    kw = [-(x(idx,2)-xc(nbod,2)), x(idx,1)-xc(nbod,1)];
    kw = kw'; kw = kw(:);
    w2 = [wt' ; wt']; w2=w2(:); 
    KW =  kw * (kw.*w2)' / cW(nbod);

    Kcomp = KT + KW; 
    
    Ksub(1:2,:) = kron(wt',eye(2,2))/ cU(nbod);
    Ksub(3,:) = (kw.*w2)' / cW(nbod);
    
    
    
    Kp(idx2,idx2) = Kcomp;
    K(3*nbod-2:3*nbod,idx2) = Ksub;
    
    nstart = nstart+npt;
    
end

A = DLmatrix - Kp;

KDmat = K/A;


fname = sprintf(sprintf('%s/%s_KD_B%02dN%03d.mat',output_dir,geom,Nbody,npts(1)));
save(fname,'KDmat')



% % Check Lambda
% u_stokeslet = eval_pStokeslet2d(x, Ntot, xc, Nbody, FT(:,1:2), L, TOL, NB);
% u_rotlet    = eval_pRotlet2d(x, Ntot, xc, Nbody, FT(:,3), L, TOL, NB);
% 
% % right hand side vector
% rhs = -(u_stokeslet+u_rotlet);
% rhs = rhs';
% rhs = rhs(:);
% 
% Lambda*rhs;



% 
% KT = kron(wgts',eye(2,2)) / cU;
% KT = repmat(KT,[Ntot,1]);
% 
% kw = [-(x(:,2)-xc(2)), x(:,1)-xc(1)];
% kw = kw'; kw = kw(:);
% wgts2 = [wgts' ; wgts']; wgts2=wgts2(:); 
% KW =  kw * (kw.*wgts2)' / cW;
% 
% Kp = KT + KW; 
% 
% A = DLmatrix - Kp;
% % % qdir = A\rhs;
% % % 
% % % 
% K = zeros(3,2*Nb);
% K(1:2,:) = kron(wgts',eye(2,2))/ cU;
% K(3,:) = (kw.*wgts2)' / cW;
% % 
% % %Y = K*qdir;
% % 
% Lambda = K*inv(A);
% 
% u_stokeslet = eval_pStokeslet2d(x, sum(Nb), xc, Nbody, [1,1], L, TOL, NB);
% u_rotlet    = eval_pRotlet2d(x, sum(Nb), xc, Nbody, 1, L, TOL, NB);
% 
% % right hand side vector
% rhs = -(u_stokeslet+u_rotlet);
% rhs = rhs';
% rhs = rhs(:);
% 
% F1 = [1,0];
% Q1 = eval_pStokeslet2d(x, sum(Nb), xc, Nbody, F1, L, TOL, NB);
% Q1 = Q1';
% F2 = [0,1];
% Q2 = eval_pStokeslet2d(x, sum(Nb), xc, Nbody, F2, L, TOL, NB);
% Q2 = Q2';
% T = 1;
% Q3 = eval_pRotlet2d(x, sum(Nb), xc, Nbody, T, L, TOL, NB);
% Q3 = Q3';
% 
% Q = [Q1(:),Q2(:),Q3(:)];
% 
%  fname = sprintf('../SecondKindform/%s_Lambda_B%d.txt',inputgeometry,Nb);
%  save(fname,'Lambda','-ascii', '-double')
%  
%  fname = sprintf('../SecondKindform/%s_Q_B%d.txt',inputgeometry,Nb);
%  save(fname,'Q','-ascii', '-double')
% 
% %  fname = sprintf('%s_Lambda_B%d.txt',inputf,Nb);
% %  save(fname,'Lambda','-ascii', '-double')
