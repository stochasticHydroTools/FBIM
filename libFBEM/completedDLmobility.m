% -------------------------------------------------------------------------
% completedDLmobility.m
%
% Bill Bao
% Nov 4, 2015
% 
% matrix-vector multiplication routine for the mobility problem
% see Manas' TractionAnalysis.pdf for mathematical details
%
% q - density, dim = [2*N,1]
% xtar    - target pts
% Ntar    - number of target pts
% nv_tar  - unit normal on the targets
% tau_tar - unit tangent on the targets
% kappa   - curvature
% wgts    - weights of quadrature (trapezoidal rule in 2D)
% R       - a list of radii of bodies (assuming circles for now)
% xc      - centers of each body
% TOL     - tolerance parameter for Fast Ewald
% L       - domain size
% NB      - number of boxes in [0,L]
% -------------------------------------------------------------------------


function u = completedDLmobility(q, xtar, Ntar, nv_tar, tau_tar, kappa, wgts,...
                                 xc, TOL, L, NB, boxes,cU,cW, StressletRsum_sparse)
% number of bodies
Nbody = length(Ntar);

% total boundary pts
Ntot = sum(Ntar);
                             
qq = [q(1:2:end), q(2:2:end)];
qw = qq .* repmat(wgts,[1,2]);

% Q = 0;
% for j = 1:Nbody
% 
%     if j == 1
%         Nstart = 0;
%     else
%         Nstart = sum(Ntar(1:j-1));
%     end
% 
%     Nj = Ntar(j);
%     i0 = Nstart + (1:Nj);
%     qt = sum( (qq(i0,:)).^2 .* repmat(wgts(i0,:),[1,2]));
%     Q = Q + sum(qt);
% end 
% S = sum(2*pi*R);
% Q = Q*S;
Q = 1;

u = zeros(Ntot,2);


% evaluate (0.5+DLpv)*q
if nargin == 15
    uDL = evalpDL(qq,xtar,Ntot,nv_tar,tau_tar,kappa,wgts,TOL,L,NB,Q,boxes,StressletRsum_sparse);
else
    uDL = evalpDL(qq,xtar,Ntot,nv_tar,tau_tar,kappa,wgts,TOL,L,NB,Q,boxes);
end

% loop through each body, compute body velocity & angular velocity using q
for j = 1 : Nbody 
    
    if j == 1
        Nstart = 0;
    else
        Nstart = sum(Ntar(1:j-1)); 
    end
    
    Nj = Ntar(j); % number of pts on the jth body
    idx = Nstart+1:Nstart+Nj; % index in the large vector
    xj = xtar(idx,:); % point coordinates on the jth body 
    xj0 = xc(j,:); % center of the jth body
    xperp = [-xj(:,2)+xj0(2), xj(:,1)-xj0(1)];
    
    % compute velocity and angular velocity of jth body from density q
    U = 1/cU(j) * sum( qw(idx,:)  ); 
    W = 1/cW(j) * sum( qw(idx,1).*xperp(:,1)+qw(idx,2).*xperp(:,2));
    
    % matrix*vector = DL - Mobility, see notes
    u(idx,1) = uDL(idx,1) - (U(1) + W*xperp(:,1));
    u(idx,2) = uDL(idx,2) - (U(2) + W*xperp(:,2));

end

u = u';
u = u(:);



         
    