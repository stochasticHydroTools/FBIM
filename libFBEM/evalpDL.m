% -------------------------------------------------------------------------
% evalpDL.m
% 
% Bill Bao
% Nov 4, 2015
% 
% evaluate the periodic DL integral operator: (1/2+DLpv)*q at target pts
%
% qq      - density, dim = [N,2]
% xtar    - target pts
% Ntar    - number of target pts
% nv_tar  - unit normal on the targets
% tau_tar - unit tangent on the targets
% kappa   - curvature
% wgts    - weights of quadrature (trapezoidal rule in 2D)
% TOL     - tolerance parameter for Fast Ewald
% L       - domain size
% NB      - number of boxes in [0,L]
% -------------------------------------------------------------------------

function u = evalpDL(qq, xtar, Ntar, nv_tar, tau_tar, kappa, wgts,...
                    TOL, L, NB, Q, boxes, StressletRsum_sparse)
         
V = L^2;

% quadrature weighted density
qw = qq .* repmat(wgts,[1,2]);

% evaluate contribution from Stresslets using Fast Ewald, excluding self-term
if nargin == 13
    u_str = eval_pStresslet2d(xtar, Ntar, xtar, Ntar, qw, nv_tar, L, TOL,NB, Q, boxes, StressletRsum_sparse);
else
    u_str = eval_pStresslet2d(xtar, Ntar, xtar, Ntar, qw, nv_tar, L, TOL,NB, Q, boxes);
end

% non-periodic part of stresslets
qnv = sum(qw.*nv_tar,2);
u0 = [ sum(qnv.*xtar(:,1)) ,...
       sum(qnv.*xtar(:,2)) ]/V;
u0 = repmat(u0,[Ntar,1]);

% evalute self contribution analytically in 2D
tq = sum(tau_tar.*qw,2);
uD(:,1) = (2*kappa) .* tq .* tau_tar(:,1) / (4*pi);
uD(:,2) = (2*kappa) .* tq .* tau_tar(:,2) / (4*pi);

uDLpv = -(u_str + u0 + uD);

u = 0.5*qq + uDLpv;

