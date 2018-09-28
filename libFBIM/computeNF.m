
function [UW_1st, res_1st, t_GRMES] = ...
         computeNF(Body, F, paramEwald , Msing_alpert_ref, pinvM, varargin) 
     
x = Body.x_pos;
q = Body.q;
npts = Body.np;
Nbody = Body.nb;
L = Body.L;
NB = paramEwald.nbox;
eTOL = paramEwald.Etol;
gmresTOL = paramEwald.Gtol;
ntot = Body.ntot;

FT = F'; FT = FT(:);

if nargin == 5
    v = zeros(2*ntot,1);
    rhs = - [v ; FT];
    evalMat = @(z) eval_nbody1stkind(z, x, q, npts, Nbody, Msing_alpert_ref, L, paramEwald);      
elseif nargin == 6
    if length(varargin{1}) == ntot
        vslip = varargin{1};
        v = vslip';
        v = v(:);
        rhs = - [v ; FT];
        evalMat = @(z) eval_nbody1stkind(z, x, q, npts, Nbody, Msing_alpert_ref, L, paramEwald);
    else
        StokesletRsum2d_sparse = varargin{1};
        v = zeros(2*ntot,1);
        rhs = - [v ; FT];
        evalMat = @(z) eval_nbody1stkind(z, x, q, npts, Nbody, Msing_alpert_ref, L, paramEwald, StokesletRsum2d_sparse);
    end
elseif nargin == 7
    vslip = varargin{1};
    StokesletRsum2d_sparse = varargin{2};
    v = vslip';
    v = v(:);
    rhs = - [v ; FT];
    evalMat = @(z) eval_nbody1stkind(z, x, q, npts, Nbody, Msing_alpert_ref, L, paramEwald, StokesletRsum2d_sparse);
end
    
applyP  = @(z) blkdiag_precond_1stkind(z, pinvM, x, q, Nbody, npts );



maxit = 80; gmresRestart = [];

t0 = tic;
[y1, flg1, relres1, iter1, resvec1] = gmres(evalMat, rhs, gmresRestart, gmresTOL, maxit, applyP); 
t_GRMES = toc(t0);

res_1st = resvec1 / norm(applyP(rhs));
UW_1st = y1(2*ntot+1:end);
