function [using,using_iter_err] = sqrtMsingW(paramComp, paramEwald, Msingalpert_ref, pMsing, pMsing_inv, MStokesletRsum_sparse)


x = paramComp.x_pos;
xc = paramComp.q;
npts = paramComp.np;
Nbody = paramComp.nb;
L = paramComp.L;
gTOL = paramEwald.Gtol;
Ntot = paramComp.ntot;

v = randn(Ntot,2);
vv=v'; vv=vv(:);

if nargin < 6
    fun = @(z) eval_nbodyMsing(z, x, xc, npts, Nbody,Msingalpert_ref, ...
                           pMsing, L, paramEwald);
else
    fun = @(z) eval_nbodyMsing(z, x, xc, npts, Nbody,Msingalpert_ref, ...
                           pMsing, L, paramEwald, MStokesletRsum_sparse);
end
                                              
[ysqrt,using_iter_err]=KrylovSqrtMsing(fun,vv, gTOL);

                       
% rotate and apply Lref
nstart1 = 0;
nstart2 = 0;
using = zeros(2*Ntot,1);
for kbod = 1:Nbody
   
    % compute the angle with ref frame
    xd  = x(nstart1+1,:)-xc(kbod,:);
    [th,rr] = cart2pol(xd(1),xd(2));
        
    C = cos(th); S = sin(th);
    % rotation matrix
    % rM = [cos(th), -sin(th); sin(th), cos(th)];
    %rotMat(2*kbod-1:2*kbod,:) = rM;
    
    % rotate q
    for j = 1:npts(kbod)
%        i2 = nstart2 + (2*j-1:2*j);
%        using(i2) = rM' * ysqrt(i2);
        j1 = nstart2 + (2*j-1);
        j2 = nstart2 + (2*j);
        using(j1) =  C*ysqrt(j1) + S*ysqrt(j2);
        using(j2) = -S*ysqrt(j1) + C*ysqrt(j2);
    end
    
    % apply preconditioner
    i1 = nstart2 + (1:2*npts(kbod));
    using(i1) =  pMsing_inv*using(i1);
    
    % rotat back
    for j = 1:npts(kbod)
%        i2 = nstart2 + (2*j-1:2*j);
%        using(i2) = rM * using(i2);
        j1 = nstart2 + (2*j-1);
        j2 = nstart2 + (2*j);
        tmp1 =  C*using(j1) - S*using(j2);
        tmp2 =  S*using(j1) + C*using(j2);
        using(1) = tmp1;
        using(2) = tmp2;    

    end
    
    nstart1 = nstart1+npts(kbod);
    nstart2 = nstart2+2*npts(kbod);
end






