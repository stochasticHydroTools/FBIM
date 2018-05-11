%%
% qq: column vector of dimension 2*Ntot + 3*Nbody

function u= eval_nbody1stkind(qq, x, xc, npts, Nbody, Msingalpertref, L, paramEwald, StokesletRsum_sparse)

Ntot = sum(npts);
lambda = qq(1:2*Ntot);
U = qq(2*Ntot+1:end);

% apply apert quadrature to qq
using_alpert = zeros(2*Ntot,1);

nstart1= 0;
nstart2= 0;
for kbod = 1:Nbody
    
   xd = x(nstart1+1,:) - xc(kbod,:);
   [th,rr] = cart2pol(xd(1), xd(2));
      
   % 2x2 rotation matrix
   % rM = [cos(th), -sin(th); sin(th), cos(th)];
   C = cos(th); S = sin(th);
   
   
   % rotate lambda to reference config
   for j = 1:npts(kbod)
%       i2 = nstart2 + (2*j-1:2*j);
%       using_alpert(i2) = rM' * lambda(i2);
        j1 = nstart2 + (2*j-1);
        j2 = nstart2 + (2*j);
        using_alpert(j1) =  C*lambda(j1) + S*lambda(j2);
        using_alpert(j2) = -S*lambda(j1) + C*lambda(j2);

   end
    
   idx = (1:2*npts(kbod)) + nstart2;
   using_alpert(idx) = Msingalpertref * using_alpert(idx);
    
   % rotate back
   for j = 1:npts(kbod)
%       i2 = nstart2 + (2*j-1:2*j);
%       using_alpert(i2) = rM * using_alpert(i2);
        j1 = nstart2 + (2*j-1);
        j2 = nstart2 + (2*j);
        tmp1 =  C*using_alpert(j1) - S*using_alpert(j2);
        tmp2 =  S*using_alpert(j1) + C*using_alpert(j2);
        using_alpert(j1) = tmp1;
        using_alpert(j2) = tmp2;

   end
    
    nstart1 = nstart1 +   npts(kbod);
    nstart2 = nstart2 + 2*npts(kbod);
end

% fast Ewald
lambda2 = reshape(lambda, [2,Ntot])';

if nargin == 9
    uStokeslet = eval_pStokeslet2d(x, Ntot, x, Ntot, lambda2, L, paramEwald, StokesletRsum_sparse);
else
    uStokeslet = eval_pStokeslet2d(x, Ntot, x, Ntot, lambda2, L, paramEwald);
end


uStokeslet = uStokeslet'; uStokeslet = uStokeslet(:);



% M1st * lambda
M1st_lambda = using_alpert + uStokeslet;


% KU and K'lambda
KU = zeros(2*Ntot,1);
Klambda = zeros(Nbody*3,1);
nstart = 0; 
nstart2= 0;
kstart = 0;
for kbod = 1:Nbody
    idx   = (1:npts(kbod)) + nstart;
    idx_1 = (1:2:2*npts(kbod)) + nstart2;
    idx_2 = (2:2:2*npts(kbod)) + nstart2;
    
    x_perp = -(x(idx,2)-xc(kbod,2));
    y_perp =   x(idx,1)-xc(kbod,1);
    KU(idx_1) = U(1+kstart) + U(3+kstart) * x_perp;
    KU(idx_2) = U(2+kstart) + U(3+kstart) * y_perp;
    
    
    Klambda(1+kstart) = sum(lambda(idx_1));
    Klambda(2+kstart) = sum(lambda(idx_2));
    Klambda(3+kstart) = dot(x_perp,lambda(idx_1)) + dot(y_perp, lambda(idx_2));
    
    nstart = nstart + npts(kbod);
    nstart2=nstart2 + 2*npts(kbod);
    kstart = kstart + 3;
end

lambda = M1st_lambda - KU ; 
U = -Klambda;

u = [lambda ; U];



