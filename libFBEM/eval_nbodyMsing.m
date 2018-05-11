
%% with preconditioner: G * Msing * G^T

% function y = eval_nbodyMsing(qq, x, xc, npts, Nbody, Msingalpertref, Gref, L, NB, P, TOL,  StokesletRsum_sparse) %,xi, r_cutoff,idx,boxes,M,P,m_shape)
function y = eval_nbodyMsing(qq, x, xc, npts, Nbody, Msingalpertref, Gref, L, paramEwald,  StokesletRsum_sparse) %,xi, r_cutoff,idx,boxes,M,P,m_shape)

Ntot = length(x);

%% action of G*Msing,alpert*G' 
nstart1 = 0;
nstart2 = 0;
Gqq = zeros(2*Ntot,1);
y2 = zeros(2*Ntot,1);
% rotMat = zeros(2*Nbody,2);
for nbod = 1:Nbody
   
    % compute the angle with ref frame
    xd  = x(nstart1+1,:)-xc(nbod,:);
    [th,rr] = cart2pol(xd(1),xd(2));
     
    % rotation matrix
    % rM = [cos(th), -sin(th); sin(th), cos(th)];
    % rotMat(2*nbod-1:2*nbod,:) = rM;
    
    C = cos(th); S = sin(th); 
    % rM = [C -S ; S, C]
    
    
    % rotate q
    for j = 1:npts(nbod)
        % i2 = nstart2 + (2*j-1:2*j);
        % Gqq(i2) = rM' * qq(i2);
         j1 = nstart2 + (2*j-1);
         j2 = nstart2 + (2*j);
         Gqq(j1) =  C*qq(j1) + S*qq(j2);
         Gqq(j2) = -S*qq(j1) + C*qq(j2);
    end
    
    % apply preconditioner' to q
    i1 = nstart2 + (1:2*npts(nbod));
    Gqq(i1) =  Gref'*Gqq(i1); % Lchol'\ Gqq(i1);
       
    % apply alpert
    y2(i1) = sparse(Msingalpertref)*Gqq(i1);
    
    % apply preconditioner again
    y2(i1) =  Gref * y2(i1); % Lchol \ y2(i1);
    
    
    % rotate back
    for j = 1:npts(nbod)
        %i2 = nstart2 + (2*j-1:2*j);
        %y2(i2) = rM * y2(i2);
         j1 = nstart2 + (2*j-1);
         j2 = nstart2 + (2*j);
         tmp1 = C*y2(j1) - S*y2(j2);
         tmp2 = S*y2(j1) + C*y2(j2);
         y2(j1) = tmp1;
         y2(j2) = tmp2;
    end
    
    
    nstart1 = nstart1 + npts(nbod);
    nstart2 = nstart2 + npts(nbod)*2;
     
end

%% action of G*Msing,trap*G'

% apply G' to q
nstart1 = 0;
nstart2 = 0;
for nbod = 1:Nbody
    
    % rM = rotMat(2*nbod-1:2*nbod,:);
    % compute the angle with ref frame
    xd  = x(nstart1+1,:)-xc(nbod,:);
    [th,rr] = cart2pol(xd(1),xd(2));
    
    C = cos(th); S = sin(th); 

    
    %rotate q
    for j = 1:npts(nbod)
        j1 = nstart2 + (2*j-1);
        j2 = nstart2 + (2*j);
        Gqq(j1) =  C*qq(j1) + S*qq(j2);
        Gqq(j2) = -S*qq(j1) + C*qq(j2);
        % i2 = nstart2 + (2*j-1:2*j);
        % Gqq(i2) = rM' * qq(i2);
         
    end
    
    % apply preconditioner' to q
    i1 = nstart2 + (1:2*npts(nbod));
    Gqq(i1) = Gref'* Gqq(i1); %Lchol'\ Gqq(i1);    
    
    %rotate back
    for j = 1:npts(nbod)
        % i2 = nstart2 + (2*j-1:2*j);
        % Gqq(i2) = rM * Gqq(i2);
        j1 = nstart2 + (2*j-1);
        j2 = nstart2 + (2*j);
        tmp1 =  C*Gqq(j1) - S*Gqq(j2);
        tmp2 =  S*Gqq(j1) + C*Gqq(j2);
        Gqq(j1) = tmp1;
        Gqq(j2) = tmp2;
    end
    
    nstart1 = nstart1 + npts(nbod);
    nstart2 = nstart2 + npts(nbod)*2;
end

% apply Msing,trap
Gq = reshape(Gqq,[2,Ntot])';

if nargin == 10
    y1 = eval_pStokeslet2dwgt(x,Ntot,x,Ntot,Gq, L, paramEwald, 1,0, StokesletRsum_sparse);
else 
    y1 = eval_pStokeslet2dwgt(x,Ntot,x,Ntot,Gq, L, paramEwald, 1,0); %;, xi, r_cutoff, idx,boxes, M,P,m_shape, 1 ,0);
end

y1 = y1'; y1 = y1(:);


% apply G to y1
nstart1 = 0;
nstart2 = 0;
for nbod = 1:Nbody
    
    % rM = rotMat(2*nbod-1:2*nbod,:);
    xd  = x(nstart1+1,:)-xc(nbod,:);
    [th,rr] = cart2pol(xd(1),xd(2));
    
    C = cos(th); S = sin(th); 
    
    %rotate y1
    for j = 1:npts(nbod)
%        i2 = nstart2 + (2*j-1:2*j);
%        y1(i2) = rM' * y1(i2);
        j1 = nstart2 + (2*j-1);
        j2 = nstart2 + (2*j);
        tmp1 =  C*y1(j1) + S*y1(j2);
        tmp2 = -S*y1(j1) + C*y1(j2);
        y1(j1) = tmp1;
        y1(j2) = tmp2;
    end
    
    % apply preconditioner to y1
    i1 = nstart2 + (1:2*npts(nbod));
    y1(i1) = Gref * y1(i1); % Lchol \ y1(i1);    
    
    %rotate back
    for j = 1:npts(nbod)
        j1 = nstart2 + (2*j-1);
        j2 = nstart2 + (2*j);
        tmp1 =  C*y1(j1) - S*y1(j2);
        tmp2 =  S*y1(j1) + C*y1(j2);
        y1(j1) = tmp1;
        y1(j2) = tmp2;
        % i2 = nstart2 + (2*j-1:2*j);
        % y1(i2) = rM * y1(i2);
    end
    
    nstart1 = nstart1 + npts(nbod);
    nstart2 = nstart2 + npts(nbod)*2;
end

y = y1 + y2;
