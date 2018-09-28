%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fillconj.m
% Yuanxun Bill Bao
% January, 2016
% 
% set entries of W with anti-symmetric conjugates to enforce conjugate symmetry
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function W = fillanticonj(W,M)

% even grid
if mod(M,2) == 0
    
    ii = 1:M/2-1; jj = 1:M/2-1;
    W(mod(M-ii,M)+1,mod(M-jj,M)+1) = -conj(W(ii+1,jj+1));


    ii =1:M/2-1; jj = M/2+1:M-1;
    W(mod(M-ii,M)+1,mod(M-jj,M)+1) = -conj(W(ii+1,jj+1));

    ii = 0; jj = 1:M/2-1;
    W(mod(M-ii,M)+1,mod(M-jj,M)+1) = -conj(W(ii+1,jj+1));

    ii = M/2; jj = 1:M/2-1;
    W(mod(M-ii,M)+1,mod(M-jj,M)+1) = -conj(W(ii+1,jj+1));

    ii = 1:M/2-1; jj = 0;
    W(mod(M-ii,M)+1,mod(M-jj,M)+1) = -conj(W(ii+1,jj+1));

    ii = 1:M/2-1; jj = M/2;
    W(mod(M-ii,M)+1,mod(M-jj,M)+1) = -conj(W(ii+1,jj+1));

else % odd grid
    i1 = 1:(M-1)/2; i2 = 0;
    W(mod(M-i1,M)+1, i2+1) = -conj(W(i1+1,i2+1));
    
    i1 = 0; i2 = 1:(M-1)/2;
    W(i1+1,  mod(M-i2,M)+1) = -conj(W(i1+1,i2+1));
    
    
    i1 = 1:(M-1)/2; i2 = 1:(M-1)/2;
    W(mod(M-i1,M)+1, mod(M-i2,M)+1) = -conj(W(i1+1,i2+1));
    
    i1 = 1:(M-1)/2; i2 = (M+1)/2:M-1;
    W(mod(M-i1,M)+1, mod(M-i2,M)+1) = -conj(W(i1+1,i2+1));
    
end
