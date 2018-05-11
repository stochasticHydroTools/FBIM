function Wk = generateWk(M)


n = M/2-1;
Wk = zeros(M,M);

s = sqrt(1/2);

%%
k1 = 1:M/2-1; k2 = 1:M/2-1; % mode
j1 = k1+1; j2 = k2+1; % index of mode
a = randn(n,n)*s; 
b = randn(n,n)*s;
Wk(j1,j2) = a+b*1i;

%%
k1 =1:M/2-1; k2 = -M/2+1:-1;
j1 = k1+1; j2 = k2+M+1;
a = randn(n,n)*s; 
b = randn(n,n)*s;
Wk(j1,j2) = a+b*1i;

%% 
k1 = 0; k2 = 1:M/2-1;
j1 = k1+1; j2 = k2+1;
a = randn(1,n)*s; 
b = randn(1,n)*s;
Wk(j1,j2) = a+b*1i;

%%
k1 = -M/2; k2 = 1:M/2-1;
j1 = k1+M+1; j2 = k2+1;
a = randn(1,n)*s; 
b = randn(1,n)*s;
Wk(j1,j2) = a+b*1i;

%%
k1= 1:M/2-1; k2 = 0;
j1 = k1+1; j2 = k2+1;
a = randn(n,1)*s; 
b = randn(n,1)*s;
Wk(j1,j2) = a+b*1i;

%%
k1 = 1:M/2-1; k2 = -M/2;
j1 = k1+1; j2 = k2+M+1;
a = randn(n,1)*s; 
b = randn(n,1)*s;
Wk(j1,j2) = a+b*1i;


%% special real-valued mode
k1 = [0,M/2]; k2 = [0,M/2];
j1 = k1+1; j2 = k2+1;
Wk(j1,j2) = randn(2,2);

%% set (0,0) mode
Wk(1,1) = 0;


% anti conj symmetry
Wk = fillanticonj(Wk,M);





