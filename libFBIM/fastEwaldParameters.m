function param = fastEwaldParameters(x_src,L,NB,eTOL,gTOL,kernalname,P)

switch kernalname
    
    case 'Stokeslet'
        cR = 100;
        cF = 1;
        
    case 'Rotlet'
        cR = 1000;
        cF = 10;
        
    case 'Stresslet'
        cR = 100;
        cF = 100;

end



% Ewald splitting parameter
xi = (NB/L) * sqrt(log(1/eTOL)+ log(cR));

% real-pace cut-off radius
rc = L/NB;

% boxing the physical domain
N_src = size(x_src,1);
idx = @(j) j(1)+(j(2)-1)*NB;
boxes = makeBox2d(x_src,N_src,L,NB,idx);



% k-space Fourier grid size
Mmax = (xi*L/pi) * sqrt(log(1/eTOL)+ log(cF));
Mmax = round(Mmax);
M  = 2*Mmax;

% gridding parameter
%P = 11;

% shape parameter of the 2nd Gaussian split in k-space
m = sqrt(P*pi);

param = struct('xi', xi, 'r_cutoff',rc, 'pts_idx',idx, 'box_arr', boxes,...
               'Mfourier', M, 'Pspread', P, 'm_shape', m, 'nbox', NB, 'Etol', eTOL, 'Gtol', gTOL);



