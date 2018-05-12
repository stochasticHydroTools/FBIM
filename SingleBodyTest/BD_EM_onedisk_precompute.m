clear;
addpath('../libFBEM/')
addpath('../libMatCode/')

geom = 'ndisk';
geom_input = 'ndisk_input';
geom_input_ref = 'ndisk_input_ref';
input_dir = './';
output_dir = './';
body_ref = load_config(geom, geom_input_ref); % reference config

%%
%-------------------------------------------------------------------------
%  PRECOMPUTATION:
%  precompute M^(r)_alpert, M^(r)_trapezoidal, M^(w) for ref. config.
%------------------------------------------------------------------------

% user specified
norder = 8; % order of Alpert quadrature
Pspread = 21; % # of grid points to spread in NUFFT
ETOL = 1E-07; % Error Tol in Ewald
GTOL = 1E-07; % GMRES Tol
RFDtol=1E-07; % RFD Tol
NB = 4; % number of boxes for real space sum in each direction


% set parameters for Spectral Ewald
paramEwald_ref = fastEwaldParameters(body_ref.x_pos, body_ref.L, NB , ETOL, GTOL, 'Stokeslet', Pspread);
% set parameters for RFD
paramEwald_refRFD = fastEwaldParameters(body_ref.x_pos, body_ref.L, NB , RFDtol, RFDtol, 'Stokeslet', Pspread);
% set parameters for IO
paramIO = struct('geom',geom,'geom_input',geom_input,...
                 'ref_input', geom_input_ref, 'input_dir', input_dir, ...
                 'output_dir', output_dir);

% precompute conditioners
[pMsing, pMsing_inv, pM1st, Msingalpert_ref] = precondFBEM(paramIO, paramEwald_ref, body_ref, norder);
[pMsingRFD, pMsing_invRFD, pM1stRFD, Msingalpert_refRFD] = precondFBEM(paramIO, paramEwald_refRFD, body_ref, norder);