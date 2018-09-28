function [PrecondMsing,PrecondMsing_inv,PrecondM1st, Msingalpert_ref] = precondFBEM(paramIO, paramEwald, body_config, AlpertOrder)


Msingtrap_ref = exportMsing_trap_ref(paramIO.geom, paramIO.ref_input, ...
                paramEwald.Etol, paramEwald.nbox, paramEwald.Pspread, paramIO.output_dir ) / body_config.ds(1);
            
sprintf('%s_MsingAlpert_B01N%03d_order%02d_xi%.2f.dat', ...
                  paramIO.ref_input, body_config.np(1), AlpertOrder, ...
                  paramEwald.xi) 
Msingalpert_ref = importdata(sprintf('%s_MsingAlpert_B01N%03d_order%02d_xi%.2f.dat', ...
                  paramIO.ref_input, body_config.np(1), AlpertOrder, ...
                  paramEwald.xi)) / body_config.ds(1);
            
Mreg_ref = exportMreg(paramIO.geom, paramIO.ref_input, paramEwald.Etol, ...
                      paramEwald.nbox, paramEwald.Pspread, paramIO.output_dir) / body_config.ds(1);
                  
                  
Msingalpert_ref=(Msingalpert_ref+Msingalpert_ref')/2;
Msing_ref = (Msingtrap_ref + Msingalpert_ref);
M1st_ref = (Msing_ref + Mreg_ref);



% preconditioner for generating Msing^{1/2}
[eV,eD] = eig(Msing_ref);
eDinv = 1./diag(eD);
for j = 1:length(diag(eD))
    if eD(j,j) < 0 || abs(eD(j,j)) < paramEwald.Etol % set negative or spurious eigenvalues to zero
        eDinv(j) = 0;
    end
end

PrecondMsing = diag(sqrt(eDinv)) * eV';
PrecondMsing_inv = eV * sqrt(eD);

% preconditioner for M1st_ref
[eV,eD] = eig(M1st_ref);
eDinv = 1./diag(eD);
for j = 1:length(diag(eD))
    if eD(j,j) < 0 || abs(eD(j,j)) < paramEwald.Etol % set negative or spurious eigenvalues to zero
        eDinv(j) = 0;
        eD(j,j) = 0; 
    end
end

PrecondM1st = real(eV * diag(eDinv) / eV);
            


            
