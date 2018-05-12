%% First, run BD_EM_onedisk_precompute

seed = rng('shuffle');
kBT = 1;
eta = 1;

% number of trajectories to generate
N_Traj = 1; 

%%
%-------------------------------------------------------------------------
%  LOOP OVER NUMBER OF TRAJECTORIES
%-------------------------------------------------------------------------
for ntraj = 1 : N_Traj

    % random initial position and orientation
    theta = rand*2*pi;
    x0 = rand*body_ref.L;
    y0 = rand*body_ref.L;
    nb = 1; % single body

    fID = fopen(geom_input,'w');
    fprintf(fID, '%d %d %.12f \n', 1, nb, body_ref.L);
    fprintf(fID, '%d %f %f %f %f', body_ref.np, body_ref.r, x0, y0, theta); 
    fclose(fID);

    body = load_config(geom, geom_input);

    % Stokes-Einstein relation for a disk in 2D
    XI = kBT*(4*pi*eta)^(-1)*log(body_ref.L/3.708/body_ref.r);

    % translational diffusive time scale;
    t_trans = (body_ref.r)^2 / XI;
    dt = t_trans * 1e-1;
    delta = body_ref.r * RFDtol^(1/3);


    fprintf('finish SETUP --------------------------------\n')


    %--------------------------------------------------------------------------
    % LOOP OVER TIMESTEPS 
    %--------------------------------------------------------------------------
    nrun = 1000;
    Q = zeros(nrun+1,3);
    Q(1,:) = [body.q, body.theta];

    tic;
    for n = 1:nrun

        % Store near-field contribution of M = M^{(r)} as a sparse matrix
        paramEwald = fastEwaldParameters(body.x_pos, body.L, NB , ETOL, GTOL, 'Stokeslet', Pspread);
        MStokeslet_sparse = exportStokesletRsum2d_BD(paramEwald, body);

        % Random surface vel = {M^{(r)}}^(1/2)*W^{(r)} + {M^{(w)}}^(1/2)*W^{(w)}
        vsing = sqrtMsingW(body, paramEwald, Msingalpert_ref, pMsing, pMsing_inv, MStokeslet_sparse);
        vreg  = genSqrtMreg(body.x_pos, paramEwald.Etol, paramEwald.nbox, body.L);
        vslip = real(vsing+vreg);
        vslip = reshape(vslip, [2,body.ntot])';
        vslip = sqrt(2*kBT/dt) * vslip;

        % Evaluate force
        F = zeros(1,3); 

        % Deterministic translational and angular vel: U and omega
        [NF, res_1st, t_GMRES] = ...
                 computeNF(body, F, paramEwald , Msingalpert_ref, pM1st,  vslip, MStokeslet_sparse);
        U = NF(1:2)';
        omega = NF(3);

        % RFD contribution
        Wtrans = randn(body.nb ,2) / body.L;
        Wrot   = randn(body.nb, 1);
        W = [Wtrans, Wrot];
        q_plus = body.q + (delta*body.L/2) * W(1:2);
        theta_plus = body.theta + (delta/2) * W(3);
        q_minus = body.q - (delta*body.L/2) * W(1:2);
        theta_minus = body.theta - (delta/2) * W(3);
        x_plus  = update_pos(body, q_plus, theta_plus);
        x_minus = update_pos(body, q_minus, theta_minus);

        Qplus = body; 
        Qplus.x_pos = x_plus; Qplus.q = q_plus; Qplus.theta = theta_plus;
        Qminus = body;
        Qminus.x_pos = x_minus; Qminus.q = q_minus; Qminus.theta = theta_minus;


        % Solve two additional mobility problems due to RFD piece
        paramEwald_plus = fastEwaldParameters( Qplus.x_pos, Qplus.L, NB , RFDtol, RFDtol, 'Stokeslet', Pspread);
        MStokeslet_sparse_plus = exportStokesletRsum2d_BD(paramEwald_plus, Qplus);
        [NWplus, res_1st, t_GMRES] = ...
                 computeNF(Qplus, W, paramEwald_plus , Msingalpert_refRFD,  pM1stRFD, MStokeslet_sparse_plus);

        paramEwald_minus = fastEwaldParameters( Qminus.x_pos, Qminus.L, NB , RFDtol, RFDtol, 'Stokeslet', Pspread);
        MStokeslet_sparse_minus = exportStokesletRsum2d_BD(paramEwald_minus, Qminus);
        [NWminus, res_1st, t_GMRES] = ...
                 computeNF(Qminus, W, paramEwald_minus , Msingalpert_refRFD, pM1stRFD,  MStokeslet_sparse_minus);


        DivN_U = [ NWplus(1)-NWminus(1), NWplus(2)-NWminus(2) ]; 
        DivN_omega =  NWplus(3)-NWminus(3) ;

        % Update poistion and orientation
        new_q     =  body.q     + U*dt  + dt*kBT/delta*DivN_U ;
        new_theta =  body.theta + omega*dt+ dt*kBT/delta*DivN_omega ;

        Q(n+1,:) = [new_q, new_theta];

        xnew = update_pos(body_ref, new_q , new_theta);

        body.q = new_q;
        body.theta = new_theta;
        body.x_pos = xnew;

        % print msg & plot
        if mod(n,100) == 0

            fprintf('n = %d \n', n);
            fprintf('numGMRES = %d \n', length(res_1st));

        end

        plot(body.x_pos(:,1), body.x_pos(:,2), '-'); hold on;
        hold off;
        axis equal
        axis(3*[-body.L,body.L, -body.L, body.L])
        drawnow;

    end

    toc

    %rid = randi([10000,99999]);
    %save(sprintf('ellipse%02d_EM_nrun%.0e_id%d.mat',body.nb,nrun,rid),'Q','kq','ktheta', 'dt','delta', 'seed');
    % save(sprintf('%s/%s%02dfree_EM-RFD_dt%d_nrun%.0e_id%d.mat',output_dir, geom,body.nb,th,nrun,rid), ...
    %     'Q','body_ref', 'Mobility','t_trans','t_rot',...
    %     'dt','delta','seed','ETOL','GTOL','RFDtol');


end  

