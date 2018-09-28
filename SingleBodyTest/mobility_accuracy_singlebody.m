clear;
set(0,'defaultaxesfontsize',18)
set(0,'defaultlinelinewidth',2)
addpath('../libFBEM/')
addpath('../libMatCode/')


% load input file
geom = 'ndisk';
geom_input = 'ndisk_input_ref';

% geom = 'nellipse';
% geom_input = 'nellipse_input_ref';

% geom = 'nstarfish';
% geom_input = 'nstarfish_input_ref';

input_dir = './';
output_dir = './';
body = load_config(geom, geom_input);

Nbox = 6;     % number of boxes in each dim for Ewald sum.  
ETOL = 1e-16; % error tol
GTOL = 1e-14; % gmres tol

% numbers of points per body
np_list = [16,24,32,48,64,96,128];
Pspread = [11,13,13,21,21,23, 23];

% order of Alpert quadrature
AlpertOrders = [4,8,16];

% choose a large N and compute the mobility to high precision with 2nd-kind
% formulation
Nps = 512;
f_tmp = 'tmp_input';
fid = fopen(f_tmp,'w');
if strcmp(geom, 'ndisk')
    
    fprintf(fid, sprintf("1 %d %.6f \n", body.nb, body.L ));
    fprintf(fid, sprintf("%d %.6f %.6f %.6f %.6f", Nps, ... 
            body.r, body.q(1), body.q(2), body.theta ));

elseif strcmp(geom, 'nellipse')

    fprintf(fid, sprintf("1 %d %.6f \n", body.nb, body.L ));
    fprintf(fid, sprintf("%d %.6f %.6f %.6f %.6f %.6f", Nps, ... 
            body.as, body.bs, body.q(1), body.q(2), body.theta ));

elseif strcmp(geom, 'nstarfish')
    
    fprintf(fid, sprintf("1 %d %.6f \n", body.nb, body.L ));
    fprintf(fid, sprintf("%d %.6f %.6f %d %.6f %.6f %.6f", Nps, ... 
         body.r, body.a, body.freq, body.q(1), body.q(2), body.theta ));
     
else
    error('Invalid input body shape.');
end

% compute mobility accurately using 2nd-kind formulation. 
I = eye(3,3);
for j = 1:3
    F = I(j,:);
    Mob2nd(:,j) = mobilityDL_manybody(geom, f_tmp , ETOL, Nbox*2, F);
end

% loop over different points per body
for k = 1:length(AlpertOrders)
for l = 1:length(np_list)
    
    % write to temporary input file
    fid = fopen(f_tmp,'w');    
    if strcmp(geom, 'ndisk')
    
        fprintf(fid, sprintf("1 %d %.6f \n", body.nb, body.L ));
        fprintf(fid, sprintf("%d %.6f %.6f %.6f %.6f", np_list(l), ... 
                body.r, body.q(1), body.q(2), body.theta ));

    elseif strcmp(geom, 'nellipse')

        fprintf(fid, sprintf("1 %d %.6f \n", body.nb, body.L ));
        fprintf(fid, sprintf("%d %.6f %.6f %.6f %.6f %.6f", np_list(l), ... 
                body.as, body.bs, body.q(1), body.q(2), body.theta ));

    elseif strcmp(geom, 'nstarfish')
    
        fprintf(fid, sprintf("1 %d %.6f \n", body.nb, body.L ));
        fprintf(fid, sprintf("%d %.6f %.6f %d %.6f %.6f %.6f", np_list(l), ... 
                body.r, body.a, body.freq, body.q(1), body.q(2), body.theta ));
     
    else
        error('Invalid input body shape...');
    end

    body_tmp = load_config(geom, f_tmp);

    paramEwald = fastEwaldParameters(body_tmp.x_pos, body_tmp.L, Nbox , ETOL, GTOL, 'Stokeslet', Pspread(l));
        
    system(sprintf('./export%sMsing_alpert %s %.16f %d', geom, f_tmp, paramEwald.xi, AlpertOrders(k)));
    Msingalpert = importdata(sprintf('%s_MsingAlpert_B01N%03d_order%02d_xi%.2f.dat', ...
                  f_tmp, body_tmp.np(1), AlpertOrders(k), ...
                  paramEwald.xi)) / body_tmp.ds(1); 
    Msingalpert = (Msingalpert + Msingalpert')/2;
          
              
    for j = 1 : 3
        F = I(j,:);
        [NF, res_1st, t_GMRES] = computeNF(body_tmp, F, paramEwald , Msingalpert, eye(2*body_tmp.np));
        Mob1st(:,j,l) = NF;
    end
    
    MobErr(k,l) = norm(Mob1st(:,:,l) - Mob2nd,2) / norm(Mob2nd,2); 
    

end
end

% clean temp data
system('rm -rf tmp_input_*.dat');

figure(2); clf;
semilogy(np_list, MobErr, 's-')
xlabel('number of points per body')
ylabel('$\|\mathbf{N}^{1st} -\mathbf{N}^{2nd}\|_2 / \|\mathbf{N}^{2nd}\|_2$', 'interpreter','latex')



