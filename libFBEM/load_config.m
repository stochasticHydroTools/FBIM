function config = load_config(geom, geom_dir)


if strcmp(geom,'ndisk') 

    fid = fopen(geom_dir);
    tline =  fgetl(fid);
    td = str2num(tline);
    Nbody = td(2);
    L = td(3);

    npts = zeros(1,Nbody);
    rs = zeros(1,Nbody);
    xc = zeros(Nbody,2);
    ht = zeros(1,Nbody);
    phi = zeros(Nbody,1);
    
    for nbod = 1:Nbody
        tline = fgetl(fid);
        td = str2num(tline);
            
        npts(nbod) = td(1);
        rs(nbod) = td(2);
        xc(nbod,1) = td(3);
        xc(nbod,2) = td(4);
        phi(nbod) = td(5);
        ht(nbod) = 2*pi/npts(nbod);
    end
    
    
    if nargin == 1
        fid2 = fopen(sprintf('%s_FT',geom));
        FT = zeros(Nbody,3);
        for nbod = 1:Nbody
            td2 = fgetl(fid2);
            FT(nbod,:) = str2num(td2);
            
        end
        fclose(fid2);
    end

    Ntot = sum(npts);
    x = zeros(Ntot,2);
    tau = zeros(Ntot,2);
    nv = zeros(Ntot,2);
    kappa = zeros(Ntot,1);
    wgts = zeros(Ntot,1);
    cU = zeros(1,Nbody);
    cW = zeros(1,Nbody);

    nstart=0;
    for nbod = 1:Nbody
        
        i0 = nstart + (1:npts(nbod));
        t = (0:npts(nbod)-1)' * ht(nbod);

        xx = rs(nbod)*cos(t+phi(nbod)); 
        yy = rs(nbod)*sin(t+phi(nbod));

        x(i0,:) = [xx+xc(nbod,1), yy+xc(nbod,2)];
        xp = -yy;
        yp =  xx;
        nx = yp; 
        ny = -xp;
        dsdt = sqrt(nx.^2+ny.^2);
        tau(i0,:) = [xp./dsdt, yp./dsdt];
        nv(i0,:) = [nx./dsdt, ny./dsdt];

        xpp = -xx;
        ypp = -yy;
        kappa(i0,:) = (xp.*ypp-yp.*xpp) ./ dsdt.^3;
        wgts(i0,:) = dsdt*ht(nbod);
        
        cU(nbod) = sum(dsdt*ht(nbod));
        cW(nbod) = sum((xx.^2+yy.^2).*dsdt*ht(nbod)); 

        nstart = nstart+ npts(nbod);

    end

    config = struct('q', xc, 'theta', phi, 'x_pos', x, ...
                    'r', rs, 'nv', nv, 'tv', tau, ...
                    'wgts', wgts, 'kappa', kappa, 'cU', cU, 'cW', cW, ...
                    'np', npts, 'nb', Nbody, 'ntot', sum(npts) ,'L', L, 'ds', ht);



elseif strcmp(geom,'nstarfish')
    
    
    fid = fopen(geom_dir);
    tline =  fgetl(fid);
    td = str2num(tline);
    Nbody = td(2);
    L = td(3);

    
    npts = zeros(1,Nbody);
    rs = zeros(1,Nbody);
    as = zeros(1,Nbody);
    omegas = zeros(1,Nbody);
    xc = zeros(Nbody,2);
    ht = zeros(1,Nbody);
    phi = zeros(Nbody,1);

    for nbod = 1:Nbody
        tline = fgetl(fid);
        td = str2num(tline);
        
        
        npts(nbod) = td(1);
        rs(nbod) = td(2);
        as(nbod) = td(3);
        omegas(nbod) = td(4);
        xc(nbod,1) = td(5);
        xc(nbod,2) = td(6);
        phi(nbod) = td(7);

        ht(nbod) = 2*pi/npts(nbod);
    end 
    
    if nargin == 1
        fid2 = fopen(sprintf('%s_FT',geom));
        FT = zeros(Nbody,3);
        for nbod = 1:Nbody
            td2 = fgetl(fid2);
            FT(nbod,:) = str2num(td2);
            
        end
        fclose(fid2);
    end
    
    
    Ntot = sum(npts);
    x = zeros(Ntot,2);
    tau = zeros(Ntot,2);
    nv = zeros(Ntot,2);
    kappa = zeros(Ntot,1);
    wgts = zeros(Ntot,1);
    cU = zeros(1,Nbody);
    cW = zeros(1,Nbody);
    
    nstart=0;
    for nbod = 1:Nbody
        
        i0 = nstart + (1:npts(nbod));
        t = (0:npts(nbod)-1)' * ht(nbod);
        omega=omegas(nbod);
        a=as(nbod);
        r1=rs(nbod);
        
        rho  = 1+a*cos(omega*t);
        rhop = -a*omega*sin(omega*t);
        rhopp= -a*omega^2*cos(omega*t);

        xx = r1*cos(t+phi(nbod)); yy = r1*sin(t+phi(nbod));

        x(i0,:) = [xc(nbod,1) + rho.*xx, xc(nbod,2)+rho.*yy];
        xp = rhop.*xx + rho.*(-yy);
        yp = rhop.*yy + rho.*xx;
        nx = yp;
        ny = -xp;
        dsdt = sqrt(nx.^2+ny.^2);
        tau(i0,:) = [xp./dsdt,yp./dsdt];
        nv(i0,:) = [nx./dsdt,ny./dsdt];
        xpp = rhopp.*xx - 2*rhop.*yy - rho.*xx;
        ypp = rhopp.*yy + 2*rhop.*xx - rho.*yy;

        kappa(i0,:) = (xp.*ypp-yp.*xpp) ./ dsdt.^3;
        wgts(i0,:) = dsdt*ht(nbod);
        
        cU(nbod) = sum(dsdt*ht(nbod));
        cW(nbod) = sum((xx.^2+yy.^2).*dsdt*ht(nbod));        
        
        nstart = nstart + npts(nbod);
        
    end
    config = struct('q', xc, 'theta', phi, 'x_pos', x, ...
                    'r', rs, 'a', as, 'freq', omegas, 'nv', nv, 'tv', tau, ...
                    'wgts', wgts, 'kappa', kappa, 'cU', cU, 'cW', cW, ...
                    'np', npts, 'nb', Nbody, 'ntot', sum(npts) ,'L', L, 'ds', ht);


elseif strcmp(geom,'nellipse')
    
    
    fid = fopen(geom_dir);
    tline =  fgetl(fid);
    td = str2num(tline);
    Nbody = td(2);
    L = td(3);

    
    npts = zeros(1,Nbody);
    as = zeros(1,Nbody);
    bs = zeros(1,Nbody);
    xc = zeros(Nbody,2);
    ht = zeros(1,Nbody);
    phi = zeros(Nbody,1);


    for nbod = 1:Nbody
        tline = fgetl(fid);
        td = str2num(tline);
        
        npts(nbod) = td(1);
        as(nbod) = td(2);
        bs(nbod) = td(3);
        xc(nbod,1) = td(4);
        xc(nbod,2) = td(5);
        phi(nbod) = td(6);

        ht(nbod) = 2*pi/npts(nbod);
    end 
    
    if nargin == 1
        fid2 = fopen(sprintf('%s_FT',geom));
        FT = zeros(Nbody,3);
        for nbod = 1:Nbody
            td2 = fgetl(fid2);
            FT(nbod,:) = str2num(td2);
            
        end
        fclose(fid2);
    end
    
    
    Ntot = sum(npts);
    x_pos = zeros(Ntot,2);
    tau = zeros(Ntot,2);
    nv = zeros(Ntot,2);
    kappa = zeros(Ntot,1);
    wgts = zeros(Ntot,1);
    cU = zeros(1,Nbody);
    cW = zeros(1,Nbody);
    
    nstart  = 0;
    nstart2 = 0;
    for nbod = 1:Nbody
        
        i1 = nstart  + (1:npts(nbod));
        
        t = (0:npts(nbod)-1) * ht(nbod);
        
        a=as(nbod);
        b=bs(nbod);
        
        RM = [ cos(phi(nbod)), -sin(phi(nbod)) ; 
               sin(phi(nbod)),  cos(phi(nbod)) ];
           
        x  = [ a*cos(t) ; b*sin(t)];
        x = x(:);
        xp = [-a*sin(t) ; b*cos(t)];
        xp = xp(:);
        xpp= -x;  
        xpp = xpp(:);
        
        for j = 1:npts(nbod)
            x(2*j-1:2*j)   = RM * x(2*j-1:2*j);
            xp(2*j-1:2*j)  = RM * xp(2*j-1:2*j);
            xpp(2*j-1:2*j) = RM * xpp(2*j-1:2*j);
        end
        
        
        x_1 = x(1:2:end-1);
        x_2 = x(2:2:end);
        xp_1 = xp(1:2:end-1);
        xp_2 = xp(2:2:end);
        xpp_1= xpp(1:2:end-1);
        xpp_2= xpp(2:2:end);
        
        x_pos(i1,:) = [xc(nbod,1) + x_1, ...
                       xc(nbod,2) + x_2 ];
        
        dsdt = sqrt(xp_1.^2+xp_2.^2);
        tau(i1,:) = [xp_1./dsdt,  xp_2./dsdt];
        nv(i1,:)  = [xp_2./dsdt, -xp_1./dsdt];
        
        kappa(i1,:) = (xp_1.*xpp_2-xp_2.*xpp_1) ./ dsdt.^3;
        
        wgts(i1,:) = dsdt*ht(nbod);
        
        cU(nbod) = sum(dsdt*ht(nbod));
        cW(nbod) = sum((x_1.^2+x_2.^2).*dsdt*ht(nbod));        
        
        nstart = nstart + npts(nbod);
        nstart2 = nstart2 + 2*npts(nbod);
    end
    
    config = struct('q', xc, 'theta', phi, 'x_pos', x_pos, ...
                    'as', as, 'bs', bs, 'nv', nv, 'tv', tau, ...
                    'wgts', wgts, 'kappa', kappa, 'cU', cU, 'cW', cW, ...
                    'np', npts, 'nb', Nbody, 'ntot', sum(npts) ,'L', L, 'ds', ht);


    
end

fclose(fid);
