program nstarfishMsing_alpert

    implicit none
    integer nmax
    parameter (nmax = 4000)
    integer :: npts,nsys
    integer :: i, j, istart, kbod, k0, k, order, nd(0:10)
    real *8 :: r1, cx, cy, h(0:10), h2=1.0d0, xi, ds
    real *8 :: a, omega, rho, rhop, rhopp, xp, xpp, yp, ypp, xx,yy,phi
    real *8 :: pi, thetai, dwgt, swgt, cwgt
    real *8 xs(nmax), ys(nmax), dsdt(nmax), rnx(nmax), rny(nmax), rkappa(nmax)
    real *8, dimension(:,:), allocatable :: SL

    character (len=100) :: filename
    character (len=3) :: nb
    character (len=2) :: norder
    character (len=6) :: xi_str
    character (len=2) :: nbody
    character (len=30) :: inputfile
    character (len=20) :: xi_input
    
    
    pi = 4*datan(1.0d0)
    
    call get_command_argument(1, inputfile)
    call get_command_argument(2, xi_input)
    call get_command_argument(3, norder)
    read(xi_input,*) xi
    read(norder,*) order


    open(2, file=inputfile )

    read(2,*) k0,k
    
    nsys = 0
    npts = 0
    istart = 0
    do kbod = 1,k
        read(2,*) nd(kbod), r1, a, omega, cx, cy, phi
        do i=1,nd(kbod)
            thetai = (i-1)*2.0d0*pi/nd(1)
            rho   = 1.0d0+a*dcos(omega*thetai)
            rhop  = -a*omega*dsin(omega*thetai)
            rhopp = -a*omega**2*dcos(omega*thetai)
            xx = r1*dcos(thetai+phi)
            yy = r1*dsin(thetai+phi)
            xs(i+istart) = rho*xx + cx
            ys(i+istart) = rho*yy + cy
            xp = rhop*xx+rho*(-yy)
            yp = rhop*yy+rho* xx
            rnx(i+istart) = yp
            rny(i+istart) = -xp
            ds = dsqrt(xp**2+yp**2)
            dsdt(i+istart)= 1.0d0
            rnx(i+istart) = rnx(i+istart) / ds
            rny(i+istart) = rny(i+istart) / ds
            xpp = rhopp*xx - 2*rhop*yy - rho*xx
            ypp = rhopp*yy + 2*rhop*xx - rho*yy
            rkappa(i+istart) =  dabs(xp*ypp-xpp*yp)/ds**3
        enddo
        h(kbod) = 2.0d0*pi/nd(kbod)
        istart = istart + nd(kbod)
        nsys = nsys + 2*nd(kbod)
        npts = npts + nd(kbod)
    enddo

    allocate(SL(nsys,nsys))
 
    dwgt = 1.0d0
    swgt = 1.0d0
    cwgt = swgt/4.0d0


    call assembleMsing_alpert(k0,k,nsys,nd,xs,ys,dsdt,rnx,rny,rkappa,h, &
                        swgt,order,SL,xi)


    ! swap signs of odd rows & columns
    !xmatSL(1:nsys-1:2,:) = -xmatSL(1:nsys-1:2,:)
    !xmatSL(:,1:nsys-1:2) = -xmatSL(:,1:nsys-1:2)

    ! convert to more natural ordering
    !do i=1,N+1
    !do j=1,N+1
    !    SL(2*i-1, 2*j-1) = xmatSL(2*i  , 2*j  )
    !    SL(2*i  , 2*j  ) = xmatSL(2*i-1, 2*j-1)
    !    SL(2*i-1, 2*j  ) = xmatSL(2*i  , 2*j-1)
    !    SL(2*i  , 2*j-1) = xmatSL(2*i-1, 2*j  )
    !end do
    !end do
    write(norder,'(I2.2)') order
    write(nb, '(I3.3)') nd(1)
    write(nbody, '(I2.2)') k
    write(xi_str, '(F6.2)') xi


    !filename = 'ndiskMsing_alpert_N'//nb//'_order'//norder
    !filename =
    !'./twobodytest/twodiskMsing_alpert_N'//nb//'_order'//norder//'_dist3'
    !filename = './testSqrtMsing/ndiskMsing_alpert_ref_N'//nb//'_order'//norder
    filename=trim(inputfile)//'_MsingAlpert_B'//nbody//'N'//nb//'_order'//norder//'_xi'//trim(adjustl(xi_str))//'.dat'
    print *, filename

    open(4, file= filename)
    do i=1,nsys
        write(4, '(10000g18.9E3)') (SL(i,j), j=1,nsys)
    end do

    deallocate(SL)



end program
