      Program gausselim

c     This program tests the SL + DL representation for Stokes
c     flow for solving interior/exterior simply and multiply
c     connected problems in 2D for velocity boundary conditions
C     using Gaussian elimination.
c     The total size of the problem should be less than 4000
c     and the number of boundary components should not exceed 10.

      implicit none
      integer nmax
      parameter (nmax = 4000)
      integer i,j,k,l,nbod,kbod,istart,jstart,k0,ifw
      integer itmax,itype, nd(0:10),npts,nsys,order,ii

      integer ifbipot

      real *8 a,b,cx,cy,thetai,theta,xpos,ypos,xp,yp

      real *8 h(0:10),xs(nmax),ys(nmax),dsdt(nmax)
      real *8 rnx(nmax),rny(nmax),rkappa(nmax)

      real *8 cxsource,cysource
      real *8 xvel(nmax),yvel(nmax),xvt,yvt,vort,pressure
      real *8 bipot, bipotexact

      real *8 swgt,dwgt,cwgt,pwgt,pi,h2

      real *8 t1,t2,uinf,vinf,w

      complex *16 wdens(nmax),zk(0:10),zdis,zx,ratio,eye
      complex *16 zvel,zvort

c     The weight of the modified single layer in th rerpesentation
c     The modified single layer corresponds to
c     SL - 1/ (2 \pi) w(t) ds(t)
c     where w(t) is the density
      swgt = 1.0d0

c     The weight of the double layer in the representation
      dwgt = 1.0d0
c       dwgt = 0.0d0

c     The weigt of 1/(2 \pi ) \int w(t) ds(t)
      cwgt = swgt/4.0d0

c     The weight corresponding to the generalized 1's matrix
c     Note that the generalized 1's matrix is automatically activated
c     for interior problems. 
      pwgt = 0.0d0

c     Order of alpert
      order = 16

c     itype decides the type of boundary condition at infinity
c     If itype =0, then the solution is bounded at infinity
c     If itype = 1, then the boundary conditions for the exterior
c     problem has log growth conditions uinf = 0 and vinf = -1.0d0
c     where uinf is the x component of veloicty and vinf the y component
c     Note that itype= 1 should be used only for exterior problems
      itype = 0

      call prini(6,13)

c     Open relevant files: bhinmc contains the geometry description
c     bhout is a file for writing information to
      open(2, file='geom_input')
      open(3, file='bhout')
     
c      Read type of problem to be solved. The first line in
c      bhinmc should be in the format k0,k
c      If k0=0 -- interior problem
c      k0 =1 -- exterior problem
c      k -- number of additional geometry components in domain
      read(2,*) k0,k

      if(k0.eq.0) pwgt = 1.0d0
      pwgt = 0.0d0

c     For solving exterior problems, the representation has
c     2 additional unknowns. See the routine, stokesge
c     for more details
      nsys = 0
      if(k0.eq.1) nsys = 2

      pi = 4.0d0*datan(1.0d0)
      eye = dcmplx(0.0d0,1.0d0)
      ifw = 0
c     Construct outer bounadry discretization in shape of an ellipse
      istart = 0
      if(k0.eq.0) then
         read(2,*) nd(0),a,b,cx,cy
         h(0) = 2*pi/nd(0)
         zk(0) = dcmplx(0.0d0,0.0d0)
c        Length of axes of the ellipse         
         a = 1.0d0
         b = 1.2d0
         do i = 1,nd(0)
            thetai = (i-1)*2.0d0*pi/nd(0)
            xs(i) = a*dcos(thetai)
            ys(i) = b*dsin(thetai)
            rnx(i) = b*dcos(thetai)
            rny(i) = a*dsin(thetai)
            dsdt(i) = dsqrt(rnx(i)**2 + rny(i)**2)
            rnx(i) = rnx(i)/dsdt(i)
            rny(i) = rny(i)/dsdt(i)
            rkappa(i) = a*b/dsdt(i)**3
         enddo
         istart = nd(0)
         npts = npts + nd(0)
         nsys = nsys + 2*nd(0)
      endif

      do kbod = 1,k
         read(2,*),nd(kbod),a,b,cx,cy
         zk(kbod) = dcmplx(cx,cy)
         do i = 1,nd(kbod)
            thetai = (i-1)*2.0d0*pi/nd(kbod)
            xs(istart+i) = cx + a*dcos(thetai)
            ys(istart+i) = cy + b*dsin(thetai)
            rnx(istart+i) = -b*dcos(thetai)
            rny(istart+i) = -a*dsin(thetai)
            dsdt(istart+i) = dsqrt(rnx(istart+i)**2 + 
     1                             rny(istart+i)**2)
            rnx(istart+i) = rnx(istart+i)/dsdt(istart+i)
            rny(istart+i) = rny(istart+i)/dsdt(istart+i)
            rkappa(istart+i) = -a*b/dsdt(istart+i)**3
         enddo
         h(kbod) = 2.0d0*pi/nd(kbod)
         istart = istart + nd(kbod)
         nsys = nsys + 2*nd(kbod)
         npts = npts + nd(kbod)
      enddo
c     End of creating geomety
      print *, 'Finished plotting domain'

c     The exact solution is prescribed by a stokeslet/stresslet
c     placed outside the domain. The parameters cxsource,cysource
c     correspond to the location of this source.

      if(k0.eq.0) then
         cxsource = -3.05d0 + dreal(zk(0))
         cysource = -0.02d0 + dimag(zk(0))
      else
         cxsource = dreal(zk(1)) + 0.05d0
         cysource = dimag(zk(1)) - 0.04d0
      endif

c     Based on the type of boundary conditions, set the value
c     for uinf and vinf
      if(itype.eq.1) then
          vinf = -1.0d0
          uinf = 0.0d0
      else
          uinf = 0.0d0
          vinf = 0.0d0
      endif

c     Loop over boundaries to grab exact solution data
      istart = 0
      do nbod = k0,k
         do i = 1,nd(nbod)
            call uexact(itype,xs(istart+i),ys(istart+i),cxsource,
     1           cysource,xvel(istart+i),yvel(istart+i),vort,
     2           pressure,ifbipot,bipotexact)       
         enddo
         istart = istart + nd(nbod)
      enddo

c     Compute scaling parameter for boundary condition at
c     \infty
      h2 = 1000.00
      if(k0.eq.1) then
         do nbod = k0,k
            if(h(nbod).lt.h2) h2 = h(nbod)
         enddo
      endif
      h2 = 1.0d0

c     Call the solver
      t1 = second()
      call stokesGE(k0,k,nsys,nd,xs,ys,dsdt,rnx,rny,rkappa,
     1             h,xvel,yvel,uinf,vinf,wdens,swgt,dwgt,cwgt,
     2             pwgt,h2,order)

      t2 = second()
      call prin2('Total time to compute solution = *',t2-t1,1)
      write(6,*) 'Enter x'
      read(5,*) xp
      write(6,*) 'Enter y'
      read(5,*) yp
      zx = dcmplx(xp,yp)
      call direval(zx,zvort,zvel,k0,k,nd,nsys,xs,ys,rnx,rny,
     1            rkappa,dsdt,h,wdens,swgt,dwgt,cwgt,pwgt)
      write(6,*) 'zvel = ',-eye*zvel
      call uexact(itype,xp,yp,cxsource,cysource,xvt,yvt,vort,
     1            pressure,ifw,w) 
      write(6,*) 'xvt,yvt',xvt,yvt

      return
      end
c---------------------------------------------------------------------      
      
      subroutine uexact(itype,xpos,ypos,cx,cy,xvel,yvel,vort,
     1                  pressure,ifw,w)
c      This subroutine evaluates a particular stokes field
c      corresponding to W = x log r where r = sqrt(x**2 + y**2)
c      if itype = 1 and if itype =0, it returns the field
c      phi = 1/zdis and psi = 3/zdis where zdis is the
c      distance from the source located at (cx,cy)
c      ARGUMENTS
c      ------------------------------------------------------
c      itype     IN: integer
c                type of boundary data to be picked up
c
c      xpos      In: real *8
c                x coordinate of target location where field
c                is to be evaluated
c
c      ypos      In: real *8
c                y coordinate of target location where field
c                is to be evaluated
c
c      cx        In: real *8
c                x coordinate of source location
c
c      cy        In: real *8
c                y coordinate of source location
c
c      ifw       In: integer
c                flag to compute the biharmonic potential
c--------------------------------------------------------------
c      OUTPUT
c      xvel      Out: real *8
c                x component of the velocity at the target
c      yvel      Out: real *8
c                y component of the velocity at the target
c      w         Out: real *8
c                biharomnic potential at the target
c -------------------------------------------------------------

      implicit none

      integer itype,ifw
      real *8 xpos,ypos,cx,cy,xvel,yvel,vort,pressure,w
      complex *16 ztarg,zdis,zdis2
      complex *16 phi,phip,psi,gw,chi,bipot

      ztarg = dcmplx(xpos,ypos)
      zdis = ztarg - dcmplx(cx,cy)


c     Option 0: decaying solution at infinity
      if(itype.eq.0) then
         phi = 1.0d0/zdis
         phip = -1.0d0/(zdis*zdis)
         psi = 3.0d0/zdis
         chi = 3.0d0*log(zdis)

         gw = phi + ztarg*dconjg(phip)+dconjg(psi)
         vort = 4*dreal(phip)
         pressure = 4*dimag(phip)

         xvel = dimag(gw)
         yvel = -dreal(gw)
         if(ifw.eq.1)  w = dreal(dconjg(ztarg)*phi + chi)
      else
         yvel = -log(cdabs(ztarg)) - xpos**2/cdabs(ztarg)**2
         xvel = xpos*ypos/cdabs(ztarg)**2
         vort = -2.0d0*xpos/cdabs(ztarg)**2
         if(ifw.eq.1) w = xpos*log(cdabs(ztarg))
      endif
         

      return
      end
c ----------------------------------------------------------------      

      subroutine direval(zx,zvort,zvel,k0,k,nd,nsys,xs,ys,rnx,rny,
     1                   rkappa,dsdt,h,wdens,swgt,dwgt,cwgt,pwgt)
c     Given complex density wdens compute the velocity corresponding
c     to the SL+DL representation for Stokes equation. 
c  
c     All the integrals in the representation below are computed
c     using trapezoidal rule.
c
C     The stream function W  is represented
C     by the complex variables formula
C
C           W = Re( CONJG(Z)*PHI + XI)
C     with 
c           W_x + i W_y = PHI + Z*CONJG(PHI') + CONJG(PSI)
c     where u = W_y and v = -W_x are the x and the y components
c     of velocity respectively
c
c     where PSI = XI'. The modified Sherman-Lauricella equation
c     is based on the following representation.
c
c     \phi(z) = -swgt( \int [w(t) log (t-z) ds(t)])/ (8 \pi) 
c               -dwgt( \int [w(t) / (t-z) dt])/ (4 \pi i)
c               + cwgt (\int [w(t) ds(t) ]). (2 \pi) 
c
c     where z,t are complex variables, w(t) is a complex density,
c
c     for any variable x, let xb be dconjg(x)
c     \psi(z) = -swgt(\int [wb log (t-z) ds(t) +
c                    tb*w*ds(t)/(t-z)])/(8\pi)
c               - dwgt(\int [(wb*dt + w*dtb)/(t-z) 
c                      - tb*w(t)/(t-z)**2])/(4 \pi i)
c
c-----------------------------------------------------------------
c     INPUT ARGUMENTS
c     zx       In: complex *16
c              zx = xt + i yt. Location of target in complex
c              notation
C
C     k0,k      in: INTEGER
C
C               k0 = 0 for interior and 1 for exterior problems.
C               k  = number of obstacles.
c 
C     nd(0:k)  in: integer
C
C               ND(J) = number of points in discretization of
C                       Jth boundary component.
C
C**********************************************************************
C     For future reference we define ntot = SUM_{J=K0}^K ND(J)
C**********************************************************************
C
C     The x and y coordinates of the points defining the boundary
C     components are given by
C
C     xs(ntot),ys(NTOT) IN: REAL *8 
C
C               xs(1,...,ND(K0)) xcoordinates of K0th boundary component.
C               ys(1,...,ND(K0)) xcoordinates of K0th boundary component.
C
C               ...
C
C               xs(NTOT-ND(K)+1,...,NTOT) xcoords of Kth boundary.
C               ys(NTOT-ND(K)+1,...,NTOT) xcoords of Kth boundary.
C
C**********************************************************************
C     They are assumed to be equispaced with respect to
C     some parameter theta (see dsdtheta below).
C**********************************************************************
C
C     dsdt(NTOT) IN: REAL *8 
C
C               ds/d(theta) where s denotes arclength and theta is the 
C               parametrization of the curve in which the points are 
C               equispaced (using the same ordering as the X,Y array).
c
C     rnx(NTOT),rny(NTOT) IN: REAL *8 
C
C               x and y components of outward normal at corresponding
C               point in xs,ys array.
C
C     rkappa(NTOT) IN: REAL *8 
C
C               curvature at corresponding point in xs,ys array.
C
C
C     h(0:K)   IN: REAL *8 
C
C               weight of trapezoidal rule, that is h(J) is length of 
C               Jth boundary measured in parameter theta, divided by 
C               nd(J).
C     wdens(1:ntot)   In: complex *16
c               The density for SL + DL in the above representation
c
C     swgt      in: real *8
c               Weight of the single layer in the representation for
c               velocity
c
c     dwgt      in: real *8
c               Weight of the double layer in the representation for
c               velocity
c
c     cwgt      in: real *8
c               Weight of scaled integral of density in the
c               representation for velocity (Check documentation above)
c      
c     pwgt      in: real *8
c               Weight of generalized 1's matrix to compensate
c               for the nullspace in the representation for interior
c               flows
c-------------------------------------------------------------------
c     OUTPUT arguments
c     zvort     out: complex*16
c               zvort = vort + i pressure at the target
c 
c     zvel      out: complex *16
c               complex velocity  -v + i u at the target
c-------------------------------------------------------------------

      implicit none
      integer k0,k,nd(0:k)
      integer i,ii,istart,nbod,nsys

      real *8 xs(*),ys(*),dsdt(*),rnx(*),rny(*),rkappa(*)
      real *8 h(0:k),swgt,dwgt,cwgt,pwgt

      real *8 pi

      complex *16 zn,zx,zvel,zvort,zdis,eye,zsc,zu,wdens(*)

      zvel = dcmplx(0.0d0,0.0d0)
      zvort = dcmplx(0.0d0,0.0d0)
      eye = dcmplx(0.0d0,1.0d0)
      pi = 4.0d0*datan(1.0d0)
      istart = 0

      do nbod = k0,k
         do i=1,nd(nbod)
            ii = i+istart
            zn = dcmplx(rnx(ii),rny(ii))
            zdis = dcmplx(xs(ii),ys(ii)) - zx
            zu = wdens(ii)
            zsc = h(nbod)*dsdt(ii)/(2.0d0*pi*eye)

            zvort = zvort - dwgt*0.5d0*zsc*eye*zn*zu/zdis**2
            zvort = zvort + swgt*0.25d0*zsc*eye*zu/zdis

            zvel = zvel - dwgt*0.5d0*zsc*eye*zn*zu/zdis
            zvel = zvel + dwgt*0.5d0*zsc*eye*dconjg(zn*zu/zdis**2)*
     1                     zdis
            zvel = zvel + dwgt*zsc*dreal(eye*zn*dconjg(zu))/
     1                     dconjg(zdis)

            zvel = zvel - swgt*0.5d0*zsc*zu*log(cdabs(zdis))*eye
            zvel = zvel - swgt*0.25d0*zsc*dconjg(zu)*eye*
     1                     zdis/dconjg(zdis)

            zvel = zvel + cwgt*zsc*eye*zu

         enddo
         istart = istart + nd(nbod)
      enddo
      zvort = zvort*4

      if(k0.eq.1) then
         zvel = zvel + wdens(istart+1)
      endif

      return
      end

c------------------------------------------------------------------
