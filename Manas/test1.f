      Program test1

c     This program tests the computation of SL using two methods.
c     Option 1: Using the complex notation
c
c     Option 2: Using the real notation

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

      real *8, allocatable :: rhsc(:), rhsr(:)
      real *8, allocatable :: velc1(:), velc2(:), velr(:)
      real *8, allocatable :: xmatc(:,:),xmatr(:,:)

      real *8 erra,ra

      call prini(6,13)

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

c     Order of alpert
      order = 16

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

      nsys = 0
      npts = 0

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
         read(2,*) nd(kbod),a,b,cx,cy
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

      allocate(rhsc(nsys),rhsr(nsys))
      allocate(xmatc(nsys,nsys),xmatr(nsys,nsys))
      allocate(velc1(nsys),velc2(nsys),velr(nsys))

      do i=1,nsys
         rhsc(i) = 0.0d0
         rhsr(i) = 0.0d0
         velc1(i) = 0.0d0
         velc2(i) = 0.0d0
         velr(i) = 0.0d0
      enddo

      do i=1,nsys
         do j=1,nsys
            xmatc(i,j) = 0.0d0
            xmatr(i,j) = 0.0d0
         enddo
      enddo

c     Set the respective rhs      
      ii = 1
      do i=1,npts
         rhsc(ii) = -yvel(i)
         rhsc(ii+1) = xvel(i)

         rhsr(ii) = xvel(i)
         rhsr(ii+1) = yvel(i)
         ii = ii + 2
      enddo

c     Form the matrices
      call assemblesl(k0,k,nsys,nd,xs,ys,dsdt,rnx,rny,rkappa,h,swgt,
     1                order,xmatc)

      call assembleintw(k0,k,nsys,nd,xs,ys,dsdt,rnx,rny,rkappa,h,cwgt,
     1                  xmatc)     
     
      call assembleslr(k0,k,nsys,nd,xs,ys,dsdt,rnx,rny,rkappa,h,swgt,
     1                order,xmatr)

c     Compute velocities
      do i=1,nsys
         do j=1,nsys
            velc1(i) = velc1(i) + xmatc(i,j)*rhsc(j)
            velr(i) = velr(i) + xmatr(i,j)*rhsr(j)
         enddo
      enddo
   
      ii = 1
      do i=1,npts
         velc2(ii) = velc1(ii+1)
         velc2(ii+1) = -velc1(ii)
         ii = ii + 2
      enddo
      
c     Compute error
      ra = 0.0d0
      erra = 0.0d0
      do i=1,nsys
         ra = velc2(i)**2
         erra = (velc2(i)-velr(i))**2
      enddo
      erra = dsqrt(erra/ra)

      call prin2('err=*',erra,1)

      stop
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

      subroutine assembleslr(k0,k,nsys,nd,xs,ys,dsdt,rnx,rny,rkappa,h,
     1                      swgt,order,xmat)
c     This subroutine computes the matrix entries for xmat
c     to include the single layer in the representation for velocity.
c     The rhs of the system is unpacked as
c     xvel(1),yvel(1).....xvel(n),yvel(n) and the unknowns are
c     assumed to be re(wdens(1)), im(wdens(1)....re(wdens(nsys)),
c     im(wdens(nsys)). The single layer is implemented using alpert
c     quadratures. The log part of the single layer is given by the
c     subroutine ker. ker should be defined in the following format
c
c     subroutine ker(par1,par2,par3,xt,yt,xs,ys,u)
c
c     where par1,par2,par3 are real,complex and complex parameters.
c     They are not relevant for the given kernel but are parameters
c     for the general subroutine which computes integrals using
c     alpert quadratures. 
c
c     xt,yt are the x and y coordinates of the
c     target location. 
c
c     xs and ys are x and y coordinates of the source
c
C     ARGUMENTS 
C----------------------------------------------------------------
C
C     k0,k      in: INTEGER
C
C               k0 = 0 for interior and 1 for exterior problems.
C               k  = number of obstacles.
c 
c     nsys     in: integer
c              Size of linear system to be solved using gaussian
c              elimination. nsys = 2*(\sum_{i=k0}^k nd(i)) if k0=0
c                           nsys = 2*(\sum_{i=k0}^k nd(i)+1) if k0=1
C
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
c 
c     swgt     In: real *8
c              weight of single layer in representation
c
c     order    In: integer
c              order of alpert quadratures used. The allowed values
c              for order are 0,4,8,16
c
c------------------------------------------------------------------
c    OUTPUT
c    xmat      In: real *8
c              Updated matrix to incorporate effect of single layer
c              The rhs of the system is assumed to be unpacked as
c              xvel(1),yvel(1).....xvel(n),yvel(n) and the unknowns
c              are assumed to be re(wdens(1)), im(wdens(1)...
c              re(wdens(isys),im(wdens(nsys)).
c------------------------------------------------------------------

      implicit none
      integer k0,k,nsys,nd(0:k),order,ier
      integer i,j,l,istart,jstart,isource,itarg

      integer nbod,kbod

      real *8 xs(*),ys(*),dsdt(*),rnx(*),rny(*),rkappa(*)
      real *8 h(0:k),swgt,xmat(nsys,nsys),qwt

      real *8 xx,yy,xy,rr

      complex *16 zdis,zn,zx,zw

      real *8 par1,pi
      complex *16 par2,par3
      complex *16, allocatable :: xmattmp(:,:)
      external ker


      pi = 4.0d0*datan(1.0d0)
      jstart = 0
      itarg = 0
      zx = dcmplx(0.0d0,0.0d0)
c     In the loop below, zx complex weight in the representation
c     for wb(t) and zw shall be the complex weight
c     for w(t)
      do kbod = k0,k
         do j=1,nd(kbod)
            itarg = itarg + 1
            istart = 0
            isource = 0
            do nbod = k0,k
               do i=1,nd(nbod)
                   isource = isource + 1
                   qwt = h(nbod)*dsdt(isource)
c                  Contribution of (x-y \cdot a(y)) (x-y)/|x-y|^2
                   if(isource.eq.itarg) then                   
                       xmat(2*itarg-1,2*isource-1) = 
     1                 xmat(2*itarg-1,2*isource-1) +  
     2                 swgt*rny(isource)**2*qwt/(4.0d0*pi)

                       xmat(2*itarg-1,2*isource) = 
     1                 xmat(2*itarg-1,2*isource) -
     2                 swgt*rnx(isource)*rny(isource)*qwt/(4.0d0*pi)

                       xmat(2*itarg,2*isource-1) = 
     1                 xmat(2*itarg,2*isource-1) - 
     2                 swgt*rnx(isource)*rny(isource)*qwt/(4.0d0*pi)

                       xmat(2*itarg,2*isource) = 
     1                 xmat(2*itarg,2*isource) +  
     2                 swgt*rnx(isource)**2*qwt/(4.0d0*pi)
                   endif

                   if(isource.ne.itarg) then
                       rr = (xs(isource)-xs(itarg))**2 + 
     1                      (ys(isource)-ys(itarg))**2
                       xx = (xs(isource)-xs(itarg))**2
                       xy =
     1                 (xs(isource)-xs(itarg))*(ys(isource)-ys(itarg))
                       yy = (ys(isource)-ys(itarg))**2

                       xmat(2*itarg-1,2*isource-1) = 
     1                 xmat(2*itarg-1,2*isource-1) +  
     2                 swgt*xx*qwt/(4.0d0*pi)/rr

                       xmat(2*itarg-1,2*isource) = 
     1                 xmat(2*itarg-1,2*isource) +
     2                 swgt*xy*qwt/(4.0d0*pi)/rr

                       xmat(2*itarg,2*isource-1) = 
     1                 xmat(2*itarg,2*isource-1) + 
     2                 swgt*xy*qwt/(4.0d0*pi)/rr

                       xmat(2*itarg,2*isource) = 
     1                 xmat(2*itarg,2*isource) +  
     2                 swgt*yy*qwt/(4.0d0*pi)/rr
                   endif

c                  Contribution corresponding to -1/(4\pi) w log |t-z|
c                  corresponding to different boundary terms
 
                   if(kbod.ne.nbod) then
                       rr = dsqrt((xs(isource)-xs(itarg))**2 + 
     1                      (ys(isource)-ys(itarg))**2)
                       xmat(2*itarg-1,2*isource-1) = 
     1                 xmat(2*itarg-1,2*isource-1) - 
     2                 swgt*dlog(rr)*qwt/(4.0d0*pi)

                       xmat(2*itarg,2*isource) = 
     1                 xmat(2*itarg,2*isource) -  
     2                 swgt*dlog(rr)*qwt/(4.0d0*pi)
                   endif
               enddo
               istart = istart + nd(nbod)
            enddo
         enddo
         jstart = jstart + nd(kbod)
      enddo
c     Form the matrix corresponding to -1/(4\pi) \int w log|t-z| ds(t)
c     for source and target living on the same boundary component
c     using formslpmatbac, the subroutine to compute matrix
c     entries using alpert quadrature

      istart = 0
      do nbod = k0,k
         allocate(xmattmp(nd(nbod),nd(nbod)))
         do i=1,nd(nbod)
            do j=1,nd(nbod)
               xmattmp(i,j) = dcmplx(0.0d0,0.0d0)
            enddo
         enddo
         call formslpmatbac(ier,xmattmp,order,xs(istart+1),
     1                 ys(istart+1),dsdt(istart+1),h(nbod),nd(nbod),
     2                 ker,par1,par2,par3)

         do i=1,nd(nbod)
            itarg = istart + i
            do j=1,nd(nbod)
               isource = istart + j
               xmat(2*itarg-1,2*isource-1) = 
     1         xmat(2*itarg-1,2*isource-1)+
     1         swgt*dreal(xmattmp(i,j))

               xmat(2*itarg,2*isource) = 
     1         xmat(2*itarg,2*isource)+
     1         swgt*dreal(xmattmp(i,j))
            enddo
         enddo
         istart = istart + nd(nbod)
         deallocate(xmattmp)
      enddo

      return
      end
c--------------------------------------------------------------------      

      subroutine assemblesl(k0,k,nsys,nd,xs,ys,dsdt,rnx,rny,rkappa,h,
     1                      swgt,order,xmat)
c     This subroutine updates the matrix entries for xmat
c     to include the single layer in the representation for velocity.
c     The rhs of the system is assumed to be unpacked as
c     -yvel(1),xvel(1).....-yvel(n),xvel(n) and the unknowns are
c     assumed to be -im(wdens(1)), re(wdens(1)....-im(wdens(nsys)),
c     re(wdens(nsys)). The single layer is implemented using alpert
c     quadratures. The log part of the single layer is given by the
c     subroutine ker. ker should be defined in the following format
c
c     subroutine ker(par1,par2,par3,xt,yt,xs,ys,u)
c
c     where par1,par2,par3 are real,complex and complex parameters.
c     They are not relevant for the given kernel but are parameters
c     for the general subroutine which computes integrals using
c     alpert quadratures. 
c
c     xt,yt are the x and y coordinates of the
c     target location. 
c
c     xs and ys are x and y coordinates of the source
c
c     The subroutine returns the potential at the target due to the
c     source in u (complex)
C     ARGUMENTS 
C----------------------------------------------------------------
C
C     k0,k      in: INTEGER
C
C               k0 = 0 for interior and 1 for exterior problems.
C               k  = number of obstacles.
c 
c     nsys     in: integer
c              Size of linear system to be solved using gaussian
c              elimination. nsys = 2*(\sum_{i=k0}^k nd(i)) if k0=0
c                           nsys = 2*(\sum_{i=k0}^k nd(i)+1) if k0=1
C
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
c 
c     swgt     In: real *8
c              weight of single layer in representation
c
c     order    In: integer
c              order of alpert quadratures used. The allowed values
c              for order are 0,4,8,16
c
c------------------------------------------------------------------
c    OUTPUT
c    xmat      In: real *8
c              Updated matrix to incorporate effect of single layer
c              The rhs of the system is assumed to be unpacked as
c              -yvel(1),xvel(1).....-yvel(n),xvel(n) and the unknowns
c              are assumed to be -im(wdens(1)), re(wdens(1)...
c              -im(wdens(isys),re(wdens(nsys)).
c------------------------------------------------------------------

      implicit none
      integer k0,k,nsys,nd(0:k),order,ier
      integer i,j,l,istart,jstart,isource,itarg

      integer nbod,kbod

      real *8 xs(*),ys(*),dsdt(*),rnx(*),rny(*),rkappa(*)
      real *8 h(0:k),swgt,xmat(nsys,nsys)

      complex *16 zdis,zn,zx,zw

      real *8 par1,pi
      complex *16 par2,par3, xmattmp(nsys*nsys)
      external ker


      pi = 4.0d0*datan(1.0d0)
      jstart = 0
      itarg = 0
      zx = dcmplx(0.0d0,0.0d0)
c     In the loop below, zx complex weight in the representation
c     for wb(t) and zw shall be the complex weight
c     for w(t)
      do kbod = k0,k
         do j=1,nd(kbod)
            itarg = itarg + 1
            istart = 0
            isource = 0
            do nbod = k0,k
               do i=1,nd(nbod)
                   isource = isource + 1
                   if(itarg.eq.isource) then
                      zn = dcmplx(rnx(isource),rny(isource))
c   Term corresponding to -1/(8\pi) wb (t-z)/(tb-zb) ds(t)
                      zx = 0.125d0*dsdt(isource)*h(nbod)*
     1                      zn/(dconjg(zn)*pi)
                      zw = 0.0d0
                   else
                      zdis = dcmplx(xs(isource)-xs(itarg),
     1                              ys(isource)-ys(itarg))
c   Term corresponding to -1/(8\pi) wb (t-z)/(tb-zb) ds(t)
                      zx = -0.125d0*h(nbod)*dsdt(isource)*
     1                      zdis/(dconjg(zdis)*pi)
                      zw = 0.0d0
c   Term corresponding to -1/(4\pi) w log |t-z| ds(t) for target
c   not on the same boundary component as the source
                      if(nbod.ne.kbod) then
                         zw = - 0.25d0*log(cdabs(zdis))*h(nbod)*
     1                          dsdt(isource)/pi                    
                      endif
                   endif
                   xmat(2*itarg-1,2*isource-1) = 
     1             xmat(2*itarg-1,2*isource-1) + swgt*dreal(zx+zw)

                   xmat(2*itarg,2*isource-1) = 
     1             xmat(2*itarg,2*isource-1) + swgt*dimag(zx+zw)

                   xmat(2*itarg-1,2*isource) = 
     1             xmat(2*itarg-1,2*isource) + swgt*dimag(zx-zw)

                   xmat(2*itarg,2*isource) = 
     1             xmat(2*itarg,2*isource) + swgt*dreal(zw-zx) 


               enddo
               istart = istart + nd(nbod)
            enddo
         enddo
         jstart = jstart + nd(kbod)
      enddo
c     Form the matrix corresponding to -1/(4\pi) \int w log|t-z| ds(t)
c     for source and target living on the same boundary component
c     using formslpmatbac, the subroutine to compute matrix
c     entries using alpert quadrature

      istart = 0
      do nbod = k0,k
         do i=1,nd(nbod)*nd(nbod)
             xmattmp(i) = dcmplx(0.0d0,0.0d0)
         enddo
         call formslpmatbac(ier,xmattmp,order,xs(istart+1),
     1                 ys(istart+1),dsdt(istart+1),h(nbod),nd(nbod),
     2                 ker,par1,par2,par3)

         do i=1,nd(nbod)
            itarg = istart + i
            do j=1,nd(nbod)
               isource = istart + j
               xmat(2*itarg-1,2*isource-1) = 
     1         xmat(2*itarg-1,2*isource-1)+
     1         swgt*dreal(xmattmp((j-1)*nd(nbod)+i))

               xmat(2*itarg,2*isource) = 
     1         xmat(2*itarg,2*isource)+
     1         swgt*dreal(xmattmp((j-1)*nd(nbod)+i))          
            enddo
         enddo
         istart = istart + nd(nbod)
      enddo

      return
      end

c--------------------------------------------------------------------      

      subroutine ker(par1,par2,par3,xt,yt,xs,ys,u)
         real *8 par1,xt,yt,xs,ys,pi
         complex *16 zdis,par2,par3,u

         pi = 4.0d0*datan(1.0d0)
         zdis = dcmplx(xt-xs,yt-ys)
         u = -0.25d0*log(cdabs(zdis))/pi
      return
      end
c-------------------------------------------------------------------      
      subroutine assembleintw(k0,k,nsys,nd,xs,ys,dsdt,rnx,rny,rkappa,h,
     1                      cwgt,xmat)
c     This subroutine updates the matrix entries for xmat
c     to include the integral of density in the representation for velocity.
c     The rhs of the system is assumed to be unpacked as
c     -yvel(1),xvel(1).....-yvel(n),xvel(n) and the unknowns are
c     assumed to be -im(wdens(1)), re(wdens(1)....-im(wdens(nsys)),
c     re(wdens(nsys)). The integral of density is implemented using 
c     trapezoidal rule. Note that by the intergral of density we
c     mean the following integral
c     \int_{Boundary} w(t) ds(t)/(2 \pi) and not the complex integral.
c     The boundary includes all boundary components
C     ARGUMENTS 
C----------------------------------------------------------------
C
C     k0,k      in: INTEGER
C
C               k0 = 0 for interior and 1 for exterior problems.
C               k  = number of obstacles.
c 
c     nsys     in: integer
c              Size of linear system to be solved using gaussian
c              elimination. nsys = 2*(\sum_{i=k0}^k nd(i)) if k0=0
c                           nsys = 2*(\sum_{i=k0}^k nd(i)+1) if k0=1
C
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
c 
c     cwgt     In: real *8
c              weight of integral of density
c
c------------------------------------------------------------------
c    OUTPUT
c    xmat      In: real *8
c              Updated matrix to incorporate effect of integral of
c              density. The rhs of the system is assumed to be unpacked
c              as -yvel(1),xvel(1).....-yvel(n),xvel(n) and the unknowns
c              are assumed to be -im(wdens(1)), re(wdens(1)...
c              -im(wdens(isys),re(wdens(nsys)).
c------------------------------------------------------------------

      implicit none
      integer k0,k,nsys,nd(0:k),ier
      integer i,j,l,istart,jstart,isource,itarg
      integer nbod,kbod

      real *8 xs(*),ys(*),dsdt(*),rnx(*),rny(*),rkappa(*)
      real *8 h(0:k),cwgt,xmat(nsys,nsys)

      complex *16 zdis,zn,zw,zx
      real *8 pi

      pi = 4.0d0*datan(1.0d0)

c     Target loop
      jstart = 0
      itarg = 0
      do kbod = k0,k
         do j=1,nd(kbod)
            itarg = itarg + 1
            istart = 0
            isource = 0
c           Source loop
            do nbod = k0,k
               do i=1,nd(nbod)
                  isource = isource + 1
                  zw = 0.5d0*h(nbod)*dsdt(isource)/pi
                  xmat(2*itarg-1,2*isource-1) = 
     1            xmat(2*itarg-1,2*isource-1) + cwgt*dreal(zw)

                  xmat(2*itarg,2*isource-1) = 
     1            xmat(2*itarg,2*isource-1) + cwgt*dimag(zw)

                  xmat(2*itarg-1,2*isource) = 
     1            xmat(2*itarg-1,2*isource) + cwgt*dimag(-zw)

                  xmat(2*itarg,2*isource) = 
     1            xmat(2*itarg,2*isource) + cwgt*dreal(zw) 
               enddo
               istart = istart + nd(nbod)
            enddo
         enddo
         jstart = jstart + nd(kbod)
      enddo

      return
      end
