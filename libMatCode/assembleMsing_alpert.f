c ----------------------------------------------------------------      

      subroutine assembleMsing_alpert(k0,k,nsys,nd,xs,ys,dsdt,rnx,rny,
     1                                rkappa,h,swgt,order,xmat,xi)

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
      real *8 xi      

      integer nbod,kbod

      real *8 xs(*),ys(*),dsdt(*),rnx(*),rny(*),rkappa(*)
      real *8 h(0:k),swgt,xmat(nsys,nsys),qwt

      real *8 xx,yy,xy,rr,xir2

      complex *16 zdis,zn,zx,zw

      real *8 par1,pi
      complex *16 par2,par3
      complex *16, allocatable :: xmattmp(:,:)
      external ker
      real *8 eone


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
     2                 swgt*(rny(isource)**2-1.0d0)*qwt/(4.0d0*pi)

                       xmat(2*itarg-1,2*isource) =
     1                 xmat(2*itarg-1,2*isource) -
     2                 swgt*rnx(isource)*rny(isource)*qwt/(4.0d0*pi)

                       xmat(2*itarg,2*isource-1) =
     1                 xmat(2*itarg,2*isource-1) -
     2                 swgt*rnx(isource)*rny(isource)*qwt/(4.0d0*pi)

                       xmat(2*itarg,2*isource) =
     1                 xmat(2*itarg,2*isource) +
     2                 swgt*(rnx(isource)**2-1.0d0)*qwt/(4.0d0*pi)
                   endif

c                   if(isource.ne.itarg) then
c                       rr = (xs(isource)-xs(itarg))**2 + 
c     1                      (ys(isource)-ys(itarg))**2
c                       xir2 = xi**2*rr
c                       xx = (xs(isource)-xs(itarg))**2
c                       xy =
c     1                 (xs(isource)-xs(itarg))*(ys(isource)-ys(itarg))
c                       yy = (ys(isource)-ys(itarg))**2

c                       xmat(2*itarg-1,2*isource-1) =
c     1                 xmat(2*itarg-1,2*isource-1) -
c     2                 swgt*0.5d0*eone(xir2)*qwt/(4.0d0*pi)

c                       xmat(2*itarg-1,2*isource) =
c     1                 xmat(2*itarg-1,2*isource) +
c     2                 swgt*(xy/rr)*dexp(-xir2)*qwt/(4.0d0*pi)

c                       xmat(2*itarg,2*isource-1) = 
c     1                 xmat(2*itarg,2*isource-1) + 
c     2                 swgt*(xy/rr)*dexp(-xir2)*qwt/(4.0d0*pi)

c                       xmat(2*itarg,2*isource) =
c     1                 xmat(2*itarg,2*isource) -
c     2                 swgt*0.5d0*eone(xir2)*qwt/(4.0d0*pi)
c                   endif

c                  Contribution corresponding to -1/(4\pi) w log |t-z|
c                  corresponding to different boundary terms
c                   if(kbod.ne.nbod) then
c                       rr = dsqrt((xs(isource)-xs(itarg))**2 + 
c     1                      (ys(isource)-ys(itarg))**2)
c                       xir2 = xi**2*rr
c                       xmat(2*itarg-1,2*isource-1) = 
c     1                 xmat(2*itarg-1,2*isource-1) + 
c     2                 swgt*0.5d0*eone(xir2)*qwt/(4.0d0*pi)

c                       xmat(2*itarg,2*isource) = 
c     1                 xmat(2*itarg,2*isource) +  
c     2                 swgt*0.5d0*eone(xir2)*qwt/(4.0d0*pi)
c                   endif
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


c         call formslpmatbac(ier,xmattmp,order,xs(istart+1),
c     1                 ys(istart+1),dsdt(istart+1),h(nbod),nd(nbod),
c     2                 ker,xi,par2,par3)

         call formslpmatcorr(ier,xmattmp,order,xs(istart+1),
     1                 ys(istart+1),dsdt(istart+1),h(nbod),nd(nbod),
     2                 ker,xi,par2,par3)


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

      subroutine ker(par1,par2,par3,xt,yt,xs,ys,u)
         real *8 par1,xt,yt,xs,ys,pi
         real *8 rr,xi,xir2,utmp
         complex *16 zdis,par2,par3,u
         real *8 eone


         pi = 4.0d0*datan(1.0d0)
         rr = (xt-xs)**2 + (yt-ys)**2
         xi = par1
         xir2 = (xi**2)*rr
         utmp = (0.25d0/pi)*0.50d0*eone(xir2)
         u = dcmplx(utmp, 0.0d0)
                  
      return
      end   

c ----------------------------------------------------------------      

      subroutine assembleMsing_alpert2(k0,k,nsys,nd,xs,ys,dsdt,rnx,rny,
     1                                rkappa,h,swgt,order,xmat,xi)

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
      real *8 xi
      real *8, parameter :: eulergamma = 0.577215664901533d0

      integer nbod,kbod

      real *8 xs(*),ys(*),dsdt(*),rnx(*),rny(*),rkappa(*)
      real *8 h(0:k),swgt,xmat(nsys,nsys),qwt

      real *8 xx,yy,xy,rr,xir2

      complex *16 zdis,zn,zx,zw

      real *8 par1,pi
      complex *16 par2,par3
      complex *16, allocatable :: xmattmp(:,:)
      external ker2
      real *8 eone


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
     2                 swgt*(rny(isource)**2-1.0d0)*qwt/(4.0d0*pi) +
     3                 swgt*(-eulergamma/2.0d0 - dlog(xi))*qwt/(4.0d0*pi) 

                       xmat(2*itarg-1,2*isource) =
     1                 xmat(2*itarg-1,2*isource) -
     2                 swgt*rnx(isource)*rny(isource)*qwt/(4.0d0*pi)

                       xmat(2*itarg,2*isource-1) =
     1                 xmat(2*itarg,2*isource-1) -
     2                 swgt*rnx(isource)*rny(isource)*qwt/(4.0d0*pi)

                       xmat(2*itarg,2*isource) =
     1                 xmat(2*itarg,2*isource) +
     2                 swgt*(rnx(isource)**2-1.0d0)*qwt/(4.0d0*pi) + 
     3                 swgt*(-eulergamma/2.0d0 - dlog(xi))*qwt/(4.0d0*pi) 
                   endif

c                   if(isource.ne.itarg) then
c                       rr = (xs(isource)-xs(itarg))**2 + 
c     1                      (ys(isource)-ys(itarg))**2
c                       xir2 = xi**2*rr
c                       xx = (xs(isource)-xs(itarg))**2
c                       xy =
c     1                 (xs(isource)-xs(itarg))*(ys(isource)-ys(itarg))
c                       yy = (ys(isource)-ys(itarg))**2

c                       xmat(2*itarg-1,2*isource-1) =
c     1                 xmat(2*itarg-1,2*isource-1) -
c     2                 swgt*0.5d0*eone(xir2)*qwt/(4.0d0*pi)

c                       xmat(2*itarg-1,2*isource) =
c     1                 xmat(2*itarg-1,2*isource) +
c     2                 swgt*(xy/rr)*dexp(-xir2)*qwt/(4.0d0*pi)

c                       xmat(2*itarg,2*isource-1) = 
c     1                 xmat(2*itarg,2*isource-1) + 
c     2                 swgt*(xy/rr)*dexp(-xir2)*qwt/(4.0d0*pi)

c                       xmat(2*itarg,2*isource) =
c     1                 xmat(2*itarg,2*isource) -
c     2                 swgt*0.5d0*eone(xir2)*qwt/(4.0d0*pi)
c                   endif

c                  Contribution corresponding to -1/(4\pi) w log |t-z|
c                  corresponding to different boundary terms
c                   if(kbod.ne.nbod) then
c                       rr = dsqrt((xs(isource)-xs(itarg))**2 + 
c     1                      (ys(isource)-ys(itarg))**2)
c                       xir2 = xi**2*rr
c                       xmat(2*itarg-1,2*isource-1) = 
c     1                 xmat(2*itarg-1,2*isource-1) + 
c     2                 swgt*0.5d0*eone(xir2)*qwt/(4.0d0*pi)

c                       xmat(2*itarg,2*isource) = 
c     1                 xmat(2*itarg,2*isource) +  
c     2                 swgt*0.5d0*eone(xir2)*qwt/(4.0d0*pi)
c                   endif
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


c         call formslpmatbac(ier,xmattmp,order,xs(istart+1),
c     1                 ys(istart+1),dsdt(istart+1),h(nbod),nd(nbod),
c     2                 ker,xi,par2,par3)

         call formslpmatcorr(ier,xmattmp,order,xs(istart+1),
     1                 ys(istart+1),dsdt(istart+1),h(nbod),nd(nbod),
     2                 ker2,xi,par2,par3)


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

      subroutine ker2(par1,par2,par3,xt,yt,xs,ys,u)
         real *8 par1,xt,yt,xs,ys,pi
         complex *16 zdis,par2,par3,u

         pi = 4.0d0*datan(1.0d0)
         zdis = dcmplx(xt-xs,yt-ys)
         u = -0.25d0*log(cdabs(zdis))/pi
      return
      end
c-------------------------------------------------------------------   


