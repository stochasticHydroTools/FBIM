c-----------------------------------------------------------------
c     This file has the subroutines to generate the matrix
c     for solving Stokes equation in multiply connected interior
c     or exterior domains and solve them using Gaussian elimination
c------------------------------------------------------------------
      subroutine stokesGE(k0,k,nsys,nd,xs,ys,dsdt,rnx,rny,rkappa,
     1             h,xvel,yvel,uinf,vinf,wdens,swgt,dwgt,cwgt,
     2             pwgt,h2,order)
C
c     This subroutine solves the 
C     equation for Stokes flow in multiply connected interior
C     or exterior domains using combination of SL and DL. 
c     The type of problem is determined
C     by the flag k0. k0 = 0 implies the problem is an interior
C     one, while k0 = 1 implies the problem is an exterior one.
C     The enclosing boundary, if it exists, is numbered 0. 
C     The obstacle boundaries are numbered 1,...,k. 
C     
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
C    
c     As we take the limit of this representation we get the system 
c     of equations Aw = b. 
c 
C     For the exterior problem, we do Phi(z) = Phi(z) + d for
c     some complex constant d and using the conditions at infinity
c     we set up the block system
c            [A    B]    [w] = [f] 
c            [F    0]    [d] = [C]
c     where C is the log growth of the components of velocity
c     at infinity
c
c     The linear system is solved using Gaussian elimination
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
C
C     xvel(NTOT),yvel(NTOT)   IN: REAL *8 
C
C               x and y components of velocity at
C               corresponding boundary point.
C
C     uinf, vinf  IN: REAL*8  
C
C               UINF,VINF  determine the leading order behavior of 
C               the velocity at infinity. In particular, (U,V) approaches
C               (UINF log(r),VINF log(r)) + O(1) as r  -> \infty.  
C               For interior flows, set UINF = VINF = 0.0D0.
C
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
c
c     h2        in: real *8
c               scaling factor for implementing growth conditions
c               at infinity
c               The condition is
c               \int w(t) ds (t) / h2 = Log Growth/h2
c
c     order     in: integer
c               Order of alpert quadrature for single layer in the
c               representation. Allowed values of order are 0,4,8 and 16
c
C--------------------------------------------------------------------
C     OUTPUT VARIABLES: 
C
C     wdens(nsys) OUT: COMPLEX*16
C
C               complex density for Sherman-Lauricella representation 
C               of flow field (with ordering corresponding to input
C               X,Y array). 
C
C--------------------------------------------------------------------

      implicit none
      integer i,j,l,istart,jstart,job,nbod,ii
      integer k0,k,nsys,nd(0:k),order
      integer ifeigvec,ier

      real *8 xs(*),ys(*),dsdt(*),rnx(*),rny(*),rkappa(*)
      real *8 h(0:k),xvel(*),yvel(*),uinf,vinf
      real *8 swgt,dwgt,pwgt,cwgt

      complex *16 wdens(*), zk(0:k)

      real *8 rhs(nsys),solution(nsys),xmat(nsys,nsys)
      real *8 xmatsvd(nsys,nsys),xmatrg(nsys,nsys)
      real *8 sing(nsys),workspace(nsys*nsys)

      real *8 eigi(nsys),eigr(nsys)
      real *8 eigvec(nsys,nsys)
      real *8 eigtemp(nsys*nsys), eigtemp1(nsys*nsys)

      real *8 rsingvec(nsys,nsys),lsingvec(nsys,nsys)
      real *8 rcond,work(5*nsys)
      real *8 maxsing,minsing,maxeig,mineig

      real *8 h2
      logical lvec,rvec
      character(len=5) nsysstring
      character(len=25) filename

c     Initialize solution and matrix entries to 0
      do i=1,nsys
         solution(i) = 0.0d0
      enddo

      do i=1,nsys
         do j=1,nsys
            xmat(i,j) = 0.0d0
            xmatsvd(i,j) = 0.0d0
            xmatrg(i,j) = 0.0d0
            lsingvec(i,j) = 0.0d0
            rsingvec(i,j) = 0.0d0
            eigvec(i,j) = 0.0d0
            workspace((i-1)*nsys+j) = 0.0d0
            eigtemp((i-1)*nsys+j) = 0.0d0
            eigtemp1((i-1)*nsys+j) = 0.0d0
         enddo
         sing(i) = 0.0d0
         eigi(i) = 0.0d0
         eigr(i) = 0.0d0
      enddo
      
      job = 0
      do i =1,5*nsys
         work(i) = 0.0d0
      enddo
c     Setting up the rhs for the solve
      istart = 0
      do nbod = k0,k
         do i =1,nd(nbod)
            rhs(2*(istart+i)-1) = -yvel(istart+i)
            rhs(2*(istart+i)) = xvel(istart+i)
         enddo
         istart = istart + nd(nbod)
      enddo

      if(k0.eq.1) then
         rhs(nsys - 1) = -vinf/h2
         rhs(nsys) = uinf/h2
      endif

c     Form the matrix
c     Single layer
      if(abs(swgt).ge.1.0d-15) then
         call assemblesl(k0,k,nsys,nd,xs,ys,dsdt,rnx,rny,rkappa,h,
     1                   swgt,order,xmat)
cc         call prin2('xmat=*',xmat,12)
      endif

      if(abs(dwgt).ge.1.0d-15) then
         call assembledl(k0,k,nsys,nd,xs,ys,dsdt,rnx,rny,rkappa,h,
     1                  dwgt,xmat)
      endif

      if(abs(cwgt).ge.1.0d-15) then
        call assembleintw(k0,k,nsys,nd,xs,ys,dsdt,rnx,rny,rkappa,h,
     1                    cwgt,xmat)
      endif

      if(k0.eq.1) then
          call assembleconst(k0,k,nsys,nd,xs,ys,dsdt,rnx,rny,rkappa,
     1                       h,h2,swgt,xmat)
      endif

      if(k0.eq.0) then
          call assemblepert(k0,k,nsys,nd,xs,ys,dsdt,rnx,rny,rkappa,
     1                      h,pwgt,xmat)     
      endif

      do i=1,nsys
         do j=1,nsys
            solution(i) = solution(i) + xmat(i,j)*rhs(j) 
         enddo
         write(108,*) solution(i)
      enddo

cc    Uncommment this section of the code to
c     compute either the svd or the eigenvalue decomposition.
c     This section creates copies of the matrix xmat
c     which will then be used for the computation of the svd
c     or the eigenvalue decomposition

      do i=1,nsys
         do j=1,nsys
            xmatsvd(i,j) = xmat(i,j)
            xmatrg(i,j) = xmat(i,j)
            write(105,*) xmat(i,j)
         enddo
         write(106,*) rhs(i)
      enddo
c    End uncommenting here

c     SVD of the matrix
c     Uncomment this section to compute the svd of the
c     matrix. If the total number of points is <1000, then
c     the singular values are stored in <Numpts>sing.dat
c     To compute the right or the left singular vectors
c     allocate necessary arrays and set lvec or rvec to be
c     true

      lvec = .false.
      rvec = .false.
      call svd(nsys,nsys,nsys,xmatsvd,sing,lvec,lsingvec,rvec,
     1         rsingvec,ier,workspace)




c      maxsing = 0.0d0
c      minsing = 1.0d2
c      write(nsysstring,'(I3)') nsys
c      filename = trim(nsysstring)//'sing.dat'
c      open(unit=34,file=filename)
      do i=1,nsys 
          write(34,*) sing(i)
          if(abs(sing(i)).gt.maxsing) maxsing = abs(sing(i))
          if(abs(sing(i)).lt.minsing) minsing = abs(sing(i))
      enddo

c     End uncommenting here
c     End of svd


c     Eigenvalue decomposition of the matrix
c     Uncomment this section to compute the eig decomposition of
c     the matrix. If the total number of points is <1000, then
c     the eigenvalues values are stored in <Numpts>eig.dat
c     To compute the eigenvectors,
c     allocate necessary array and set ifeigvec 1
c     true

      ifeigvec = 0
      call rg(nsys,nsys,xmatrg,eigr,eigi,ifeigvec,eigvec,
     1       eigtemp,eigtemp1,ier)

c      maxeig = 0.0d0
c      mineig = 1.0d2
c      open(unit=36,file=trim(nsysstring)//'eig.dat')
      do i=1,nsys
         write(36,*) eigr(i)
         if(abs(eigr(i)).gt.maxeig) maxeig = abs(eigr(i))
         if(abs(eigr(i)).lt.mineig) mineig = abs(eigr(i))
      enddo
c      End uncommenting here
c      End of computing eigenvalues
 
      call dgeco(xmat,nsys,nsys,workspace,rcond,work)
      call prin2('rcond=*',rcond,1)

     

      call dgesl(xmat,nsys,nsys,workspace,rhs,job)
cc      call prin2('solution=*',rhs,24)

      do i=1,nsys/2
         wdens(i) = dcmplx(rhs(2*i-1),rhs(2*i))
      enddo

      do i=1,nsys
         write(120,*) rhs(i)
      enddo

      return
      end
c---------------------------------------------------------------      

      subroutine ker(par1,par2,par3,xt,yt,xs,ys,u)
         real *8 par1,xt,yt,xs,ys,pi
         complex *16 zdis,par2,par3,u

         pi = 4.0d0*datan(1.0d0)
         zdis = dcmplx(xt-xs,yt-ys)
         u = -0.25d0*log(cdabs(zdis))/pi
      return
      end
c-------------------------------------------------------------------      

      subroutine assemblesl(k0,k,nsys,nd,xs,ys,dsdt,rnx,rny,rkappa,h,
     1                      swgt,order,xmat)
c     This subroutine updates the matrix entries for xmat
c     to include the single layer in the representation for velocity.
c     The rhs of the system is assumed to be unpacked as
c     -yvel(1),xvel(1).....-yvel(n),xvel(n) and the unknowns are
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
c              are assumed to be re(wdens(1)), im(wdens(1)...
c              re(wdens(isys),im(wdens(nsys)).
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

      subroutine assembledl(k0,k,nsys,nd,xs,ys,dsdt,rnx,rny,rkappa,h,
     1                      dwgt,xmat)
c     This subroutine updates the matrix entries for xmat
c     to include the double layer in the representation for velocity.
c     The rhs of the system is assumed to be unpacked as
c     -yvel(1),xvel(1).....-yvel(n),xvel(n) and the unknowns are
c     assumed to be re(wdens(1)), im(wdens(1)....re(wdens(nsys)),
c     im(wdens(nsys)). The double layer is implemented using 
c     trapezoidal rule. 
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
c     dwgt     In: real *8
c              weight of double layer in representation
c
c------------------------------------------------------------------
c    OUTPUT
c    xmat      In: real *8
c              Updated matrix to incorporate effect of double layer
c              The rhs of the system is assumed to be unpacked as
c              -yvel(1),xvel(1).....-yvel(n),xvel(n) and the unknowns
c              are assumed to be re(wdens(1)), im(wdens(1)...
c              re(wdens(isys),im(wdens(nsys)).
c------------------------------------------------------------------

      implicit none
      integer k0,k,nsys,nd(0:k),ier
      integer i,j,l,istart,jstart,isource,itarg
      integer nbod,kbod

      real *8 xs(*),ys(*),dsdt(*),rnx(*),rny(*),rkappa(*)
      real *8 h(0:k),dwgt,xmat(nsys,nsys),wgt

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
                  wgt = h(nbod)*dsdt(isource)/pi
                  zn = dcmplx(rnx(isource),rny(isource))
                  if(itarg.eq.isource) then
                     zw = -0.5d0 - 0.25d0*wgt*rkappa(isource)
                     zx = -0.25d0*wgt*rkappa(isource)*zn*zn
                  else
                     zdis = dcmplx(xs(isource)-xs(itarg),
     1                             ys(isource)-ys(itarg))
                     zw = -0.5d0*wgt*dreal(zn/zdis)
                     zx = 0.5d0*wgt*dreal(zn/zdis)*(zdis/dconjg(zdis))
                  endif
                  xmat(2*itarg-1,2*isource-1) = 
     1            xmat(2*itarg-1,2*isource-1) + dwgt*dreal(zx+zw)

                  xmat(2*itarg,2*isource-1) = 
     1            xmat(2*itarg,2*isource-1) + dwgt*dimag(zx+zw)

                  xmat(2*itarg-1,2*isource) = 
     1            xmat(2*itarg-1,2*isource) + dwgt*dimag(zx-zw)

                  xmat(2*itarg,2*isource) = 
     1            xmat(2*itarg,2*isource) + dwgt*dreal(zw-zx) 
               enddo
               istart = istart + nd(nbod)
            enddo
         enddo
         jstart = jstart + nd(kbod)
      enddo

      return
      end

      subroutine assembleintw(k0,k,nsys,nd,xs,ys,dsdt,rnx,rny,rkappa,h,
     1                      cwgt,xmat)
c     This subroutine updates the matrix entries for xmat
c     to include the integral of density in the representation for velocity.
c     The rhs of the system is assumed to be unpacked as
c     -yvel(1),xvel(1).....-yvel(n),xvel(n) and the unknowns are
c     assumed to be re(wdens(1)), im(wdens(1)....re(wdens(nsys)),
c     im(wdens(nsys)). The integral of density is implemented using 
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
c              are assumed to be re(wdens(1)), im(wdens(1)...
c              re(wdens(isys),im(wdens(nsys)).
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

      subroutine assembleconst(k0,k,nsys,nd,xs,ys,dsdt,rnx,rny,rkappa,h,
     1                      h2,swgt,xmat)
c     This subroutine updates the matrix entries for xmat
c     to enforce the growth conditions at infinity and include the
c     constants in the repesentation to enforce the condition
c     The rhs of the system is assumed to be unpacked as
c     -yvel(1),xvel(1).....-yvel(n),xvel(n) and the unknowns are
c     assumed to be re(wdens(1)), im(wdens(1)....re(wdens(nsys)),
c     im(wdens(nsys)).  
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
c     h2       In: real *8
c              Scaling parameter for enforcing growth conditions at
c              \infty
c              The condition is
c              \int w(t) ds (t) / h2 = Log Growth/h2
c 
c     swgt     In: real *8
c              weight of single layer
c
c------------------------------------------------------------------
c    OUTPUT
c    xmat      In: real *8
c              Updated matrix to enforce the growth conditions at
c              infinity and include the constant term in the 
c              representation. 
c              The rhs of the system is assumed to be unpacked
c              as -yvel(1),xvel(1).....-yvel(n),xvel(n) and the unknowns
c              are assumed to be re(wdens(1)), im(wdens(1)...
c              re(wdens(isys),im(wdens(nsys)).
c------------------------------------------------------------------

      implicit none
      integer k0,k,nsys,nd(0:k),ier
      integer i,j,l,istart,jstart,isource,itarg
      integer nbod,kbod

      real *8 xs(*),ys(*),dsdt(*),rnx(*),rny(*),rkappa(*)
      real *8 h(0:k),swgt,xmat(nsys,nsys),h2

      complex *16 zdis,zn,zw,zx
      real *8 pi

      pi = 4.0d0*datan(1.0d0)

      istart = 0
      isource = 0
c     Source loop
      do nbod = k0,k
         do i=1,nd(nbod)
            isource = isource + 1
c           Enforce log growth at infinity            
            zw = -0.25d0*h(nbod)*dsdt(isource)/(pi*h2)
            xmat(nsys-1,2*isource-1) = 
     1      xmat(nsys-1,2*isource-1) + swgt*dreal(zw)

            xmat(nsys,2*isource-1) = 
     1      xmat(nsys,2*isource-1) + swgt*dimag(zw)

            xmat(nsys-1,2*isource) = 
     1      xmat(nsys-1,2*isource) + swgt*dimag(-zw)

            xmat(nsys,2*isource) = 
     1      xmat(nsys,2*isource) + swgt*dreal(zw) 
c           Include constant term in representation
            xmat(2*isource-1,nsys - 1) = 
     1      xmat(2*isource-1,nsys - 1) + 1.0d0

            xmat(2*isource-1,nsys) = 
     1      xmat(2*isource-1,nsys) + 0.0d0

            xmat(2*isource,nsys - 1) = 
     1      xmat(2*isource,nsys - 1) + 0.0d0

            xmat(2*isource,nsys) = 
     1      xmat(2*isource,nsys) + 1.0d0
         enddo
         istart = istart + nd(nbod)
      enddo

      return
      end


      subroutine assemblepert(k0,k,nsys,nd,xs,ys,dsdt,rnx,rny,rkappa,h,
     1                      pwgt,xmat)
c     This subroutine updates the matrix entries for xmat
c     to incorportate the generalized 1's matrix to complete the
c     represenation for the interior problems.
c     The term added to the complex velocity is
c     zn(target)*\int w(t) dt/(2\pi) where zn(target) is the
c     complex normal nx + i ny at the target location.
c     The rhs of the system is assumed to be unpacked as
c     -yvel(1),xvel(1).....-yvel(n),xvel(n) and the unknowns are
c     assumed to be re(wdens(1)), im(wdens(1)....re(wdens(nsys)),
c     im(wdens(nsys)).  
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
c     pwgt     In: real *8
c              weight of generalized 1's matrix contribution
c
c------------------------------------------------------------------
c    OUTPUT
c    xmat      In: real *8
c              Updated matrix to include the generalized 1's
c              matrix to complete the representation for the
c              interior problem. The following term is added
c              to the complex velocity zn(target)*\int w(t) dt/(2\pi)
c              The rhs of the system is assumed to be unpacked
c              as -yvel(1),xvel(1).....-yvel(n),xvel(n) and the unknowns
c              are assumed to be re(wdens(1)), im(wdens(1)...
c              re(wdens(isys),im(wdens(nsys)).
c------------------------------------------------------------------

      implicit none
      integer k0,k,nsys,nd(0:k),ier
      integer i,j,l,istart,jstart,isource,itarg
      integer nbod,kbod

      real *8 xs(*),ys(*),dsdt(*),rnx(*),rny(*),rkappa(*)
      real *8 h(0:k),pwgt,xmat(nsys,nsys),wgt

      complex *16 zdis,zt,zw,zx,zn
      real *8 pi

      pi = 4.0d0*datan(1.0d0)
      itarg = 0
      jstart = 0
      do kbod = k0,k
         do j=1,nd(kbod)
            itarg = itarg + 1
            istart = 0
            isource = 0
c           Source loop
            do nbod = k0,k
               do i=1,nd(nbod)
                  isource = isource + 1
                  wgt = 0.5d0*h(nbod)*dsdt(isource)/pi
                  xmat(2*itarg-1,2*isource-1) = 
     1            xmat(2*itarg-1,2*isource-1)
     2            +pwgt*wgt*rny(isource)*rny(itarg)

                  xmat(2*itarg,2*isource-1) = 
     1            xmat(2*itarg,2*isource-1)
     2            - pwgt*wgt*rnx(itarg)*rny(isource)             

                  xmat(2*itarg-1,2*isource) = 
     1            xmat(2*itarg-1,2*isource)
     2            - pwgt*wgt*rny(itarg)*rnx(isource)             

                  xmat(2*itarg,2*isource) = 
     1            xmat(2*itarg,2*isource)
     2            + pwgt*wgt*rnx(itarg)*rnx(isource)              
               enddo
               istart = istart + nd(nbod)
            enddo
         enddo
         jstart = jstart + nd(kbod)
      enddo

      return
      end


