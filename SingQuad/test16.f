cc Copyright (C) 2010: Leslie Greengard and Mike O'Neil
cc Contact: greengard@cims.nyu.edu
cc Contact: oneil@cims.nyu.edu
c
c    This file contains matrix generators for the single layer potential,
c    the principal value of the double layer potential and the principal
c    value of the normal derivative of the single layer potential.
c
c    These are sometimes denoted S_k (SLP), D_k (DLP) and 
c    S_k' (normal deriv of SLP). 
c
c    formsprimematbac creates S_k'
c    formdlpmatbac    creates D_k
c    formslpmatbac    creates S_k
c
c    The bac suffix refers to the fact that we use Bradley Alpert
c    quadrature and that the matrices are complex.
c
c     The quadrature rule used is taken from:
c
c     B. Alpert,
c     Hybrid Gauss-Trapezoidal Quadrature Rules 
c     SIAM J. Sci. Comput. 20 (1999), pp. 1551-1584. 
c'
c
c
        subroutine formsprimematbac(ier,amat,norder,xs,ys,
     1      rnx,rny,dsdt,h,ns,gfun,eps,zk,mode)
        implicit real *8 (a-h,o-z)
        real *8 xs(1),ys(1),rnx(1),rny(1),dsdt(1),tpts(100),
     1      xpts(100),ypts(100),spts(100),
     2      txtra(100),coefs(100)
        complex *16 amat(ns,ns),zk,ima,u
        integer *4 its(100),its2(100)
        dimension extranodes4(6), extraweights4(6)
        dimension extranodes8(14), extraweights8(14)
        dimension extranodes16(30), extraweights16(30)
        dimension extranodes(30), extraweights(30)
c
c       this routine builds the matrix which applies the double
c       layer potential gfun to a vector using alpert quadrature
c
c     input:
c       ier - error code, anything other than 0 is bad
c       norder - the order of alpert quadrature to use, 0, 4, 8, 16
c       xs,ys - the x,y coordinates of points on the curve
c       rnx,rny - the unit normal derivative at the points xs,ys
c       dsdt - derivative of parameterization
c       h - step size in parameterization variable
c       ns - length of xs and ys
c       gfun - routine that evaluates the kernel, must have a
c           calling sequence of the form
c
c           gfun(par1,par2,par3,xtarg,ytarg,xsrc,ysrc,xp,yp,uval)
c
c           where par1,par2,par3 - arbitrary parameters
c                 xtarg,ytarg - target x,y values
c                 xsrc,ysrc - source x,y values
c                 xp,yp - unit normal direction to calculate derivative in
c                 uval - the complex *16 value of the kernel 
c                     at ((xtarg,ytarg),(xsrc,ysrc))
c
c       eps - the precision with which to evaluate the kernel, gets
c           passed as par1 above (usually)
c       zk - complex helmholtz parameter
c       mode - fourier mode to compute, used only in passing to gfun
c
c     output:
c       amat - the ns x ns matrix that will apply the normal derivative
c           of the single layer potential, i.e. S'
c'
c
      data extranodes4/
     1 -1.023715124251890D+00,
     1 -2.935370741501914D-01,
     1 -2.379647284118974D-02, 
     1 2.379647284118974D-02, 
     1 2.935370741501914D-01,
     1 1.023715124251890D+00/ 
      data extraweights4/
     1 9.131388579526912D-01,
     1 4.989017152913699D-01,
     1 8.795942675593887D-02,
     1 8.795942675593887D-02,
     1 4.989017152913699D-01,
     1 9.131388579526912D-01/
c
      data extranodes8/
     1 -3.998861349951123D+00,
     1 -2.980147933889640D+00,
     1 -1.945288592909266D+00,
     1 -1.027856640525646D+00,
     1 -3.967966533375878D-01,
     1 -9.086744584657729D-02,
     1 -6.531815708567918D-03,
     1 6.531815708567918D-03,
     1 9.086744584657729D-02,
     1 3.967966533375878D-01,
     1 1.027856640525646D+00,
     1 1.945288592909266D+00,
     1 2.980147933889640D+00,
     1 3.998861349951123D+00/
      data extraweights8/
     1 1.004787656533285D+00,
     1 1.036093649726216D+00,
     1 1.008710414337933D+00,
     1 7.947291148621894D-01,
     1 4.609256358650077D-01,
     1 1.701315866854178D-01,
     1 2.462194198995203D-02,
     1 2.462194198995203D-02,
     1 1.701315866854178D-01,
     1 4.609256358650077D-01,
     1 7.947291148621894D-01,
     1 1.008710414337933D+00,
     1 1.036093649726216D+00,
     1 1.004787656533285D+00/
c
      data extranodes16/
     1 -8.999998754306120d+00,
     1 -7.999888757524622d+00,
     1 -6.997957704791519d+00,
     1 -5.986360113977494d+00,
     1 -4.957203086563112d+00,
     1 -3.928129252248612d+00,
     1 -2.947904939031494d+00,
     1 -2.073471660264395d+00,
     1 -1.348993882467059d+00,
     1 -7.964747731112430d-01,
     1 -4.142832599028031d-01,
     1 -1.805991249601928d-01,
     1 -6.009290785739468d-02,
     1 -1.239382725542637d-02,
     1 -8.371529832014113d-04,
     1 8.371529832014113d-04,
     1 1.239382725542637d-02,
     1 6.009290785739468d-02,
     1 1.805991249601928d-01,
     1 4.142832599028031d-01,
     1 7.964747731112430d-01,
     1 1.348993882467059d+00,
     1 2.073471660264395d+00,
     1 2.947904939031494d+00,
     1 3.928129252248612d+00,
     1 4.957203086563112d+00,
     1 5.986360113977494d+00,
     1 6.997957704791519d+00,
     1 7.999888757524622d+00,
     1 8.999998754306120d+00/
      data extraweights16/
     1 1.000007149422537d+00,
     1 1.000395017352309d+00,
     1 1.004798397441514d+00,
     1 1.020308624984610d+00,
     1 1.035167721053657d+00,
     1 1.014359775369075d+00,
     1 9.362411945698647d-01,
     1 8.051212946181061d-01,
     1 6.401489637096768d-01,
     1 4.652220834914617d-01,
     1 3.029123478511309d-01,
     1 1.704889420286369d-01,
     1 7.740135521653088d-02,
     1 2.423621380426338d-02,
     1 3.190919086626234d-03,
     1 3.190919086626234d-03,
     1 2.423621380426338d-02,
     1 7.740135521653088d-02,
     1 1.704889420286369d-01,
     1 3.029123478511309d-01,
     1 4.652220834914617d-01,
     1 6.401489637096768d-01,
     1 8.051212946181061d-01,
     1 9.362411945698647d-01,
     1 1.014359775369075d+00,
     1 1.035167721053657d+00,
     1 1.020308624984610d+00,
     1 1.004798397441514d+00,
     1 1.000395017352309d+00,
     1 1.000007149422537d+00/
c
        ima=(0,1)
        done=1
        pi=4*atan(done)
c
        if (norder.eq. 0) then
        nskip=1
        nextra=0
        goto 1111
        endif
c
      if (norder.eq.4) then
         nskip = 2
         nextra = 6
	 do i = 1,nextra
	    extranodes(i) = extranodes4(i)
	    extraweights(i) = extraweights4(i)
         enddo
         goto 1111
      endif
c
      if (norder.eq.8) then 
         nskip = 5
         nextra = 14
	 do i = 1,nextra
	    extranodes(i) = extranodes8(i)
	    extraweights(i) = extraweights8(i)
         enddo
         goto 1111
      endif
c
      if (norder.eq.16) then 
         nskip = 10
         nextra = 30
	 do i = 1,nextra
	    extranodes(i) = extranodes16(i)
	    extraweights(i) = extraweights16(i)
         enddo
         goto 1111
      endif      

c
c       quadrature code not available for any other order -> return
c
        ier = 1
        call prinf('wrong order for quadrature, norder=*',norder,1)
        stop
        return
 1111   continue
        ier = 0

c
c       carry out "punctured" trapezoidal rule and fill in matrix
c       entries, skipping entries within nskip of the diagonal
c
        do 1400 j=1,ns
        do 1200 i=1,ns
        amat(i,j)=0
 1200   continue
 1400   continue

c
        n=ns-2*nskip+1
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,iii,k,u) 
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 

        do 1800 i=1,ns
        iii=i-1+nskip
c
        do 2000 k=0,n-1
c
        iii=iii+1
        if (iii .gt. ns) iii=iii-ns

        call gfun(eps,zk,mode,xs(i),ys(i),xs(iii),ys(iii),
     1      rnx(i),rny(i),u)
        amat(i,iii)=u*dsdt(iii)*h
 2000   continue
 1800   continue
c
C$OMP END PARALLEL DO

        if (norder .eq. 0) return

c
c       now add in corrections and interpolated stuff for alpert
c       first determine all the interpolation coefficients
c
        ninterp=norder+2
c
        do 6000 ipt=1,ns
c

        do 2400 i=1,nextra
        txtra(i)=h*(ipt-1)+h*extranodes(i)
 2400   continue

c
        do 4000 i=1,nextra
c
c       find the closest ninterp points to each of the txtra
c
        n1=txtra(i)/h
        if (txtra(i) .lt. 0) n1=n1-1
        n2=n1+1
        nnn=n1-(ninterp-2)/2

c
        do 2800 j=1,ninterp
        its(j)=nnn+j-1
        its2(j)=its(j)+1
        if (its2(j) .le. 0) its2(j)=its2(j)+ns
        if (its2(j) .gt. ns) its2(j)=its2(j)-ns
 2800   continue

c
c       fill interpolation nodes and function values
c
        do 3000 j=1,ninterp
        tpts(j)=its(j)*h
        xpts(j)=xs(its2(j))
        ypts(j)=ys(its2(j))
        spts(j)=dsdt(its2(j))
 3000   continue

c
c       now compute the values of xs, ys, dsdt at ttt using barycentric
c       interpolation
c
        ttt=txtra(i)
        call bary1_coefs(ninterp,tpts,ttt,coefs)

        xxx=0
        yyy=0
        sss=0
c
        do 3200 j=1,ninterp
        xxx=xxx+xpts(j)*coefs(j)
        yyy=yyy+ypts(j)*coefs(j)
        sss=sss+spts(j)*coefs(j)
 3200   continue

c
c       evaluate the kernel at the new quadrature point xxx,yyy and
c       add its contribution to the matrix at its interpolation points
c

        call gfun(eps,zk,mode,xs(ipt),ys(ipt),xxx,yyy,
     1      rnx(ipt),rny(ipt),u)
c
        do 3600 j=1,ninterp
        jjj=its2(j)
        amat(ipt,jjj)=amat(ipt,jjj)+u*sss*h*extraweights(i)*coefs(j)
 3600   continue
 4000   continue
c
 6000   continue
c
        return
        end
c
c
c
c
c
        subroutine formdlpmatbac(ier,amat,norder,xs,ys,rnx,rny,dsdt,
     1      h,ns,gfun,eps,zk,mode)
        implicit real *8 (a-h,o-z)
        real *8 xs(1),ys(1),rnx(1),rny(1),dsdt(1),tpts(100),
     1      xpts(100),ypts(100),spts(100),rnxpts(100),rnypts(100),
     2      txtra(100),coefs(100)
        complex *16 amat(ns,ns),zk,ima,u
        integer *4 its(100),its2(100)
        dimension extranodes4(6), extraweights4(6)
        dimension extranodes8(14), extraweights8(14)
        dimension extranodes16(30), extraweights16(30)
        dimension extranodes(30), extraweights(30)
c
c       this routine builds the matrix which applies the double
c       layer potential gfun to a vector using alpert quadrature
c
c     input:
c       ier - error code, anything other than 0 is bad
c       norder - the order of alpert quadrature to use, 0, 4, 8, 16
c       xs,ys - the x,y coordinates of points on the curve
c       rnx,rny - the unit normal derivatives at the points xs,ys
c       dsdt - derivative of parameterization
c       h - step size in parameterization variable
c       ns - length of xs and ys
c       gfun - routine that evaluates the kernel, must have a
c           calling sequence of the form
c
c           gfun(par1,par2,par3,xtarg,ytarg,xsrc,ysrc,xp,yp,uval)
c
c           where par1,par2,par3 - arbitrary parameters
c                 xtarg,ytarg - target x,y values
c                 xsrc,ysrc - source x,y values
c                 xp,yp - unit normal direction to calculate derivative in
c                 uval - the complex *16 value of the kernel 
c                     at ((xtarg,ytarg),(xsrc,ysrc))
c
c       eps - the precision with which to evaluate the kernel, gets
c           passed as par1 above (usually)
c       zk - complex helmholtz parameter
c       mode - fourier mode to compute, used only in passing to gfun
c
c     output:
c       amat - the ns x ns matrix that will apply the double layer
c           potential
c
c
      data extranodes4/
     1 -1.023715124251890D+00,
     1 -2.935370741501914D-01,
     1 -2.379647284118974D-02, 
     1 2.379647284118974D-02, 
     1 2.935370741501914D-01,
     1 1.023715124251890D+00/ 
      data extraweights4/
     1 9.131388579526912D-01,
     1 4.989017152913699D-01,
     1 8.795942675593887D-02,
     1 8.795942675593887D-02,
     1 4.989017152913699D-01,
     1 9.131388579526912D-01/
c
      data extranodes8/
     1 -3.998861349951123D+00,
     1 -2.980147933889640D+00,
     1 -1.945288592909266D+00,
     1 -1.027856640525646D+00,
     1 -3.967966533375878D-01,
     1 -9.086744584657729D-02,
     1 -6.531815708567918D-03,
     1 6.531815708567918D-03,
     1 9.086744584657729D-02,
     1 3.967966533375878D-01,
     1 1.027856640525646D+00,
     1 1.945288592909266D+00,
     1 2.980147933889640D+00,
     1 3.998861349951123D+00/
      data extraweights8/
     1 1.004787656533285D+00,
     1 1.036093649726216D+00,
     1 1.008710414337933D+00,
     1 7.947291148621894D-01,
     1 4.609256358650077D-01,
     1 1.701315866854178D-01,
     1 2.462194198995203D-02,
     1 2.462194198995203D-02,
     1 1.701315866854178D-01,
     1 4.609256358650077D-01,
     1 7.947291148621894D-01,
     1 1.008710414337933D+00,
     1 1.036093649726216D+00,
     1 1.004787656533285D+00/
c
      data extranodes16/
     1 -8.999998754306120d+00,
     1 -7.999888757524622d+00,
     1 -6.997957704791519d+00,
     1 -5.986360113977494d+00,
     1 -4.957203086563112d+00,
     1 -3.928129252248612d+00,
     1 -2.947904939031494d+00,
     1 -2.073471660264395d+00,
     1 -1.348993882467059d+00,
     1 -7.964747731112430d-01,
     1 -4.142832599028031d-01,
     1 -1.805991249601928d-01,
     1 -6.009290785739468d-02,
     1 -1.239382725542637d-02,
     1 -8.371529832014113d-04,
     1 8.371529832014113d-04,
     1 1.239382725542637d-02,
     1 6.009290785739468d-02,
     1 1.805991249601928d-01,
     1 4.142832599028031d-01,
     1 7.964747731112430d-01,
     1 1.348993882467059d+00,
     1 2.073471660264395d+00,
     1 2.947904939031494d+00,
     1 3.928129252248612d+00,
     1 4.957203086563112d+00,
     1 5.986360113977494d+00,
     1 6.997957704791519d+00,
     1 7.999888757524622d+00,
     1 8.999998754306120d+00/
      data extraweights16/
     1 1.000007149422537d+00,
     1 1.000395017352309d+00,
     1 1.004798397441514d+00,
     1 1.020308624984610d+00,
     1 1.035167721053657d+00,
     1 1.014359775369075d+00,
     1 9.362411945698647d-01,
     1 8.051212946181061d-01,
     1 6.401489637096768d-01,
     1 4.652220834914617d-01,
     1 3.029123478511309d-01,
     1 1.704889420286369d-01,
     1 7.740135521653088d-02,
     1 2.423621380426338d-02,
     1 3.190919086626234d-03,
     1 3.190919086626234d-03,
     1 2.423621380426338d-02,
     1 7.740135521653088d-02,
     1 1.704889420286369d-01,
     1 3.029123478511309d-01,
     1 4.652220834914617d-01,
     1 6.401489637096768d-01,
     1 8.051212946181061d-01,
     1 9.362411945698647d-01,
     1 1.014359775369075d+00,
     1 1.035167721053657d+00,
     1 1.020308624984610d+00,
     1 1.004798397441514d+00,
     1 1.000395017352309d+00,
     1 1.000007149422537d+00/
c
        ima=(0,1)
        done=1
        pi=4*atan(done)
c
        if (norder.eq. 0) then
        nskip=1
        nextra=0
        goto 1111
        endif
c
      if (norder.eq.4) then
         nskip = 2
         nextra = 6
	 do i = 1,nextra
	    extranodes(i) = extranodes4(i)
	    extraweights(i) = extraweights4(i)
         enddo
         goto 1111
      endif
c
      if (norder.eq.8) then 
         nskip = 5
         nextra = 14
	 do i = 1,nextra
	    extranodes(i) = extranodes8(i)
	    extraweights(i) = extraweights8(i)
         enddo
         goto 1111
      endif
c
      if (norder.eq.16) then 
         nskip = 10
         nextra = 30
	 do i = 1,nextra
	    extranodes(i) = extranodes16(i)
	    extraweights(i) = extraweights16(i)
         enddo
         goto 1111
      endif      

c
c       quadrature code not available for any other order -> return
c
        ier = 1
        call prinf('wrong order for quadrature, norder=*',norder,1)
        stop
        return
 1111   continue
        ier = 0

c
c       carry out "punctured" trapezoidal rule and fill in matrix
c       entries, skipping entries within nskip of the diagonal
c
        do 1400 j=1,ns
        do 1200 i=1,ns
        amat(i,j)=0
 1200   continue
 1400   continue

c
        n=ns-2*nskip+1
c

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,iii,k,u) 
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 

        do 1800 i=1,ns
        iii=i-1+nskip
c
        do 2000 k=0,n-1
c
        iii=iii+1
        if (iii .gt. ns) iii=iii-ns

        call gfun(eps,zk,mode,xs(i),ys(i),xs(iii),ys(iii),
     1      rnx(iii),rny(iii),u)
        amat(i,iii)=u*dsdt(iii)*h
 2000   continue
 1800   continue
c
C$OMP END PARALLEL DO


        if (norder .eq. 0) return

c
c       now add in corrections and interpolated stuff for alpert
c       first determine all the interpolation coefficients
c
        ninterp=norder+2
c
        do 6000 ipt=1,ns
c

        do 2400 i=1,nextra
        txtra(i)=h*(ipt-1)+h*extranodes(i)
 2400   continue

c
        do 4000 i=1,nextra
c
c       find the closest ninterp points to each of the txtra
c
        n1=txtra(i)/h
        if (txtra(i) .lt. 0) n1=n1-1
        n2=n1+1
        nnn=n1-(ninterp-2)/2

c
        do 2800 j=1,ninterp
        its(j)=nnn+j-1
        its2(j)=its(j)+1
        if (its2(j) .le. 0) its2(j)=its2(j)+ns
        if (its2(j) .gt. ns) its2(j)=its2(j)-ns
 2800   continue

c
c       fill interpolation nodes and function values
c
        do 3000 j=1,ninterp
        tpts(j)=its(j)*h
        xpts(j)=xs(its2(j))
        ypts(j)=ys(its2(j))
        spts(j)=dsdt(its2(j))
        rnxpts(j)=rnx(its2(j))
        rnypts(j)=rny(its2(j))
 3000   continue

c
c       now compute the values of xs, ys, dsdt at ttt using barycentric
c       interpolation
c
        ttt=txtra(i)
        call bary1_coefs(ninterp,tpts,ttt,coefs)

        xxx=0
        yyy=0
        sss=0
        rnxxx=0
        rnyyy=0
c
        do 3200 j=1,ninterp
        xxx=xxx+xpts(j)*coefs(j)
        yyy=yyy+ypts(j)*coefs(j)
        sss=sss+spts(j)*coefs(j)
        rnxxx=rnxxx+rnxpts(j)*coefs(j)
        rnyyy=rnyyy+rnypts(j)*coefs(j)
 3200   continue

c
c       evaluate the kernel at the new quadrature point xxx,yyy and
c       add its contribution to the matrix at its interpolation points
c

        call gfun(eps,zk,mode,xs(ipt),ys(ipt),xxx,yyy,rnxxx,rnyyy,u)
c
        do 3600 j=1,ninterp
        jjj=its2(j)
        amat(ipt,jjj)=amat(ipt,jjj)+u*sss*h*extraweights(i)*coefs(j)
 3600   continue
 4000   continue
c
 6000   continue
c
        return
        end
c
c
c
c
c
        subroutine formslpmatbac(ier,amat,norder,xs,ys,dsdt,
     1      h,ns,gfun,eps,zk,mode)
        implicit real *8 (a-h,o-z)
        real *8 xs(1),ys(1),dsdt(1),tpts(100),xpts(100),
     1      ypts(100),spts(100),txtra(100),coefs(100)
        complex *16 ima,amat(ns,ns),zk,u
        integer *4 its(100),its2(100)
        dimension extranodes4(6), extraweights4(6)
        dimension extranodes8(14), extraweights8(14)
        dimension extranodes16(30), extraweights16(30)
        dimension extranodes(30), extraweights(30)
c
c       this routine builds the matrix which applies the single
c       layer potential gfun to a vector using alpert quadrature
c
c     input:
c       ier - error code, anything other than 0 is bad
c       norder - the order of alpert quadrature to use, 0, 4, 8, 16
c       xs,ys - the x,y coordinates of points on the curve
c       dsdt - derivative of parameterization
c       h - step size in parameterization variable
c       ns - length of xs and ys
c       gfun - routine that evaluates the kernel, must have a
c           calling sequence of the form
c
c           gfun(par1,par2,par3,xtarg,ytarg,xsrc,ysrc,uval)
c
c           where par1,par2,par3 - arbitrary parameters
c                 xtarg,ytarg - target x,y values
c                 xsrc,ysrc - source x,y values
c                 uval - the complex *16 value of the kernel 
c                     at ((xtarg,ytarg),(xsrc,ysrc))
c
c       eps - the precision with which to evaluate the kernel, gets
c           passed as par1 above (usually)
c       zk - complex helmholtz parameter
c       mode - fourier mode to compute, used only in passing to gfun
c
c     output:
c       amat - the ns x ns matrix that will apply the single layer
c           potential
c
c
      data extranodes4/
     1 -1.023715124251890D+00,
     1 -2.935370741501914D-01,
     1 -2.379647284118974D-02, 
     1 2.379647284118974D-02, 
     1 2.935370741501914D-01,
     1 1.023715124251890D+00/ 
      data extraweights4/
     1 9.131388579526912D-01,
     1 4.989017152913699D-01,
     1 8.795942675593887D-02,
     1 8.795942675593887D-02,
     1 4.989017152913699D-01,
     1 9.131388579526912D-01/
c
      data extranodes8/
     1 -3.998861349951123D+00,
     1 -2.980147933889640D+00,
     1 -1.945288592909266D+00,
     1 -1.027856640525646D+00,
     1 -3.967966533375878D-01,
     1 -9.086744584657729D-02,
     1 -6.531815708567918D-03,
     1 6.531815708567918D-03,
     1 9.086744584657729D-02,
     1 3.967966533375878D-01,
     1 1.027856640525646D+00,
     1 1.945288592909266D+00,
     1 2.980147933889640D+00,
     1 3.998861349951123D+00/
      data extraweights8/
     1 1.004787656533285D+00,
     1 1.036093649726216D+00,
     1 1.008710414337933D+00,
     1 7.947291148621894D-01,
     1 4.609256358650077D-01,
     1 1.701315866854178D-01,
     1 2.462194198995203D-02,
     1 2.462194198995203D-02,
     1 1.701315866854178D-01,
     1 4.609256358650077D-01,
     1 7.947291148621894D-01,
     1 1.008710414337933D+00,
     1 1.036093649726216D+00,
     1 1.004787656533285D+00/
c
      data extranodes16/
     1 -8.999998754306120d+00,
     1 -7.999888757524622d+00,
     1 -6.997957704791519d+00,
     1 -5.986360113977494d+00,
     1 -4.957203086563112d+00,
     1 -3.928129252248612d+00,
     1 -2.947904939031494d+00,
     1 -2.073471660264395d+00,
     1 -1.348993882467059d+00,
     1 -7.964747731112430d-01,
     1 -4.142832599028031d-01,
     1 -1.805991249601928d-01,
     1 -6.009290785739468d-02,
     1 -1.239382725542637d-02,
     1 -8.371529832014113d-04,
     1 8.371529832014113d-04,
     1 1.239382725542637d-02,
     1 6.009290785739468d-02,
     1 1.805991249601928d-01,
     1 4.142832599028031d-01,
     1 7.964747731112430d-01,
     1 1.348993882467059d+00,
     1 2.073471660264395d+00,
     1 2.947904939031494d+00,
     1 3.928129252248612d+00,
     1 4.957203086563112d+00,
     1 5.986360113977494d+00,
     1 6.997957704791519d+00,
     1 7.999888757524622d+00,
     1 8.999998754306120d+00/
      data extraweights16/
     1 1.000007149422537d+00,
     1 1.000395017352309d+00,
     1 1.004798397441514d+00,
     1 1.020308624984610d+00,
     1 1.035167721053657d+00,
     1 1.014359775369075d+00,
     1 9.362411945698647d-01,
     1 8.051212946181061d-01,
     1 6.401489637096768d-01,
     1 4.652220834914617d-01,
     1 3.029123478511309d-01,
     1 1.704889420286369d-01,
     1 7.740135521653088d-02,
     1 2.423621380426338d-02,
     1 3.190919086626234d-03,
     1 3.190919086626234d-03,
     1 2.423621380426338d-02,
     1 7.740135521653088d-02,
     1 1.704889420286369d-01,
     1 3.029123478511309d-01,
     1 4.652220834914617d-01,
     1 6.401489637096768d-01,
     1 8.051212946181061d-01,
     1 9.362411945698647d-01,
     1 1.014359775369075d+00,
     1 1.035167721053657d+00,
     1 1.020308624984610d+00,
     1 1.004798397441514d+00,
     1 1.000395017352309d+00,
     1 1.000007149422537d+00/
c
c        print *, "In formslp"
        ima=(0,1)
        done=1
        pi=4*atan(done)
c
        if (norder.eq. 0) then
        nskip=1
        nextra=0
        goto 1111
        endif
c
      if (norder.eq.4) then
         nskip = 2
         nextra = 6
	 do i = 1,nextra
	    extranodes(i) = extranodes4(i)
	    extraweights(i) = extraweights4(i)
         enddo
         goto 1111
      endif
c
      if (norder.eq.8) then 
         nskip = 5
         nextra = 14
	 do i = 1,nextra
	    extranodes(i) = extranodes8(i)
	    extraweights(i) = extraweights8(i)
         enddo
         goto 1111
      endif
c
      if (norder.eq.16) then 
         nskip = 10
         nextra = 30
	 do i = 1,nextra
	    extranodes(i) = extranodes16(i)
	    extraweights(i) = extraweights16(i)
         enddo
         goto 1111
      endif      
c
c       quadrature code not available for any other order -> return
c
        ier = 1
        call prinf('wrong order for quadrature, norder=*',norder,1)
        stop
        return
 1111   continue
        ier = 0
        do 1400 j=1,ns
        do 1200 i=1,ns
        amat(i,j)=0
 1200   continue
 1400   continue
c
        n=ns-2*nskip+1
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,iii,k,u) 
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do 1800 i=1,ns
        iii=i-1+nskip
c
        do 2000 k=0,n-1
c
        iii=iii+1
        if (iii .gt. ns) iii=iii-ns

        call gfun(eps,zk,mode,xs(i),ys(i),xs(iii),ys(iii),u)
        amat(i,iii)=u*dsdt(iii)*h
 2000   continue
 1800   continue
c
C$OMP END PARALLEL DO


        if (norder .eq. 0) return

c
c       now add in corrections and interpolated stuff for alpert
c       first determine all the interpolation coefficients
c
        ninterp=norder+2
c
        do 6000 ipt=1,ns
c

        do 2400 i=1,nextra
        txtra(i)=h*(ipt-1)+h*extranodes(i)
 2400   continue

c
        do 4000 i=1,nextra
c
c       find the closest ninterp points to each of the txtra
c
        n1=txtra(i)/h
        if (txtra(i) .lt. 0) n1=n1-1
        n2=n1+1
        nnn=n1-(ninterp-2)/2

c
        do 2800 j=1,ninterp
        its(j)=nnn+j-1
        its2(j)=its(j)+1
        if (its2(j) .le. 0) its2(j)=its2(j)+ns
        if (its2(j) .gt. ns) its2(j)=its2(j)-ns
 2800   continue

c
c       fill interpolation nodes and function values
c
        do 3000 j=1,ninterp
        tpts(j)=its(j)*h
        xpts(j)=xs(its2(j))
        ypts(j)=ys(its2(j))
        spts(j)=dsdt(its2(j))
 3000   continue

c
c       now compute the values of xs, ys, dsdt at ttt using barycentric
c       interpolation
c
        ttt=txtra(i)
        call bary1_coefs(ninterp,tpts,ttt,coefs)

        xxx=0
        yyy=0
        sss=0
c
        do 3200 j=1,ninterp
        xxx=xxx+xpts(j)*coefs(j)
        yyy=yyy+ypts(j)*coefs(j)
        sss=sss+spts(j)*coefs(j)
 3200   continue

c
c       evaluate the kernel at the new quadrature point xxx,yyy and
c       add its contribution to the matrix at its interpolation points
c

        call gfun(eps,zk,mode,xs(ipt),ys(ipt),xxx,yyy,u)
c
        do 3600 j=1,ninterp
        jjj=its2(j)
        amat(ipt,jjj)=amat(ipt,jjj)+u*sss*h*extraweights(i)*coefs(j)
 3600   continue
 4000   continue
c
 6000   continue
c
        return
        end
c
c
        subroutine formslpmatcorr(ier,amat,norder,xs,ys,dsdt,
     1      h,ns,gfun,eps,zk,mode)
        implicit real *8 (a-h,o-z)
        real *8 xs(1),ys(1),dsdt(1),tpts(100),xpts(100),
     1      ypts(100),spts(100),txtra(100),coefs(100)
        complex *16 ima,amat(ns,ns),zk,u,u1
        integer its(100),its2(100),nskip2,norder
        dimension extranodes4(6), extraweights4(6)
        dimension extranodes8(14), extraweights8(14)
        dimension extranodes16(30), extraweights16(30)
        dimension extranodes(30), extraweights(30)
c
c       this routine builds the matrix which applies the single
c       layer potential gfun to a vector using alpert quadrature
c
c     input:
c       ier - error code, anything other than 0 is bad
c       norder - the order of alpert quadrature to use, 0, 4, 8, 16
c       xs,ys - the x,y coordinates of points on the curve
c       dsdt - derivative of parameterization
c       h - step size in parameterization variable
c       ns - length of xs and ys
c       gfun - routine that evaluates the kernel, must have a
c           calling sequence of the form
c
c           gfun(par1,par2,par3,xtarg,ytarg,xsrc,ysrc,uval)
c
c           where par1,par2,par3 - arbitrary parameters
c                 xtarg,ytarg - target x,y values
c                 xsrc,ysrc - source x,y values
c                 uval - the complex *16 value of the kernel 
c                     at ((xtarg,ytarg),(xsrc,ysrc))
c
c       eps - the precision with which to evaluate the kernel, gets
c           passed as par1 above (usually)
c       zk - complex helmholtz parameter
c       mode - fourier mode to compute, used only in passing to gfun
c
c     output:
c       amat - the ns x ns matrix that will apply the single layer
c           potential
c
c
      data extranodes4/
     1 -1.023715124251890D+00,
     1 -2.935370741501914D-01,
     1 -2.379647284118974D-02, 
     1 2.379647284118974D-02, 
     1 2.935370741501914D-01,
     1 1.023715124251890D+00/ 
      data extraweights4/
     1 9.131388579526912D-01,
     1 4.989017152913699D-01,
     1 8.795942675593887D-02,
     1 8.795942675593887D-02,
     1 4.989017152913699D-01,
     1 9.131388579526912D-01/
c
      data extranodes8/
     1 -3.998861349951123D+00,
     1 -2.980147933889640D+00,
     1 -1.945288592909266D+00,
     1 -1.027856640525646D+00,
     1 -3.967966533375878D-01,
     1 -9.086744584657729D-02,
     1 -6.531815708567918D-03,
     1 6.531815708567918D-03,
     1 9.086744584657729D-02,
     1 3.967966533375878D-01,
     1 1.027856640525646D+00,
     1 1.945288592909266D+00,
     1 2.980147933889640D+00,
     1 3.998861349951123D+00/
      data extraweights8/
     1 1.004787656533285D+00,
     1 1.036093649726216D+00,
     1 1.008710414337933D+00,
     1 7.947291148621894D-01,
     1 4.609256358650077D-01,
     1 1.701315866854178D-01,
     1 2.462194198995203D-02,
     1 2.462194198995203D-02,
     1 1.701315866854178D-01,
     1 4.609256358650077D-01,
     1 7.947291148621894D-01,
     1 1.008710414337933D+00,
     1 1.036093649726216D+00,
     1 1.004787656533285D+00/
c
      data extranodes16/
     1 -8.999998754306120d+00,
     1 -7.999888757524622d+00,
     1 -6.997957704791519d+00,
     1 -5.986360113977494d+00,
     1 -4.957203086563112d+00,
     1 -3.928129252248612d+00,
     1 -2.947904939031494d+00,
     1 -2.073471660264395d+00,
     1 -1.348993882467059d+00,
     1 -7.964747731112430d-01,
     1 -4.142832599028031d-01,
     1 -1.805991249601928d-01,
     1 -6.009290785739468d-02,
     1 -1.239382725542637d-02,
     1 -8.371529832014113d-04,
     1 8.371529832014113d-04,
     1 1.239382725542637d-02,
     1 6.009290785739468d-02,
     1 1.805991249601928d-01,
     1 4.142832599028031d-01,
     1 7.964747731112430d-01,
     1 1.348993882467059d+00,
     1 2.073471660264395d+00,
     1 2.947904939031494d+00,
     1 3.928129252248612d+00,
     1 4.957203086563112d+00,
     1 5.986360113977494d+00,
     1 6.997957704791519d+00,
     1 7.999888757524622d+00,
     1 8.999998754306120d+00/
      data extraweights16/
     1 1.000007149422537d+00,
     1 1.000395017352309d+00,
     1 1.004798397441514d+00,
     1 1.020308624984610d+00,
     1 1.035167721053657d+00,
     1 1.014359775369075d+00,
     1 9.362411945698647d-01,
     1 8.051212946181061d-01,
     1 6.401489637096768d-01,
     1 4.652220834914617d-01,
     1 3.029123478511309d-01,
     1 1.704889420286369d-01,
     1 7.740135521653088d-02,
     1 2.423621380426338d-02,
     1 3.190919086626234d-03,
     1 3.190919086626234d-03,
     1 2.423621380426338d-02,
     1 7.740135521653088d-02,
     1 1.704889420286369d-01,
     1 3.029123478511309d-01,
     1 4.652220834914617d-01,
     1 6.401489637096768d-01,
     1 8.051212946181061d-01,
     1 9.362411945698647d-01,
     1 1.014359775369075d+00,
     1 1.035167721053657d+00,
     1 1.020308624984610d+00,
     1 1.004798397441514d+00,
     1 1.000395017352309d+00,
     1 1.000007149422537d+00/
c
        ima=(0,1)
        done=1
        pi=4*atan(done)
c
        if (norder.eq. 0) then
        nskip=1
        nextra=0
        goto 1111
        endif
c
      if (norder.eq.4) then
         nskip = 2
         nextra = 6
	 do i = 1,nextra
	    extranodes(i) = extranodes4(i)
	    extraweights(i) = extraweights4(i)
         enddo
         goto 1111
      endif
c
      if (norder.eq.8) then 
         nskip = 5
         nextra = 14
	 do i = 1,nextra
	    extranodes(i) = extranodes8(i)
	    extraweights(i) = extraweights8(i)
         enddo
         goto 1111
      endif
c
      if (norder.eq.16) then
         nskip = 10
         nextra = 30
	 do i = 1,nextra
	    extranodes(i) = extranodes16(i)
	    extraweights(i) = extraweights16(i)
         enddo
         goto 1111
      endif      
c
c       quadrature code not available for any other order -> return
c
        ier = 1
        call prinf('wrong order for quadrature, norder=*',norder,1)
        stop
        return
 1111   continue
        ier = 0
        do 1400 j=1,ns
        do 1200 i=1,ns
        amat(i,j)=0
 1200   continue
 1400   continue
c
        n=ns-2*nskip+1
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,iii,k,u) 
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
cc        do 1800 i=1,ns
cc        iii=i-1+nskip
c
cc        do 2000 k=0,n-1
c
cc        iii=iii+1
cc        if (iii .gt. ns) iii=iii-ns

cc        call gfun(eps,zk,mode,xs(i),ys(i),xs(iii),ys(iii),u)
cc        amat(i,iii)=u*dsdt(iii)*h
cc 2000   continue
cc 1800   continue
c
C$OMP END PARALLEL DO


        if (norder .eq. 0) return

c
c       now add in corrections and interpolated stuff for alpert
c       first determine all the interpolation coefficients
c
        ninterp=norder+2
c
        do 6000 ipt=1,ns
c

        do 2400 i=1,nextra
        txtra(i)=h*(ipt-1)+h*extranodes(i)
 2400   continue

c
        do 4000 i=1,nextra
c
c       find the closest ninterp points to each of the txtra
c
        n1=txtra(i)/h
        if (txtra(i) .lt. 0) n1=n1-1
        n2=n1+1
        nnn=n1-(ninterp-2)/2

c
        do 2800 j=1,ninterp
        its(j)=nnn+j-1
        its2(j)=its(j)+1
        if (its2(j) .le. 0) its2(j)=its2(j)+ns
        if (its2(j) .gt. ns) its2(j)=its2(j)-ns
 2800   continue

c
c       fill interpolation nodes and function values
c
        do 3000 j=1,ninterp
        tpts(j)=its(j)*h
        xpts(j)=xs(its2(j))
        ypts(j)=ys(its2(j))
        spts(j)=dsdt(its2(j))
 3000   continue

c
c       now compute the values of xs, ys, dsdt at ttt using barycentric
c       interpolation
c
        ttt=txtra(i)
        call bary1_coefs(ninterp,tpts,ttt,coefs)

        xxx=0
        yyy=0
        sss=0
c
        do 3200 j=1,ninterp
        xxx=xxx+xpts(j)*coefs(j)
        yyy=yyy+ypts(j)*coefs(j)
        sss=sss+spts(j)*coefs(j)
 3200   continue

c
c       evaluate the kernel at the new quadrature point xxx,yyy and
c       add its contribution to the matrix at its interpolation points
c

        call gfun(eps,zk,mode,xs(ipt),ys(ipt),xxx,yyy,u)
c
        do 3600 j=1,ninterp
           jjj=its2(j)
cc           if(ipt.ne.jjj) then
cc              call gfun(eps,zk,mode,xs(ipt),ys(ipt),xs(jjj),ys(jjj),u1)
cc              amat(ipt,jjj)=amat(ipt,jjj)-u1*h*dsdt(jjj)
cc           endif
           amat(ipt,jjj)=amat(ipt,jjj)+u*sss*h*extraweights(i)*coefs(j)
     
 3600   continue
 4000   continue
        do j=1,nskip-1
            jjj = ipt+j
            if(jjj.le.0) jjj = jjj+ns
            if(jjj.gt.ns) jjj=jjj-ns
            call gfun(eps,zk,mode,xs(ipt),ys(ipt),xs(jjj),ys(jjj),u1)
            amat(ipt,jjj)=amat(ipt,jjj)-u1*h*dsdt(jjj)
        enddo

        do j=1,nskip-1
            jjj = ipt-j
            if(jjj.le.0) jjj = jjj+ns
            if(jjj.gt.ns) jjj=jjj-ns
            call gfun(eps,zk,mode,xs(ipt),ys(ipt),xs(jjj),ys(jjj),u1)
            amat(ipt,jjj)=amat(ipt,jjj)-u1*h*dsdt(jjj)
        enddo
c
 6000   continue
c
        return
        end
c
c
c
c
c
        subroutine bary1_coefs(n,ts,ttt,coefs)
        implicit real *8 (a-h,o-z)
        real *8 ts(1),ttt,whts(1000),coefs(1)
c
c       use barycentric interpolation on ts,xs to evaluate at ttt
c       first evaluate the barycentric weights (real routine)
c
c       input:
c         n - the length of ts,xs
c         ts - nodes with which to interpolate
c         ttt - point at which to interpolate
c
c       output:
c         coefs - coefficients for the interpolation
c
c
        do 1200 i=1,n
        whts(i)=1
 1200   continue
c
        do 1600 i=1,n
        do 1400 j=1,n
        if (i .ne. j) whts(i)=whts(i)/(ts(i)-ts(j))
 1400   continue
 1600   continue
c
c       this uses the '2nd form' of barycentric interpolation, first
c       form the denominator
c
        dd=0
        do 2000 i=1,n
        dd=dd+whts(i)/(ttt-ts(i))
 2000   continue
c
c       and next the interpolation coefficients
c
        do 2400 i=1,n
        coefs(i)=whts(i)/(ttt-ts(i))/dd
 2400   continue
c
        return
        end
c
c
c
c
c
      subroutine alpertslpdirect(ier,norder,xs,ys,dsdt,
     1    h,ns,gfun,eps,zk,mode,fs,vals)
      implicit real *8 (a-h,o-z)
      real *8 xs(1),ys(1),dsdt(1)
      complex *16 fs(1),vals(1),zk
      external gfun
c
c     this routine applies the single layer kernel to the
c     vector fs directly, i.e. the entire matrix is never formed
c
c     input:
c       norder - the order of alpert quadrature to use, 0, 4, 8, 16
c       xs,ys - the x,y coordinates of points on the curve
c       dsdt - derivative of parameterization
c       h - step size in parameterization variable
c       ns - lenght of xs and ys
c       gfun - routine that evaluates the kernel, must have a
c           calling sequence of the form
c
c           gfun(par1,par2,par3,xtarg,ytarg,xsrc,ysrc,uval)
c
c           where par1,par2,par3 - arbitrary parameters
c                 xtarg,ytarg - target x,y values
c                 xsrc,ysrc - source x,y values
c                 uval - the complex *16 value of the kernel 
c                     at ((xtarg,ytarg),(xsrc,ysrc))
c
c       eps - the precision with which to evaluate the kernel, gets
c           passed as par1 above (usually)
c       zk - complex helmholtz parameter
c       fs - complex vector to apply the single layer kernel gfun to
c
c     output:
c       vals - the value obtained by applying the single layer to fs
c
c
c     just call the 1 point evaluator n times
c

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i) 
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 

      do 2000 i=1,ns
      call alpertslp1(ier,norder,xs,ys,dsdt,
     1    h,ns,gfun,eps,zk,mode,fs,i,vals(i))
 2000 continue

c$omp end parallel do
c
      return
      end
c
c
c
c
c
      subroutine alpertslp1(ier,norder,xs,ys,dsdt,
     1    h,ns,gfun,eps,zk,mode,fs,ipt,val)
      implicit real *8 (a-h,o-z)
      dimension xs(ns),ys(ns)
      dimension dsdt(ns)
      complex *16 val,fs(1),fpts(10000),fff,ima
      complex *16 gg, zk, u
      integer *4 its(10000),its2(10000)
      dimension extranodes4(6), extraweights4(6)
      dimension extranodes8(14), extraweights8(14)
      dimension extranodes16(30), extraweights16(30)
      dimension extranodes(30), extraweights(30)
      real *8 tpts(10000),xpts(10000),ypts(10000),spts(10000),
     1    txtra(10000)
c
c     this routine computes 1 single layer entry directly, at
c     the target point ipt
c
c     input:
c       norder - the order of alpert quadrature to use, 0, 4, 8, 16
c       xs,ys - the x,y coordinates of points on the curve
c       dsdt - derivative of parameterization
c       h - step size in parameterization variable
c       ns - lenght of xs and ys
c       gfun - routine that evaluates the kernel, must have a
c           calling sequence of the form
c
c           gfun(par1,par2,par3,xtarg,ytarg,xsrc,ysrc,uval)
c
c           where par1,par2,par3 - arbitrary parameters
c                 xtarg,ytarg - target x,y values
c                 xsrc,ysrc - source x,y values
c                 uval - the complex *16 value of the kernel 
c                     at ((xtarg,ytarg),(xsrc,ysrc))
c
c       eps - the precision with which to evaluate the kernel, gets
c           passed as par1 above (usually)
c       zk - complex helmholtz parameter
c       fs - complex vector to apply the single layer kernel gfun to
c       ipt - integer determine what index in xs,ys is the target point
c
c     output:
c       vals - the value obtained by applying the single layer to fs
c
c
      data extranodes4/
     1 -1.023715124251890D+00,
     1 -2.935370741501914D-01,
     1 -2.379647284118974D-02, 
     1 2.379647284118974D-02, 
     1 2.935370741501914D-01,
     1 1.023715124251890D+00/ 
      data extraweights4/
     1 9.131388579526912D-01,
     1 4.989017152913699D-01,
     1 8.795942675593887D-02,
     1 8.795942675593887D-02,
     1 4.989017152913699D-01,
     1 9.131388579526912D-01/
c
      data extranodes8/
     1 -3.998861349951123D+00,
     1 -2.980147933889640D+00,
     1 -1.945288592909266D+00,
     1 -1.027856640525646D+00,
     1 -3.967966533375878D-01,
     1 -9.086744584657729D-02,
     1 -6.531815708567918D-03,
     1 6.531815708567918D-03,
     1 9.086744584657729D-02,
     1 3.967966533375878D-01,
     1 1.027856640525646D+00,
     1 1.945288592909266D+00,
     1 2.980147933889640D+00,
     1 3.998861349951123D+00/
      data extraweights8/
     1 1.004787656533285D+00,
     1 1.036093649726216D+00,
     1 1.008710414337933D+00,
     1 7.947291148621894D-01,
     1 4.609256358650077D-01,
     1 1.701315866854178D-01,
     1 2.462194198995203D-02,
     1 2.462194198995203D-02,
     1 1.701315866854178D-01,
     1 4.609256358650077D-01,
     1 7.947291148621894D-01,
     1 1.008710414337933D+00,
     1 1.036093649726216D+00,
     1 1.004787656533285D+00/
c
      data extranodes16/
     1 -8.999998754306120d+00,
     1 -7.999888757524622d+00,
     1 -6.997957704791519d+00,
     1 -5.986360113977494d+00,
     1 -4.957203086563112d+00,
     1 -3.928129252248612d+00,
     1 -2.947904939031494d+00,
     1 -2.073471660264395d+00,
     1 -1.348993882467059d+00,
     1 -7.964747731112430d-01,
     1 -4.142832599028031d-01,
     1 -1.805991249601928d-01,
     1 -6.009290785739468d-02,
     1 -1.239382725542637d-02,
     1 -8.371529832014113d-04,
     1 8.371529832014113d-04,
     1 1.239382725542637d-02,
     1 6.009290785739468d-02,
     1 1.805991249601928d-01,
     1 4.142832599028031d-01,
     1 7.964747731112430d-01,
     1 1.348993882467059d+00,
     1 2.073471660264395d+00,
     1 2.947904939031494d+00,
     1 3.928129252248612d+00,
     1 4.957203086563112d+00,
     1 5.986360113977494d+00,
     1 6.997957704791519d+00,
     1 7.999888757524622d+00,
     1 8.999998754306120d+00/
      data extraweights16/
     1 1.000007149422537d+00,
     1 1.000395017352309d+00,
     1 1.004798397441514d+00,
     1 1.020308624984610d+00,
     1 1.035167721053657d+00,
     1 1.014359775369075d+00,
     1 9.362411945698647d-01,
     1 8.051212946181061d-01,
     1 6.401489637096768d-01,
     1 4.652220834914617d-01,
     1 3.029123478511309d-01,
     1 1.704889420286369d-01,
     1 7.740135521653088d-02,
     1 2.423621380426338d-02,
     1 3.190919086626234d-03,
     1 3.190919086626234d-03,
     1 2.423621380426338d-02,
     1 7.740135521653088d-02,
     1 1.704889420286369d-01,
     1 3.029123478511309d-01,
     1 4.652220834914617d-01,
     1 6.401489637096768d-01,
     1 8.051212946181061d-01,
     1 9.362411945698647d-01,
     1 1.014359775369075d+00,
     1 1.035167721053657d+00,
     1 1.020308624984610d+00,
     1 1.004798397441514d+00,
     1 1.000395017352309d+00,
     1 1.000007149422537d+00/
c
      ima=(0,1)
c
      if (norder.eq. 0) then
      nskip=1
      nextra=0
      goto 111
      endif
c
      if (norder.eq.4) then
         nskip = 2
         nextra = 6
	 do i = 1,nextra
	    extranodes(i) = extranodes4(i)
	    extraweights(i) = extraweights4(i)
         enddo
         goto 111
      endif
c
      if (norder.eq.8) then 
         nskip = 5
         nextra = 14
	 do i = 1,nextra
	    extranodes(i) = extranodes8(i)
	    extraweights(i) = extraweights8(i)
         enddo
         goto 111
      endif
c
      if (norder.eq.16) then 
         nskip = 10
         nextra = 30
	 do i = 1,nextra
	    extranodes(i) = extranodes16(i)
	    extraweights(i) = extraweights16(i)
         enddo
         goto 111
      endif      

c
c     quadrature code not available for any other order -> return
c
      ier = 1
      call prinf('wrong order for quadrature, norder=*',norder,1)
      stop
      return
111   continue
      ier = 0

c
c     carry out "punctured" trapezoidal rule and fill in matrix
c     entries, skipping entries within nskip of the diagonal
c
      pi = 4.0d0*datan(1.0d0)

c
      na=nskip
      nj=nextra/2
      n=ns-2*na+1

c
      val=0
      iii=ipt-1+na
c
      do 2000 k=0,n-1
      iii=iii+1
      if (iii .gt. ns) iii=iii-ns
      call gfun(eps,zk,mode,xs(ipt),ys(ipt),xs(iii),ys(iii),u)
      val=val+u*dsdt(iii)*h*fs(iii)
 2000 continue
c

      if (norder .eq. 0) return

c
c     now add in the extra points on either side of 
c     the singularity located at ipt, use norder+2 interpolation points
c
      ninterp=norder+2
c
      do 2400 i=1,nextra
      txtra(i)=h*(ipt-1)+h*extranodes(i)
 2400 continue

c
      do 4000 i=1,nextra
c
c     find the closest ninterp points to each of the txtra
c
      n1=txtra(i)/h
      if (txtra(i) .lt. 0) n1=n1-1
      n2=n1+1
      nnn=n1-(ninterp-2)/2

c
      do 2800 j=1,ninterp
      its(j)=nnn+j-1
      its2(j)=its(j)+1
      if (its2(j) .le. 0) its2(j)=its2(j)+ns
      if (its2(j) .gt. ns) its2(j)=its2(j)-ns
 2800 continue

c
c     fill interpolation nodes and function values
c
      do 3000 j=1,ninterp
      tpts(j)=its(j)*h
      xpts(j)=xs(its2(j))
      ypts(j)=ys(its2(j))
      spts(j)=dsdt(its2(j))
      fpts(j)=fs(its2(j))
 3000 continue

c
c     now compute the values of xs, ys, dsdt at ttt using barycentric
c     interpolation
c
      ttt=txtra(i)
      call bary1r(ninterp,tpts,xpts,ttt,xxx)
      call bary1r(ninterp,tpts,ypts,ttt,yyy)
      call bary1r(ninterp,tpts,spts,ttt,sss)
      call bary1c(ninterp,tpts,fpts,ttt,fff)

c
c     evaluate the kernel at the new quadrature point xxx,yyy and
c     add its contribution to val using the alpert weight
c
      call gfun(eps,zk,mode,xs(ipt),ys(ipt),xxx,yyy,u)
      val=val+u*sss*h*fff*extraweights(i)
c
 4000 continue
c
      return
      end
c
c
c
c
c
      subroutine bary1r(n,ts,xs,ttt,xxx)
      implicit real *8 (a-h,o-z)
      real *8 ts(1),xs(1),whts(1000)
c
c     use barycentric interpolation on ts,xs to evaluate at ttt
c     first evaluate the barycentric weights (real routine)
c
c     input:
c       n - the length of ts,xs
c       ts,xs - nodes and function values with which to interpolate
c       ttt - point at which to interpolate
c
c     output:
c       xxx - the interpolated function value at ttt
c
c
      do 1200 i=1,n
      whts(i)=1
 1200 continue
c
      do 1600 i=1,n
      do 1400 j=1,n
      if (i .ne. j) whts(i)=whts(i)/(ts(i)-ts(j))
 1400 continue
 1600 continue
c
c     this uses the '2nd form' of barycentric interpolation, first
c     form the denominator
c
      dd=0
      do 2000 i=1,n
      dd=dd+whts(i)/(ttt-ts(i))
 2000 continue
c
c     now the numerator
c
      dd1=0
      do 2400 i=1,n
      dd1=dd1+whts(i)*xs(i)/(ttt-ts(i))
 2400 continue
c
      xxx=dd1/dd
c
      return
      end
c
c
c
c
c
      subroutine bary1c(n,ts,xs,ttt,xxx)
      implicit real *8 (a-h,o-z)
      real *8 ts(1),whts(1000)
      complex *16 xs(1),xxx,dd1
c
c     use barycentric interpolation on ts,xs to evaluate at ttt
c     first evaluate the barycentric weights (complex routine)
c
c     input:
c       n - the length of ts,xs
c       ts,xs - nodes and function values with which to interpolate
c       ttt - point at which to interpolate
c
c     output:
c       xxx - the interpolated function value at ttt
c
c
      do 1200 i=1,n
      whts(i)=1
 1200 continue
c
      do 1600 i=1,n
      do 1400 j=1,n
      if (i .ne. j) whts(i)=whts(i)/(ts(i)-ts(j))
 1400 continue
 1600 continue
c
c     this uses the '2nd form' of barycentric interpolation,
c     first form the denominator
c
      dd=0
      do 2000 i=1,n
      dd=dd+whts(i)/(ttt-ts(i))
 2000 continue
c
c     now the numerator
c
      dd1=0
      do 2400 i=1,n
      dd1=dd1+whts(i)*xs(i)/(ttt-ts(i))
 2400 continue
c
      xxx=dd1/dd
c
      return
      end
