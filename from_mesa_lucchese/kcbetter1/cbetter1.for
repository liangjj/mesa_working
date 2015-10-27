      program coul
      implicit real*8 (a-h,o-z)
      real*8 mu
      character*8 iword,iw(2)
      parameter (nbig=120000,lbig=#maxltop,nsmall=30000)
      parameter (nchmx=#maxchan,nemx=200)
      parameter (nsmall2=2*nsmall)
      dimension ipow(0:lbig),echan(nchmx),energy(nemx)
      common /integ/ xstart,ystart,ypstart,al,mu,h,xmax
      common /energy/ e,eta,ak
      common/spl/rd,r(nsmall),wwr(nsmall),cs(nsmall),rd26,l
      complex*16 csc(nsmall,0:lbig),scc(nsmall)
      real*8 aj(0:(lbig+1)),ay(0:(lbig+1))
      real*8 aj1(0:(lbig+1)),ay1(0:(lbig+1))
      real*8 csl(nsmall,0:lbig),derl(nsmall,0:lbig),fp(0:lbig)
     $ ,gp(0:lbig),fc(0:lbig),gc(0:lbig),fc1(0:lbig),gc1(0:lbig)
     $ ,fc2(0:lbig),gc2(0:lbig)
      complex*16 yc1,yc2,cf1,cf2
      real*8 css(nsmall),w(nbig,0:lbig)
      data pi/3.14159265358/,iw/"bessel","coulomb"/
      complex*16 c1,c2,ai,ww(nsmall),wwl(nsmall,0:lbig)
      external zero,gg
      dimension der(nsmall),xx(nbig),yy(nbig),zz(nbig)
      dimension scr(nsmall)
      dimension wwlr(nsmall2,0:lbig),cscr(nsmall2,0:lbig)
      equivalence (wwl(1,0),wwlr(1,0)),(csc(1,0),cscr(1,0))
      ai=(0.,1.)
c      call link ("unit5=incoul,unit6=(outcoul,hc,create)
c    1,unit8=cesspl,unit14=bndmat//")
      open(5,file='incoul')
      open(6,file='outcoul')
      ipow(0)=0
      do 66 i=1,lbig
   66 ipow(i)=1
      read(5,101)iword
      znuc=0.
      if(iword.eq.iw(2))then
      open(8,file='cesspl',form='unformatted')
      open(14,file='bndmat',form='unformatted')
      write(8)iword
101   format(a8)
      read(5,*)znuc
      endif
      mu=1.
      read(5,*)lmax,nper,xstart,xmax,nint,alpha
      np=nper*xmax+1
      if(np.gt.nsmall)then
      write(6,100)
100   format(" np.gt.nsmall")
      stop
      endif
      rd=(xmax-xstart)/(np-1)
      rd26=rd*rd/6.
      do 1 i=1,np
1     r(i)=rd*(i-1.)+xstart
      aint=nint
      h=rd/aint
      nl=lmax+1
      if(iword.eq.iw(1))then
      open(8,file='besspl',form='unformatted')
      eta=0.
      go to 666
      endif
      read(14)nchan
      read(14)(echan(i),i=1,nchan)
      ignd=1
      emin=echan(1)
      do 903 i=1,nchan
         if(echan(i).lt.emin) then
            ignd=i
            emin=echan(i)
            endif
 903        continue
c      ignd=ismin(nchan,echan,1)
      write(6,309)(echan(i),i=1,nchan)
309   format(" channel energies",10f10.5)
      read(14)nener
      read(14)(energy(i),i=1,nener)
      do 200 ie=1,nener
      do 200 ic=1,nchan
      e=energy(ie)-echan(ic)+echan(ignd)
      if(e.le.0.0) then
         ak=sqrt(-2.*mu*e)
      else
         ak=sqrt(2.*mu*e)
         endif
      eta=-znuc*mu/ak
      write(8)ak,znuc
666   continue
      if(e.le.0.0) go to 200
      etac=eta
      xlmin=0
      xlmax=nl-1
      xc=xmax
c     call coulcc(xc,etac,zc,nl,fc,gc,fp,gp,sig,12,0,ifail)
      call coulfg(xmax,eta,xlmin,xlmax,fc,gc,fp,gp,2,0,ifail)
      do 33 l=0,lmax
      al=l
      ystart=xstart**(l+1)
      ypstart=(al+1.)*xstart**l
      call outward(xx,yy,index,zero)
      i1=index
      i2=index-1./h
      rat=fc(l)/yy(i1)
      do 2 i=1,index
         x=xx(i)
      w(i,l)=exp(-x*alpha/ak)*yy(i)*rat
2     continue
      j=0
      do 15 i=1,index,nint
      j=j+1
      wwr(j)=w(i,l)
 15   continue
      yp1=(wwr(2)-wwr(1))/rd
      yp2=(wwr(np)-wwr(np-1))/rd
      call spline(r,wwr,np,yp1,yp2,scr,cs)
      call outward(xx,zz,index,gg)
      x1=xx(i1)
      x2=xx(i2)
      det=yy(i1)*zz(i2)-yy(i2)*zz(i1)
c     call coulcc(x1,etac,zl,1,f,f1,fp1,gp1,sig,12,0,ifail)
      if(l.eq.0)
     $ call coulfg(x1,eta,xlmin,xlmax,fc1,gc1,fp,gp,2,0,ifail)
      cf1=gc1(l)+ai*fc1(l)
c     call coulcc(x2,etac,zl,1,f,f2,fp1,gp1,sig,12,0,ifail)
      if(l.eq.0)
     $ call coulfg(x2,eta,xlmin,xlmax,fc2,gc2,fp,gp,2,0,ifail)
      cf2=gc2(l)+ai*fc2(l)
      c1=cf1*zz(i2)-cf2*zz(i1)
      c1=c1/det
      c2=cf2*yy(i1)-cf1*yy(i2)
      c2=c2/det
      j=0
      do 3 i=1,index,nint
      j=j+1
      ww(j)=c1*yy(i)+c2*zz(i)
      ww(j)=ww(j)/xx(i)
      wwl(j,l)=ww(j)
 3    continue
      yc1=(ww(2)-ww(1))/rd
      yc2=(ww(np)-ww(np-1))/rd
      call splinec(r,ww,np,yc1,yc2,scc,csc(1,l))
      c2r=c2
      j=0
      do 5 i=1,index,nint
      j=j+1
      der(j)=w(i,l)*c2r
      der(j)=der(j)/r(j)**ipow(l)
      derl(j,l)=der(j)
 5    continue
      yp1=(der(2)-der(1))/rd
      yp2=(der(np)-der(np-1))/rd
      call spline(r,der,np,yp1,yp2,scr,csl(1,l))
33    continue
      write(8)lmax,np,xstart,rd,alpha
      write(6,*)lmax,np,xstart,rd,alpha
      write(8)(r(i),i=1,np)
      np2=2*np
      write(8)((wwlr(i,k),i=1,np2),k=0,lmax)
      write(8)((derl(i,k),i=1,np),k=0,lmax)
      write(8)((cscr(i,k),i=1,np2),k=0,lmax)
      write(8)((csl(i,k),i=1,np),k=0,lmax)
      if(iword.eq.iw(1))call exit
 200  continue
 1051 format(6e12.3)
      stop
      end
      subroutine splinec(x,y,n,yp1,ypn,scr,y2)
      implicit real*8 (a-h,o-z)
      complex*16 y,yp1,ypn,scr,y2,p,qn,un,sig
      dimension x(n),y(n),scr(n),y2(n)
      y2(1)=-.5
      scr(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      do 11 i=2,n-1
      sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
      p=sig*y2(i-1)+2.
      y2(i)=(sig-1.)/p
11    scr(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     1 /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*scr(i-1))/p
      qn=.5
      un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      y2(n)=(un-qn*scr(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
      y2(k)=y2(k)*y2(k+1)+scr(k)
12    continue
      return
      end
      function zero(x)
      implicit real*8 (a-h,o-z)
      zero=0.
      return
      end
      subroutine outward(xx,yy,index,g)
      implicit real*8 (a-h,o-z)
      real*8 k1,k2,k3
      real*8 mu
      dimension xx(*),yy(*)
      common /integ/ xstart,ystart,ypstart,al,mu,h,xmax
      common /energy/ e,eta,ak
      f(z)=((2.*eta/z-1.)+al*(al+1.)/z/z)*t-g(z)*2.
      tmu=2.*mu
c initialize solution
      itest=0
      h8=h/8.
      h2=h/2.
      yp=ypstart
      y=ystart
      x=xstart
      index=0
c begin integration
    1 t=y
      k1=h*f(x)
      index=index+1
      xx(index)=x
      yy(index)=y
      itest=itest+1
      t=y+h2*yp+h8*k1
      k2=h*f(x+h2)
      t=y+h*yp+h2*k2
      k3=h*f(x+h)
      ypp=yp+k1/6.+2.*k2/3.+k3/6.
      yn=y+h*(yp+(k1+2.*k2)/6.)
      xn=x+h
      if(xn.gt.xmax)go to 23
      yp=ypp
      y=yn
      x=xn
      go to 1
   23 continue
      return
      end
      subroutine spline(x,y,n,yp1,ypn,scr,y2)
      implicit real*8 (a-h,o-z)
      real*8 y,yp1,ypn,scr,y2,p,qn,un,sig
      dimension x(n),y(n),scr(n),y2(n)
      y2(1)=-.5
      scr(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      do 11 i=2,n-1
      sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
      p=sig*y2(i-1)+2.
      y2(i)=(sig-1.)/p
11    scr(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     1 /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*scr(i-1))/p
      qn=.5
      un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      y2(n)=(un-qn*scr(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
      y2(k)=y2(k)*y2(k+1)+scr(k)
12    continue
      return
      end
      function gg(x)
      implicit real*8 (a-h,o-z)
      real*8 mu
      parameter (nsmall=30000)
      common/integ/xstart,ystart,ypstart,al,mu,h,xmax
      common/spl/rd,r(nsmall),ww(nsmall),cs(nsmall),rd26,l
      klo=(x-xstart)/rd
      klo=klo+1
      a=(r(klo+1)-x)/rd
      b=(x-r(klo))/rd
      gg=a*ww(klo)+b*ww(klo+1)+(a*(a*a-1.)*cs(klo)+b*
     1       (b*b-1.)*cs(klo+1))*rd26
      return
      end
      subroutine coulfg(xx,eta1,xlmin,xlmax, fc,gc,fcp,gcp,
     *                  mode1,kfn,ifail)
c
c
c
c          cray version   llnl contact r.l. pexton  x2-4194
c
c          note that routine writes error messages on 6
c
c          reference  journal of computational physics v46 #2 may1982
c                     a.r.barnett
c                     continued fraction evaluation of coulomb functions
c                     f (eta,x)   g (eta,x)  and their derivatives
c                      l           l
c
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c  revised coulomb wavefunction program using steed"s method           c
c                                                                      c
c  a. r. barnett           manchester  march   1981                    c
c                                                                      c
c  original program "rcwfn"      in    cpc  8 (1974) 377-395           c
c                 + "rcwff"      in    cpc 11 (1976) 141-142           c
c  full description of algorithm in    cpc 21 (1981) 297-314           c
c  this version written up       in    cpc xx (1982) yyy-zzz           c
c                                                                      c
c  coulfg returns f,g,f",g", for real xx.gt.0,real eta1 (including 0), c
c   and real lambda(xlmin) .gt. -1 for integer-spaced lambda values    c
c   thus giving positive-energy solutions to the coulomb schrodinger   c
c   equation,to the klein-gordon equation and to suitable forms of     c
c   the dirac equation ,also spherical + cylindrical bessel equations  c
c                                                                      c
c  for a range of lambda values (xlmax - xlmin) must be an integer,    c
c  starting array element is m1 = max0(  int(xlmin+accur),0) + 1       c
c      see text for modifications for integer l-values                 c
c                                                                      c
c  if "mode" = 1  get f,g,f",g"   for integer-spaced lambda values     c
c            = 2      f,g      unused arrays must be dimensioned in    c
c            = 3      f               call to at least length (1)      c
c  if "kfn"  = 0 real        coulomb functions are returned            c
c            = 1 spherical   bessel                                    c
c            = 2 cylindrical bessel                                    c
c  the use of "mode" and "kfn" is independent                          c
c                                                                      c
c  precision#  results to within 2-3 decimals of "machine accuracy"    c
c   in oscillating region x .ge. eta1 + sqrt(eta1**2 + xlm(xlm+1))     c
c   coulfg is coded for real*8 on ibm or equivalent  accur = 10**-16   c
c   use autodbl + extended precision on hx compiler  accur = 10**-33   c
c   for mantissas of 56 + 112 bits. for single precision cdc (48 bits) c
c   reassign dsqrt=sqrt etc.  see text for complex arithmetic version  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit real*8 (a-h,o-z)
c
      dimension    fc(*),gc(*),fcp(*),gcp(*)
      logical      etane0,xlturn
      common       /steed/ paccq,nfp,npq,iexp,m1
c***  common block is for information only.  not required in code
c***  coulfg has calls to#  sqrt, abs, mod ,  int, sign, float, min1
      data zero,one,two,ten2,abort /0.0e0, 1.0e0, 2.0e0, 1.0e2, 2.0e4/
      data half,tm30 / 0.5e0, 1.0e-30 /
      data rt2epi /0.79788 45608 02865e0/
c *** this constant is  dsqrt(two/pi)#  use q0 for ibm real*16# d0 for
c ***  real*8 + cdc double p#  e0 for cdc single p; and truncate value.
c
                        accur = 1.0e-14
c ***            change accur to suit machine and precision required
      mode  = 1
      if(mode1 .eq. 2 .or. mode1 .eq. 3 ) mode = mode1
      ifail = 0
      iexp  = 1
      npq   = 0
      eta   = eta1
      gjwkb = zero
      paccq = one
      if(kfn .ne. 0) eta = zero
                 etane0  = eta .ne. zero
      acc   = accur
      acc4  = acc*ten2*ten2
      acch  = sqrt(acc)
c ***    test range of xx, exit if.le. sqrt(accur) or if negative
c
      if(xx .le. acch)                          go to 100
      x     = xx
      xlm   = xlmin
      if(kfn .eq. 2)  xlm = xlm - half
      if(xlm .le. -one .or. xlmax .lt. xlmin)   go to 105
      e2mm1 = eta*eta + xlm*xlm + xlm
      xlturn= x*(x - two*eta) .lt. xlm*xlm + xlm
      dell  = xlmax - xlmin + acc
      if( abs(mod(dell,one)) .gt. acc) write(6,2040)xlmax,xlmin,dell
      lxtra =   int(dell)
      xll   = xlm +  float(lxtra)
c ***       lxtra is number of additional lambda values to be computed
c ***       xll  is max lambda value, or 0.5 smaller for j,y bessels
c ***         determine starting array element (m1) from xlmin
      m1  = max0(  int(xlmin + acc),0) + 1
      l1  = m1 + lxtra
c
c ***    evaluate cf1  =  f   =  fprime(xl,eta,x)/f(xl,eta,x)
c
      xi  = one/x
      fcl = one
      pk  = xll + one
      px  = pk  + abort
    2 ek  = eta / pk
      f   = (ek + pk*xi)*fcl + (fcl - one)*xi
      pk1 =  pk + one
c ***   test ensures b1 .ne. zero for negative eta; fixup is exact.
             if( abs(eta*x + pk*pk1) .gt. acc)  go to 3
             fcl  = (one + ek*ek)/(one + (eta/pk1)**2)
             pk   =  two + pk
      go to 2
    3 d   =  one/((pk + pk1)*(xi + ek/pk1))
      df  = -fcl*(one + ek*ek)*d
            if(fcl .ne. one )  fcl = -one
            if(d   .lt. zero)  fcl = -fcl
      f   =  f  + df
c
c ***   begin cf1 loop on pk = k = lambda + 1
c
      p     = one
    4 pk    = pk1
        pk1 = pk1 + one
        ek  = eta / pk
        tk  = (pk + pk1)*(xi + ek/pk1)
        d   =  tk - d*(one + ek*ek)
              if( abs(d) .gt. acch)             go to 5
              write (6,1000) d,df,acch,pk,ek,eta,x
              p = p  +   one
              if( p .gt. two )                  go to 110
    5 d     = one/d
              if (d .lt. zero) fcl = -fcl
        df  = df*(d*tk - one)
        f   = f  + df
              if(pk .gt. px)                    go to 110
      if( abs(df) .ge.  abs(f)*acc)             go to 4
                  nfp = pk - xll - 1
      if(lxtra .eq. 0)                          go to 7
c
c *** downward recurrence to lambda = xlm. array gc,if present,stores rl
c
      fcl = fcl*tm30
      fpl = fcl*f
      if(mode .eq. 1) fcp(l1) = fpl
                      fc (l1) = fcl
      xl  = xll
      rl  = one
      el  = zero
      do 6  lp = 1,lxtra
         if(etane0) el = eta/xl
         if(etane0) rl =  sqrt(one + el*el)
         sl    =  el  + xl*xi
         l     =  l1  - lp
         fcl1  = (fcl *sl + fpl)/rl
         fpl   =  fcl1*sl - fcl *rl
         fcl   =  fcl1
         fc(l) =  fcl
         if(mode .eq. 1) fcp(l)  = fpl
         if(mode .ne. 3 .and. etane0) gc(l+1) = rl
    6 xl = xl - one
      if(fcl .eq. zero) fcl = acc
      f  = fpl/fcl
c ***    now we have reached lambda = xlmin = xlm
c ***    evaluate cf2 = p + i.q  again using steed"s algorithm
c ***    see text for compact complex code for sp cdc or non-ansi ibm
c
    7 if( xlturn ) call jwkb(x,eta,max(xlm,zero),fjwkb,gjwkb,iexp)
      if( iexp .gt. 1 .or. gjwkb .gt. one/(acch*ten2))  go to 9
          xlturn = .false.
      ta =  two*abort
      pk =  zero
      wi =  eta + eta
      p  =  zero
      q  =  one - eta*xi
      ar = -e2mm1
      ai =  eta
      br =  two*(x - eta)
      bi =  two
      dr =  br/(br*br + bi*bi)
      di = -bi/(br*br + bi*bi)
      dp = -xi*(ar*di + ai*dr)
      dq =  xi*(ar*dr - ai*di)
    8 p     = p  + dp
         q  = q  + dq
         pk = pk + two
         ar = ar + pk
         ai = ai + wi
         bi = bi + two
         d  = ar*dr - ai*di + br
         di = ai*dr + ar*di + bi
         c  = one/(d*d + di*di)
         dr =  c*d
         di = -c*di
         a  = br*dr - bi*di - one
         b  = bi*dr + br*di
         c  = dp*a  - dq*b
         dq = dp*b  + dq*a
         dp = c
         if(pk .gt. ta)                         go to 120
      if( abs(dp)+ abs(dq).ge.( abs(p)+ abs(q))*acc)   go to 8
                      npq   = pk/two
                      paccq = half*acc/min( abs(q),one)
                      if( abs(p) .gt.  abs(q)) paccq = paccq* abs(p)
c
c *** solve for fcm = f at lambda = xlm,then find norm factor w=w/fcm
c
      gam = (f - p)/q
            if(q .le. acc4* abs(p))             go to 130
      w   = one/ sqrt((f - p)*gam + q)
            go to 10
c *** arrive here if g(xlm) .gt. 10**6 or iexp .gt. 250 + xlturn = .true.
    9 w   = fjwkb
      gam = gjwkb*w
      p   = f
      q   = one
c
c *** normalise for spherical or cylindrical bessel functions
c
   10                     alpha = zero
          if(kfn  .eq. 1) alpha = xi
          if(kfn  .eq. 2) alpha = xi*half
                          beta  = one
          if(kfn  .eq. 1) beta  = xi
          if(kfn  .eq. 2) beta  =  sqrt(xi)*rt2epi
      fcm  =  sign(w,fcl)*beta
           fc(m1)  = fcm
                      if(mode .eq. 3)           go to 11
           if(.not. xlturn)   gcl =  fcm*gam
           if(      xlturn)   gcl =  gjwkb*beta
           if( kfn .ne. 0 )   gcl = -gcl
           gc(m1)  = gcl
           gpl =  gcl*(p - q/gam) - alpha*gcl
                      if(mode .eq. 2)           go to 11
           gcp(m1) = gpl
           fcp(m1) = fcm*(f - alpha)
   11 if(lxtra .eq. 0 ) return
c *** upward recurrence from gc(m1),gcp(m1)  stored value is rl
c *** renormalise fc,fcp at each lambda and correct regular derivative
c ***    xl   = xlm here  and rl = one , el = zero for bessels
         w    = beta*w/ abs(fcl)
         maxl = l1 - 1
      do 12 l = m1,maxl
                      if(mode .eq. 3)           go to 12
                      xl = xl + one
         if(etane0)   el = eta/xl
         if(etane0)   rl = gc(l+1)
                      sl = el + xl*xi
         gcl1     = ((sl - alpha)*gcl - gpl)/rl
         gpl      =   rl*gcl -  (sl + alpha)*gcl1
         gcl      = gcl1
         gc(l+1)  = gcl1
                      if(mode .eq. 2)           go to 12
         gcp(l+1) = gpl
         fcp(l+1) = w*(fcp(l+1) - alpha*fc(l+1))
   12 fc(l+1)     = w* fc(l+1)
      return
 1000 format(/" cf1 accuracy loss# d,df,acch,k,eta/k,eta,x = ",1p7e9.2/)
c
c ***    error messages
c
  100 ifail = -1
      write(6,2000) xx,acch
 2000 format(" for xx = ",1pe12.3," try small-x  solutions",
     *" or x negative"/ ," square root accuracy parameter =  ",e12.3/)
      return
  105 ifail = -2
      write (6,2005) xlmax,xlmin,xlm
 2005 format(/" problem with input order values#xlmax,xlmin,xlm = ",
     *1p3e15.6/)
      return
  110 ifail =  1
      write (6,2010) abort,f ,df,pk,px,acc
 2010 format(" cf1 has failed to converge after ",f10.0," iterations",/
     *" f,df,pk,px,accur =  ",1p5e12.3//)
      return
  120 ifail =  2
      write (6,2020) abort,p,q,dp,dq,acc
 2020 format(" cf2 has failed to converge after ",f7.0," iterations",/
     *" p,q,dp,dq,accur =  ",1p4e17.7,e12.3//)
      return
  130 ifail =  3
      write (6,2030) p,q,acc,dell,lxtra,m1
 2030 format(" final q.le. abs(p)*acc*10**4 , p,q,acc = ",1p3e12.3,4x,
     *" dell,lxtra,m1 = ",e12.3,2i5 /)
      return
 2040 format(" xlmax - xlmin = dell not an integer ",1p3e20.10/)
      end
c
      subroutine jwkb(xx,eta1,xl,fjwkb,gjwkb,iexp)
      implicit real*8 (a-h,o-z)
      real*8      xx,eta1,xl,fjwkb,gjwkb, zero
c *** computes jwkb approximations to coulomb functions    for xl.ge. 0
c *** as modified by biedenharn et al. phys rev 97 (1955) 542-554
c *** calls max,sqrt,log,exp,atan2,float,int        barnett feb 1981
      data   zero,half,one,six,ten/ 0.0e0, 0.5e0, 1.0e0, 6.0e0, 10.0e0 /
      data  dzero, rl35, aloge  /0.0e0, 35.0e0, 0.43429 45 e0 /
      x     = xx
      eta   = eta1
      gh2   = x*(eta + eta - x)
      xll1  = max(xl*xl + xl,dzero)
      if(gh2 + xll1 .le. zero) return
       hll  = xll1 + six/rl35
       hl   = sqrt(hll)
       sl   = eta/hl + hl/x
       rl2  = one + eta*eta/hll
       gh   = sqrt(gh2 + hll)/x
       phi  = x*gh - half*( hl*log((gh + sl)**2/rl2) - log(gh) )
          if(eta .ne. zero) phi = phi - eta*atan2(x*gh,x - eta)
      phi10 = -phi*aloge
      iexp  =  int(phi10)
      if(iexp .gt. 250) gjwkb = ten**(phi10 - float(iexp))
      if(iexp .le. 250) gjwkb = exp(-phi)
      if(iexp .le. 250) iexp  = 0
      fjwkb = half/(gh*gjwkb)
      return
      end
