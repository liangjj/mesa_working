      program asy
*
* $Id: asy.f,v 2.1 94/09/29 18:47:35 vmb Exp Locker: vmb $
*
*
*     main program for FARM
*
      implicit integer (a-z)
      real *8 x, RTC, T1,T2
      parameter (icore = 50000000, iread = 5)
      dimension x(icore), ix(icore)
      common /io/ iin,iout
      equivalence (x(1),ix(1))

      iin = iread
      open(iin, file='farm.inp', access='sequential', status='old')
cz
      t1=RTC()

      call cntrl (x,ix,icore)

      t2=RTC()

      write(iout,'(a,f12.2,a)') ' Total time  ',t2-t1,'  sec'

      end

      integer function iadtwp (i)
*
* $id: iadtwp.f,v 2.1 94/09/29 18:47:29 vmb exp locker: vmb $
*
      implicit integer (a-z)
      iadtwp = ( i + 1 ) / 2
*     iadtwp = i
      return
      end

      integer function wpadti (i)
*
* $id: wpadti.f,v 2.1 94/09/29 18:47:31 vmb exp locker: vmb $
*
      implicit integer (a-z)
      wpadti = i + i - 1
*     wpadti = i
      return
      end

      subroutine asywfn (nchan,nmax,iasy,f,fp,fc,fcp,gc,gcp,eps,ion,
     x                   nfc,ipflg,iasx,xg,xc,xd,j,lj,rafin,
     x                   skj,etj,etaj,jopen,ga1,ga2,xcp,xdp,ierfg)
*
* $Id: asywfn.f,v 2.1 94/09/29 18:43:26 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      logical jopen
      real*8 f,fp,fc,fcp,gc,gcp,eps,xg,xc,xd,rafin,skj,etj,etaj,ga1,
     x       ga2,rscalx,zero,rootk,rho,flj,f1,g1,fp1,gp1,xi,del1,del2,
     x       del3,del4,del,delp,dr1,dpr1,dr2,dpr2,aerr1,aerr2,aepr1,
     x       aepr2,zz,u,uder,xcp,xdp,fp2,gp2,fdd,one,two,dfloat,dr1p,
     x       dr2p,dpr1p,dpr2p,scalx
      parameter ( zero=0.0d0,one=1.0d0,two=2.0d0 )
      dimension fc(nfc),gc(nfc),fcp(nfc),gcp(nfc),ga1(nmax,0:iasy),
     x          f(nchan,nchan,2),fp(nchan,nchan,2),ga2(nmax,0:iasy),
     x          xg(0:iasy),xc(0:iasx),xd(0:iasx),xcp(0:iasx),
     x          xdp(0:iasx)
      common /io/ iin,iout
      common /mach/ scalx,rscalx
*
*     asywfn : asymptotic scattering functions, f, and derivatives, fp.
*
      rootk = sqrt(skj)
      rho = rafin * skj
      flj = dfloat(lj)
      llj = lj * ( lj + 1 )
      if ( jopen ) then
         call coulfg (rho,etaj,flj,flj,fc,gc,fcp,gcp,1,0,ifail)
         if (ifail .ne. 0) then
            write (iout,1000) ifail,lj,rho,etaj
         end if
         f1 = fc(lj+1) / rootk
         g1 = gc(lj+1) / rootk
         fp1 = fcp(lj+1) / rootk
         gp1 = gcp(lj+1) / rootk
         fdd = ( dfloat(llj) / rho + two * etaj ) / rho - one
         fp2 = fdd * f1
         gp2 = fdd * g1
         if ( ipflg .eq. 1 ) then
            write (iout,1000) ifail,lj,rho,etaj
            write (iout,1010) rafin,skj,f1,g1,fp1,gp1
         endif
      else
         zz = dfloat(ion)
         if ( ion .ne. 0 ) then
            call coul (lj,zz,etj,rafin,u,uder,iasy,xg,xc,xd,xcp,xdp,
     x                 ipflg,ierfg)
         else
            call decay (etj,lj,rafin,u,uder)
         endif
         f1 = u / rootk
         fp1 = uder / rootk
         fp2 = ((dfloat(llj) / rho + two * etaj) / rho + one) * skj *
     x         skj * f1
         if ( ipflg .eq. 1 ) write (iout,1020) lj,ion,rafin,etj,skj,f1,
     x                                         fp1
      endif
*     evaluate asymptotic series
      xi = one / rafin
      do 320 i = 1, nchan
         ily = 0
         ilz = 0
         do 220 p = 0, iasy
            if ( abs(ga1(i,p)) .le. scalx * 1.0d-6) ily = p
            if ( abs(ga2(i,p)) .le. scalx * 1.0d-6) ilz = p
  220    continue
         ily1 = ily + 1
         ilym = ily - 1
         if ( ily .eq. iasy ) then
            dr1 = zero
            dpr1 = zero
            aerr1 = zero
            aepr1 = zero
         elseif ( ily .lt. iasy ) then
            ncfct = iasy - ily1
            if ( ipflg .ne. 0 ) write (iout,1030) ily1,iasy
            do 240 ip1 = ily1, iasy
               xg(ip1-ily1) = ga1(i,ip1)
  240       continue
            np = ( ncfct + 1 ) / 2
            call rrfr (ncfct,np,xg,xc,xcp,xd,xdp)
            call rrfvl (np,ily1,xi,xc,xd,xcp,xdp,dr1,dpr1,aerr1,aepr1,
     x                  ierfg)
         endif
         ilz1 = ilz + 1
         ilzm = ilz - 1
         if ( ilz .eq. iasy ) then
            dr2 = zero
            dpr2 = zero
            aerr2 = zero
            aepr2 = zero
         elseif ( ilz .lt. iasy ) then
            ncfct = iasy - ilz1
            if ( ipflg .ne. 0 ) write (iout,1030) ilz1,iasy
            do 245 ip1 = ilz1, iasy
               xg(ip1-ilz1) = ga2(i,ip1)
  245       continue
            np = ( ncfct + 1 ) / 2
            call rrfr (ncfct,np,xg,xc,xcp,xd,xdp)
            call rrfvl (np,ilz1,xi,xc,xd,xcp,xdp,dr2,dpr2,aerr2,aepr2,
     x                  ierfg)
         endif
         if ( ily .ge. 1 ) then
            call rpolq (ga1,xi,dr1p,dpr1p,ilym,nmax,i)
            dr1 = dr1 + dr1p
            dpr1 = dpr1 + dpr1p
         endif
         if ( ilz .ge. 1 ) then
            call rpolq (ga2,xi,dr2p,dpr2p,ilzm,nmax,i)
            dr2 = dr2 + dr2p
            dpr2 = dpr2 + dpr2p
         endif
         dr1 = rscalx * dr1
         dpr1 = rscalx * dpr1
         dr2 = rscalx * dr2
         dpr2 = rscalx * dpr2
*     asymptotic functions and errors
         if ( jopen ) then
            f(i,j,1) =  f1 * dr1 + fp1 * dr2
            fp(i,j,1) =  skj * fp1 * dr1 + f1 * dpr1 + skj *
     x                  fp2 * dr2 + fp1 * dpr2
            f(i,j,2) =  g1 * dr1 + gp1 * dr2
            fp(i,j,2) =  skj * gp1 * dr1 + g1 * dpr1 + skj *
     x                  gp2 * dr2 + gp1 * dpr2
            if ( ierfg .ne. 0 ) then
               aerr1 = rscalx * aerr1
               aerr2 = rscalx * aerr2
               aepr1 = rscalx * aepr1
               aepr2 = rscalx * aepr2
               del1 = abs( f1 * aerr1 + fp1 * aerr2 )
               del2 = abs( skj * fp1 * aerr1 + f1 * aepr1 +
     x                skj * fp2 * aerr2 + fp1 * aepr2 )
               del = max( del1, del2 )
               if ( del .gt. eps .or. ipflg .ne. 0 ) then
                  call galerr (i,j,f(i,j,1),del1,fp(i,j,1),del2)
               endif
               del3 = abs( g1 * aerr1 + gp1 * aerr2 )
               del4 = abs( skj * gp1 * aerr1 + g1 * aepr1 +
     x                skj * gp2 * aerr2 + gp1 * aepr2 )
               delp = max( del3, del4 )
               if ( delp .gt. eps .or. ipflg .ne. 0 ) then
                  call galerr (i,j,f(i,j,2),del3,fp(i,j,2),del4)
               endif
            endif
         else
            f(i,j,1) = f1 * dr1 - fp1 * dr2 / skj
            fp(i,j,1) = fp1 * dr1 + f1 * dpr1 - ( fp2 *
     x                  dr2 + fp1 * dpr2 ) / skj
            if ( ierfg .ne. 0 ) then
               aerr1 = rscalx * aerr1
               aerr2 = rscalx * aerr2
               aepr1 = rscalx * aepr1
               aepr2 = rscalx * aepr2
               del1 = abs( f1 * aerr1 - fp1 * aerr2 / skj )
               del2 = abs( fp1 * aerr1 + f1 * aepr1 - ( fp2 *
     x                aerr2 + fp1 * aepr2 ) / skj )
               del = max( del1, del2 )
               if ( del .gt. eps .or. ipflg .ne. 0 ) then
                  call galerr (i,j,f(i,j,1),del1,fp(i,j,1),del2)
               endif
            endif
         endif
  320 continue
      return
*
 1000 format (' coulfg : ifail =',i4,' lj =',i4,' rho =',d16.8,
     x        ' etaj =',d16.8)
 1010 format (' rafin =',f12.5,' skj =',d16.8,' f =',d16.8,' g =',d16.8,
     x        /,' fp =',d16.8,' gp =',d16.8)
 1020 format (' coul/decay : lj =',i3,' ion =',i3,' rafin =',d16.8,
     x        ' etj =',d16.8,' skj =',f12.5,/,14x,' f1 =',d16.8,
     x        ' dp1 =',d16.8)
 1030 format (' coefficients',i3,' to',i3,' extrapolated ')
      end
      subroutine couerr (fx,efx,eps,fxp,efxp)
*
* $Id: couerr.f,v 2.1 94/09/29 18:43:28 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 fx,efx,eps,fxp,efxp,escat
      common /ie/ escat,ierr
*
      if ( ierr .ne. 0 ) then
         write (ierr,1000) escat,fx,efx,eps,fxp,efxp
      endif
      return
 1000 format (' coul : e =',f12.5,' fx  =',d16.8,' efx  =',d16.8,
     x        ' eps =',d16.8,/,'         fxp =',d30.17,' efxp =',d16.8)
      end
      subroutine coul (l,z,e,r,fx,fxp,nterm,xg,xc,xd,xcp,xdp,
     x                 ipflg,ierr)
*
* $Id: coul.f,v 2.1 94/09/29 18:43:29 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 z,e,r,fx,fxp,xg,xc,xd,one,half,two,eps,ak,c,x,cx,dfloat,
     x       fl,fl1,c1,q,f1,f2,f3,an,efx,fac,f,ef,epf,efxp,xcp,xdp,fp,
     x       xi
      parameter ( one=1.0d0,half=0.5d0,two=2.0d0,eps=1.0d-08 )
      dimension xg(0:nterm),xc(*),xd(*),xcp(*),xdp(*)
      common /io/ iin,iout
*
*     coul : whittaker function and derivative using the asymptotic
*            expansion converted to a continued fraction
*
      ak = sqrt(-e)
      c = ak / z
      x = z * r
      cx = ak * r
      fl = dfloat(l)
      fl1 = dfloat(l+1)
      c1 = one / c
      q = one
      xg(0) = q
      f1 = fl - c1
      f2 = - fl1 - c1
      f3 = - half * c1
      do 10 n = 1, nterm
         an = dfloat(n)
         q = ( f1 + an ) * ( f2 + an ) * f3 * q / an
         xg(n) = q
   10 continue
      if ( ipflg .ne. 0 ) write (iout,1000) nterm,(xg(i),i = 0,nterm)
*
*     continued fraction
*
      np = ( nterm + 1 ) / 2
      call rrfr (nterm,np,xg,xc,xcp,xd,xdp)
      xi = one / x
      call rrfvl (np,0,xi,xc,xd,xcp,xdp,f,fp,ef,epf,ierr)
      fac = exp(-cx) * ( two * cx ) ** c1
      fx = fac * f
      fxp = z * ( ( one / cx - c ) * fx + fac * f p)
      if ( ipflg .ne. 0 ) then
         write (iout,1010) fx,fxp
      end if
      if ( ierr .ne. 0 ) then
         efx = abs(ef*fx)
         efxp = abs(epf*fxp)
         if ( efx .gt. eps .or. efxp .gt. eps ) then
             write (iout,1020) e,r,fx,efx,eps,fxp,efxp
c            call couerr (fx,efx,eps,fxp,efxp)
         endif
      endif
      return
*
 1000 format (' nterm =',i6,/,(4d16.8))
 1010 format (' coul : fx  =',d16.8,'fxp =',d16.8)
 1020 format (' coul : e =',f12.5,'  r =', f12.5,/,
     x        '        fx  =',d16.8,' efx  =',d16.8,' eps =',d16.8,/,
     x        '        fxp =',d30.17,' efxp =',d16.8)
      end
      subroutine coulfg (xx,eta1,xlmin,xlmax, fc,gc,fcp,gcp,mode1,kfn,
     x                   ifail)
*
* $Id: coulfg.f,v 2.1 94/09/29 18:43:31 vmb Exp Locker: vmb $
*
*
************************************************************************
*                                                                      *
*     coulomb wavefunction program using steed's method                *
*                                                                      *
*  a. r. barnett           manchester  march   1981                    *
*                                                                      *
*  original program 'rcwfn'      in    cpc  8 (1974) 377-395           *
*                 + 'rcwff'      in    cpc 11 (1976) 141-142           *
*  full description of algorithm in    cpc 21 (1981) 297-314           *
*  this version written up       in    cpc 27 (1982) 147-166           *
*                                                                      *
*  coulfg returns f,g,f',g', for real xx.gt.0,real eta1 (including 0), *
*   and real lambda(xlmin) .gt. -1 for integer-spaced lambda values    *
*   thus giving positive-energy solutions to the coulomb schrodinger   *
*   equation,to the klein-gordon equation and to suitable forms of     *
*   the dirac equation ,also spherical & cylindrical bessel equations  *
*                                                                      *
*  for a range of lambda values (xlmax - xlmin) must be an integer,    *
*  starting array element is m1 = max0(idint(xlmin+accur),0) + 1       *
*      see text for modifications for integer l-values                 *
*                                                                      *
*  if 'mode' = 1  get f,g,f',g'   for integer-spaced lambda values     *
*            = 2      f,g      unused arrays must be dimensioned in    *
*            = 3      f               call to at least length (1)      *
*  if 'kfn'  = 0 real        coulomb functions are returned            *
*            = 1 spherical   bessel      '      '     '                *
*            = 2 cylindrical bessel      '      '     '                *
*  the use of 'mode' and 'kfn' is independent                          *
*                                                                      *
*  precision:  results to within 2-3 decimals of 'machine accuracy'    *
*   in oscillating region x .ge. eta1 + sqrt(eta1**2 + xlm(xlm+1))     *
*   coulfg is coded for real*8 on ibm or equivalent  accur = 10**-16   *
*   use autodbl + extended precision on hx compiler  accur = 10**-33   *
*   for mantissas of 56 & 112 bits. for single precision cdc (48 bits) *
*   reassign dsqrt=sqrt etc.  see text for complex arithmetic version  *
************************************************************************
*
      implicit real*8 (a-h,o-z)
      dimension    fc(1),gc(1),fcp(1),gcp(1)
      logical      etane0,xlturn
      parameter ( zero=0.0d0,one=1.0d0,two=2.0d0,ten2=1.0d2,abort=8.0d4,
     x            half=0.5d0,tm30=1.0d-30,ten=1.0d1,
     x            rt2dpi=0.79788 45608 02865 35587 98921 19868 76373 d0)
      common /io/ iin,iout
*
*     constant is  dsqrt(two/pi)
*
      accur = 1.0d-15
*
      mode  = 1
      if ( mode1 .eq. 2 .or. mode1 .eq. 3 ) mode = mode1
      ifail = 0
      iexp  = 1
      eta   = eta1
      gjwkb = zero
      paccq = one
      if ( kfn .ne. 0 ) eta = zero
      etane0  = eta .ne. zero
      acc   = accur * ten
      acc4  = acc*ten2*ten2
      acch  = sqrt(acc)
*
*        test range of xx, exit if .le. sqrt(accur) or if negative
*
      if ( xx .le. acch ) then
         ifail = - 1
         write (iout,1000) xx,acch
         return
      endif
      x     = xx
      xlm   = xlmin
      if ( kfn .eq. 2 ) xlm = xlm - half
      if ( xlm .le. -one .or. xlmax .lt. xlmin ) then
         ifail = - 2
         write (iout,1010) xlmax,xlmin,xlm
         return
      endif
      e2mm1 = eta * eta + xlm * xlm + xlm
      xlturn = x * ( x - two * eta ) .lt. xlm * xlm + xlm
      dell  = xlmax - xlmin + acc
      if ( abs(mod(dell,one)) .gt. 10*acc ) write (iout,1050) xlmax,xlmin,
     x                                                      dell
      lxtra = int(dell)
      xll   = xlm + dfloat(lxtra)
*
*      lxtra = number of additional lambda values to be computed
*      xll   = max lambda value, or 0.5 smaller for j,y bessels
*              determine starting array element (m1) from xlmin
*
      m1  = max(int(xlmin + acc),0) + 1
      l1  = m1 + lxtra
*
*     evaluate cf1  =  f   =  fprime(xl,eta,x) / f(xl,eta,x)
*
      xi  = one / x
      fcl = one
      pk  = xll + one
      px  = pk  + abort
      f   =  eta / pk + pk * xi
      if ( abs(f) .lt. tm30 ) f = tm30
      d = zero
      c = f
*
*     begin cf1 loop on pk = k = lambda + 1
*
   10 pk1   = pk + one
      ek  = eta / pk
      rk2 = one + ek * ek
      tk  = ( pk + pk1 ) * ( xi + ek / pk1 )
      d   =  tk - rk2 * d
      c   =  tk - rk2 / c
      if ( abs(c) .lt. tm30 ) c = tm30
      if ( abs(d) .lt. tm30 ) d = tm30
      d = one / d
      df = d * c
      f  = f * df
      if ( d .lt. zero ) fcl = - fcl
      pk = pk1
      if ( pk .gt. px ) then
         ifail = 1
         write (iout,1020) abort,f,df,pk,px,acc
         return
      endif
      if ( abs(df-one) .ge. acc ) go to 10
      if ( lxtra .eq. 0 ) go to 30
*
*     downward recurrence to lambda = xlm. array gc,if present,stores rl
*
      fcl = fcl * tm30
      fpl = fcl * f
      if ( mode .eq. 1 ) fcp(l1) = fpl
      fc (l1) = fcl
      xl  = xll
      rl  = one
      el  = zero
      do 20  lp = 1, lxtra
         if (etane0) el = eta / xl
         if (etane0) rl = sqrt(one + el*el)
         sl    =  el  + xl * xi
         l     =  l1  - lp
         fcl1  = ( fcl * sl + fpl ) / rl
         fpl   =  fcl1 * sl - fcl * rl
         fcl   =  fcl1
         fc(l) =  fcl
         if ( mode .eq. 1 ) fcp(l)  = fpl
         if ( mode .ne. 3 .and. etane0 ) gc(l+1) = rl
         xl = xl - one
   20 continue
      if ( fcl .eq. zero ) fcl = acc
      f  = fpl / fcl
*
*     now reached lambda = xlmin = xlm
*     evaluate cf2 = p + i.q  again using steed's algorithm
*
   30 if ( xlturn ) call jwkb (x,eta,max(xlm,zero),fjwkb,gjwkb,iexp)
      if ( iexp .gt. 1 .or. gjwkb .gt. one/(acch*ten2) )  go to 50
      xlturn = .false.
      ta =  two * abort
      pk =  zero
      wi =  eta + eta
      p  =  zero
      q  =  one - eta * xi
      ar = - e2mm1
      ai =  eta
      br =  two * ( x - eta )
      bi =  two
      dr =  br / ( br * br + bi * bi )
      di = - bi / ( br * br + bi * bi )
      dp = - xi * ( ar * di + ai * dr )
      dq =  xi * ( ar * dr - ai * di )
   40 p     = p  + dp
      q  = q  + dq
      pk = pk + two
      ar = ar + pk
      ai = ai + wi
      bi = bi + two
      d  = ar * dr - ai * di + br
      di = ai * dr + ar * di + bi
      c  = one / ( d * d + di * di )
      dr =  c * d
      di = - c * di
      a  = br * dr - bi * di - one
      b  = bi * dr + br * di
      c  = dp * a - dq * b
      dq = dp * b + dq * a
      dp = c
      if ( pk .gt. ta ) then
         ifail = 2
         write (iout,1030) abort,p,q,dp,dq,acc
         return
      endif
      if ( abs(dp)+abs(dq) .ge. (abs(p)+abs(q))*acc )   go to 40
      paccq = half * acc / min( abs(q), one )
      if ( abs(p) .gt. abs(q) ) paccq = paccq * abs(p)
*
*     solve for fcm = f at lambda = xlm, then find norm factor w = w/fcm
*
      gam = ( f - p ) / q
      if ( q .le. acc4*abs(p) ) then
         ifail = 3
         write (iout,1040) p,q,acc,dell,lxtra,m1
         return
      endif
      w   = one / sqrt((f - p)*gam + q)
      go to 60
*
*     here if g(xlm) .gt. 10**6 or iexp .gt. 70 and xlturn = .true.
*
   50 w   = fjwkb
      gam = gjwkb * w
      p   = f
      q   = one
*
*     normalise for spherical or cylindrical bessel functions
*
   60 alpha = zero
      if ( kfn  .eq. 1 ) alpha = xi
      if ( kfn  .eq. 2 ) alpha = xi * half
      beta  = one
      if ( kfn  .eq. 1 ) beta  = xi
      if ( kfn  .eq. 2 ) beta  = sqrt(xi) * rt2dpi
      fcm  = sign(w,fcl) * beta
      fc(m1)  = fcm
      if ( mode .eq. 3 ) go to 70
      if ( .not. xlturn ) gcl =  fcm * gam
      if ( xlturn ) gcl =  gjwkb * beta
      if ( kfn .ne. 0 ) gcl = - gcl
      gc(m1) = gcl
      gpl =  gcl * ( p - q / gam ) - alpha * gcl
      if ( mode .eq. 2 ) go to 70
      gcp(m1) = gpl
      fcp(m1) = fcm * ( f - alpha )
   70 if ( lxtra .eq. 0 ) return
*
*     upward recurrence from gc(m1),gcp(m1)  stored value is rl
*     renormalise fc,fcp at each lambda and correct regular derivative
*     xl   = xlm here and rl = one , el = zero for bessels
*
      w    = beta * w / abs(fcl)
      maxl = l1 - 1
      do 80 l = m1, maxl
         if ( mode .eq. 3 ) go to 80
         xl = xl + one
         if (etane0) el = eta / xl
         if (etane0) rl = gc(l+1)
         sl = el + xl * xi
         gcl1 = ( ( sl - alpha ) * gcl - gpl ) / rl
         gpl = rl * gcl - ( sl + alpha ) * gcl1
         gcl = gcl1
         gc(l+1) = gcl1
         if ( mode .eq. 2 ) go to 80
         gcp(l+1) = gpl
         fcp(l+1) = w * ( fcp(l+1) - alpha * fc(l+1) )
         fc(l+1) = w * fc(l+1)
   80 continue
      return
*
 1000 format (' for xx = ',1pd12.3,' try small-x  solutions',
     x        ' or x negative',/,' square root accuracy parameter =  ',
     x        d12.3/)
 1010 format (/,' problem with input order values:xlmax,xlmin,xlm = ',
     x        1p3d15.6,/)
 1020 format (' cf1 has failed to converge after ',f10.0,' iterations',
     x        /,' f,df,pk,px,accur =  ',1p5d12.3,//)
 1030 format (' cf2 has failed to converge after ',f7.0,' iterations',
     x        /,' p,q,dp,dq,accur =  ',1p4d17.7,d12.3//)
 1040 format (' final q.le.abs(p)*acc*10**4 , p,q,acc = ',1p3d12.3,4x,
     x        ' dell,lxtra,m1 = ',d12.3,2i5,/)
 1050 format (' xlmax - xlmin = dell not an integer ',1p3d20.10,/)
      end
      subroutine decay (skj,l,r,u,uder)
*
* $Id: decay.f,v 2.1 94/09/29 18:43:33 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 skj,r,u,uder,zero,one,two,a,sk,rootk,twor,factor,xlp,sum1,
     x       sum2,dfloat
      parameter ( zero=0.0d0,one=1.0d0,two=2.0d0,lmax=150 )
      dimension a(lmax+1)
      common /io/ iin,iout
*
*     decay : exponentially decaying solution at r of the equation
*               2     2    2
*              d u /dr +(-k-l(l+1)/(r*r))u =0
*                 l                       l
*             where k*k is real positive, l = angular momentum
*             solution obtained by expanding,
*             u = exp(-k*r) * ( a(1) + a(2)/r + a(3)/r**2 + ... )
*             and using recursion relation
*             a    = ( (l*l+l-n*n-n) / ( 2*k*(n+1) ) ) *a  with a =1
*              n+1                                       n       0
*             series terminates when n=l
*
      sk = - skj
      if ( sk .lt. zero ) then
          write (iout,1000) sk
          stop
      endif
      if ( l .gt. lmax ) then
         write (iout,1010) l,lmax
         stop
      endif
      rootk = sqrt(sk)
      twor = two * rootk
      factor = exp(-rootk*r)
      if ( l .gt. 0 ) then
         a(1) = one
         xlp = dfloat(l*(l+1))
         n = l
         do 10 j = 1, n
            a(j+1) = a(j) * ( xlp - dfloat((j-1)*j)) / ( twor *
     x               dfloat(j) )
   10    continue
         sum1 = one
         sum2 = zero
         do 20 jj = 1, n
            sum1 = sum1 + a(jj+1) / ( r ** jj )
            sum2 = sum2 + ( a(jj+1) * dfloat(jj) ) / ( r ** (jj+1) )
   20    continue
         u = factor * sum1
         uder = - factor * sum2 - rootk * u
      else
         u = factor
         uder = - rootk * u
      endif
      return
*
 1000 format (' decay : input error, sk =',f12.5)
 1010 format (' decay : input error, l  =',i10,' lmax =',i10)
      end
      function flminx ()
*
* $Id: flminx.f,v 2.1 94/09/29 18:43:35 vmb Exp Locker: vmb $
*
      real*8 flminx,value
      parameter ( value=1.0d-160 )
      flminx = value
      return
      end
      subroutine gail2 (j,skj,etaj,et,l2p,nchan,cf,numax,jopen,iasy,
     x                  degeny,ga,gb,gmax,fac1,fac2,fac3,fac4,
     x                  iac1,deg,sma,smb)
*
* $Id: gail2.f,v 2.1 94/09/29 18:43:36 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      logical jopen
      real*8 skj,etaj,et,cf,degeny,zero,etj,ga,gb,sma,smb,sgn,
     x       fac1,fac2,fac3,fac4,save,two,dfloat,dena,gmax,one,fp1,fp2,
     x       fp12,scalx,rscalx
      parameter ( zero=0.0d0,one=1.0d0,two=2.0d0 )
      dimension et(nchan),l2p(nchan),cf(nchan,nchan,numax),deg(nchan),
     x          ga(nchan,0:iasy+1),gb(nchan,0:iasy+1),fac1(nchan),
     x          fac2(nchan),fac3(nchan),fac4(nchan),iac1(nchan),
     x          sma(nchan),smb(nchan)
      common /io/ iin,iout
      common /mach/ scalx,rscalx
*
      iasy1 = iasy + 1
      etj = et(j)
      lj = l2p(j)
      llj = lj * ( lj + 1 )
      sgn = one
      if (.not.jopen) sgn = - sgn
      do 800 i = 1, nchan
            fac2(i) = two * skj
            iac1(i) = llj - l2p(i) * ( l2p(i) + 1 )
            fac3(i) = two * etaj
            fac4(i) = two * dfloat( llj ) / skj
            ga(i,0) = zero
            ga(i,1) = zero
            gb(i,0) = zero
            gb(i,2) = zero
         if ( abs(et(i)-etj) .gt. degeny ) then
            fac1(i) = one / ( et(i) - etj )
            deg(i) = 0
            ga(i,2) = scalx * fac1(i) * cf(i,j,1)
            gb(i,1) = zero
         else
            deg(i) = 1
            gb(i,1) = scalx * cf(i,j,1) / fac2(i)
         endif
  800 continue
      ga(j,0) = scalx
      do 500 p = 3, iasy1
         p1 = p - 1
         p2 = p - 2
         p12 = p1 * p2
         fp1 = dfloat(p1)
         fp2 = dfloat(p2)
         fp12 = dfloat(p1 + p2)
         do 300 m = 1, nchan
            sma(m) = zero
            smb(m) = zero
            save = dfloat( p12 + iac1(m) )
            gb(m,p) = - save * gb(m,p2)
            ga(m,p) =  - save * ga(m,p2) + fac3(m) * sgn *
     x                  fp12 * gb(m,p2) + sgn * fac4(m) *
     x                  fp2 * gb(m,p-3)
  300    continue
         if (numax .ge. 1) then
            nux = min( numax, p1 )
            do 35 nu = 1,nux
               call dgemv('n',nchan,nchan,one,cf(1,1,nu),nchan,
     x                    ga(1,p1-nu),1,one,sma,1)
               call dgemv('n',nchan,nchan,one,cf(1,1,nu),nchan,
     x                    gb(1,p1-nu),1,one,smb,1)
 35         continue
         end if
         do 45 m = 1,nchan
            if ( deg(m) .eq. 1 ) then
               dena = one / (fac2(m) * fp1)
               ga(m,p1) = -sgn * (gb(m,p) + smb(m)) * dena
               gb(m,p1) = (ga(m,p) + sma(m)) * dena
            else
               dena = fac2(m) * fp1
               ga(m,p) = (ga(m,p) + sma(m) - dena * gb(m,p1))*fac1(m)
               gb(m,p) = (gb(m,p) + smb(m) + sgn * dena *
     x                   ga(m,p1)) * fac1(m)
            end if
             gmax = max(abs(ga(m,p1)),abs(gb(m,p1)),gmax)
 45      continue
         if ( gmax .gt. rscalx ) then
            p2 = (p2/2) * 2
            write (iout,1050) p2
            stop
         endif
  500 continue
      return
*
 1050 format (' gail2 : scaling overflow, set iasy =',i3)
      end
      subroutine gailco (nchan,ipflg,et,l2p,ion,cf,lamax,
     x                   iasy,degeny,gmax,x,ix)
*
* $Id: gailco.f,v 2.1 94/09/29 18:43:38 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      logical jopen
      real*8 et,cf,degeny,x,zero,etj,skj,
     x       etaj,dfloat,gmax
      parameter ( zero=0.0d0 )
      dimension et(nchan),ipflg(4),
     x          l2p(nchan),cf(nchan,nchan,lamax),x(*),ix(*)
      common /glpt/ ga1,ga2,fc,fcp,gc,gcp,xg,xc,xd,xcp,xdp,fac1,fac2,
     x              fac3,fac4,iac1,deg,sma,smb
*
      iasy1 = iasy - 1
      ncha = nchan * (iasy + 1)
      ga = ga1 - ncha
      gb = ga2 - ncha
      gmax = zero
      do 80 j = 1, nchan
         etj = et(j)
         skj = sqrt(abs(etj))
         etaj = - dfloat(ion) / skj
         jopen = etj .ge. zero
         ga = ga + ncha
         gb = gb + ncha
         call gail2 (j,skj,etaj,et,l2p,nchan,cf,lamax,jopen,iasy1,
     x               degeny,x(ga),x(gb),gmax,x(fac1),x(fac2),
     x               x(fac3),x(fac4),ix(iac1),ix(deg),x(sma),x(smb))
   80 continue
*
*     print expansion coefficients if ipflg=1
*
      if (ipflg(1) .ne. 0) then
         ga = ga1 - ncha
         gb = ga2 - ncha
         do 90 j = 1, nchan
            ga = ga + ncha
            gb = gb + ncha
            call prgail (j,nchan,iasy1,x(ga),x(gb))
 90      continue
      end if
      return
      end
      subroutine gailit (f,fp,nchan,nmax,ipflg,rafin,et,l2p,ion,
     x                   iasy,eps,lmax,x)
*
* $Id: gailit.f,v 2.1 94/09/29 18:43:40 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      logical jopen
      real*8 f,fp,rafin,et,eps,x,zero,etj,skj,
     x       etaj,dfloat
      parameter ( zero=0.0d0 )
      dimension f(nchan,nchan,2),fp(nchan,nchan,2),et(nmax),ipflg(4),
     x          l2p(nmax),x(*)
      common /glpt/ ga1,ga2,fc,fcp,gc,gcp,xg,xc,xd,xcp,xdp,fac1,fac2,
     x              fac3,fac4,iac1,deg,sma,smb
*
      nfc = lmax + 10
      iasx = iasy / 2
      iasy1 = iasy - 1
      do 50 jdd = 1, nchan
         do 40 idd = 1, nchan
            f(idd,jdd,1) = zero
            fp(idd,jdd,1) = zero
            f(idd,jdd,2) = zero
            fp(idd,jdd,2) = zero
   40    continue
   50 continue
      ncha = nmax * (iasy + 1)
      ga = ga1 - ncha
      gb = ga2 - ncha
      do 80 j = 1, nchan
         jchan = j
         lj = l2p(j)
         etj = et(j)
         skj = sqrt(abs(etj))
         etaj = - dfloat(ion) / skj
         jopen = etj .ge. zero
         ga = ga + ncha
         gb = gb + ncha
         call asywfn (nchan,nmax,iasy1,f,fp,x(fc),x(fcp),x(gc),
     x                x(gcp),eps,ion,nfc,ipflg(3),iasx,x(xg),
     x                x(xc),x(xd),jchan,lj,rafin,skj,etj,etaj,jopen,
     x                x(ga),x(gb),x(xcp),x(xdp),ipflg(4))
   80 continue
      return
      end
      subroutine galerr (i,j,f1,d1,fp1,dp1)
*
* $Id: galerr.f,v 2.1 94/09/29 18:43:42 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 f1,fp1,d1,dp1,escat,e
      common /ie/ escat,ierr
*
      if ( ierr .ne. 0 ) then
         e = escat
         write (ierr,1000) e,i,j,f1,d1,fp1,dp1
      endif
      return
 1000 format (' galerr : e =',f12.5,' i =',i3,' j =',i3,' f =',f12.6,
     x        ' (',f12.6,') fp =',f12.6,' (',f12.6,')')
      end
      subroutine galptr (nchan,iasy,lmax,pt0,mcor)
*
* $Id: galptr.f,v 2.1 94/09/29 18:43:43 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      common /io/ iin,iout
      common /glpt/ ga1,ga2,fc,fcp,gc,gcp,xg,xc,xd,xcp,xdp,fac1,fac2,
     x              fac3,fac4,iac1,deg,sma,smb
*
      nfc = lmax + 10
      iasx = iasy / 2
      iasx1 = iasx + 1
      iasy1 = iasy + 1
      lega =iasy1  * nchan * nchan
*
      ga1 = pt0
      ga2 = ga1 + lega
      fc = ga2 + lega
      fcp = fc + nfc
      gc = fcp + nfc
      gcp = gc + nfc
      xg = gcp + nfc
      xc = xg + iasy1
      xd = xc + iasx1
      xcp = xd + iasy1
      xdp = xcp + iasx1
      fac1 = xdp + iasx1
      fac2 = fac1 + nchan
      fac3 = fac2 + nchan
      fac4 = fac3 + nchan
      sma = fac4 + nchan
      smb = sma + nchan
      iac1 = wpadti(smb + nchan + 1)
      deg = iac1 + nchan
      last = iadtwp(deg + nchan + 1) + 1
      if ( last .gt. mcor ) then
         write (iout,1000) last,mcor
         stop
      endif
      pt0 = last
      return
 1000 format (' storage overflow ( galptr ) : last =',i10,' core =',
     x        i10)
      end
      subroutine jwkb (xx,eta1,xl,fjwkb,gjwkb,iexp)
*
* $Id: jwkb.f,v 2.1 94/09/29 18:43:45 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 aloge,half,one,rl35,six,ten,zero,eta,eta1,fjwkb,gh,gh2,
     x       gjwkb,hl,hll,phi,phi10,rl2,sl,x,xl,xll1,xx,dfloat
      parameter ( zero=0.0d0,half=0.5d0,one=1.0d0,six=6.0d0,ten=1.0d1,
     x            rl35=35.0d0,aloge=0.4342945d0 )
*
*     jwkb approximations to coulomb functions for xl.ge. 0
*          biedenharn et al. phys rev 97 (1955) 542-554
*
      x     = xx
      eta   = eta1
      gh2   = x * ( eta + eta - x )
      xll1  = max(xl*xl + xl,zero)
      if ( gh2 + xll1 .le. zero ) return
      hll  = xll1 + six / rl35
      hl   = sqrt(hll)
      sl   = eta / hl + hl / x
      rl2  = one + eta * eta / hll
      gh   = sqrt(gh2 + hll) / x
      phi  = x * gh - half * ( hl * log((gh + sl)**2/rl2) - log(gh) )
      if ( eta .ne. zero ) phi = phi - eta * atan2(x*gh,x - eta)
      phi10 = - phi * aloge
      iexp  =  int(phi10)
      if ( iexp .gt. 70 ) gjwkb = ten ** ( phi10 - dfloat(iexp) )
      if ( iexp .le. 70 ) gjwkb = exp(-phi)
      if ( iexp .le. 70 ) iexp  = 0
      fjwkb = half / ( gh * gjwkb )
      return
      end
      subroutine prgail (j,nchan,iasy,ga,gb)
*
* $Id: prgail.f,v 2.1 94/09/29 18:43:55 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 ga,gb
      dimension ga(nchan,0:iasy+1),gb(nchan,0:iasy+1)
      common /io/ iin,iout
*
*     print expansion coefficients if ipflg=1
*
      write (iout,1010) j
      write (iout,1020)
      write (iout,1030)
      do 260 i = 0, iasy
         write (iout,1000) (ga(jj,i),jj = 1,nchan)
  260 continue
      write (iout,1040)
      do 270 i = 0, iasy
         write (iout,1000) (gb(jj,i),jj = 1,nchan)
  270 continue
*
 1000 format (1x,5d14.6)
 1010 format (' expansion coefficients with j= ',i2)
 1020 format (1x,33('-'))
 1030 format (' gailitis a-coefficients :')
 1040 format (' gailitis b-coefficients :')
      return
      end
      subroutine rpolp (cr,zr,sr,srp,n)
*
* $Id: rpolp.f,v 2.1 94/09/29 18:43:46 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 cr,zr,sr,srp,dfloat
      dimension cr(0:n)
*
*     evaluate power series and derivative, coefficients cr,
*     argument 1/zr, results in sr and srp
*
      sr = cr(n)
      srp = dfloat(n) * sr
      do 10 i = n-1, 0, -1
         sr = sr * zr + cr(i)
         srp = srp * zr + dfloat(i) * cr(i)
   10 continue
      srp = - srp * zr
      return
      end
      subroutine rpolpp (cr,zr,sr,srp,n,np)
*
* $Id: rpolpp.f,v 2.1 94/09/29 18:43:48 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 cr,zr,sr,srp,dfloat,srq
      dimension cr(0:n)
*
*     evaluate power series and derivative, coefficients cr,
*     argument 1/zr, results in sr and srp
*     leading inverse power of np
*
      sr = cr(n)
      srp = dfloat(np+n) * sr
      do 10 i = n-1, 0, -1
         sr = sr * zr + cr(i)
         srp = srp * zr + dfloat(np+i) * cr(i)
   10 continue
      srq = zr**np
      sr = srq * sr
      srp = - srq * srp * zr
      return
      end
      subroutine rpolq (cr,zr,sr,srp,n,m,m1)
*
* $Id: rpolq.f,v 2.1 94/09/29 18:43:49 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 cr,zr,sr,srp,dfloat
      dimension cr(m,0:n)
*
*     evaluate power series and derivative, coefficients cr,
*     argument zr, results in sr and srp
*
      sr = cr(m1,n)
      srp = dfloat(n) * sr
      do 10 i = n-1, 0, -1
         sr = sr * zr + cr(m1,i)
         srp = srp * zr + dfloat(i) * cr(m1,i)
   10 continue
      srp = - zr * srp
      return
      end
      subroutine rrfr (n,np,g,c,cp,d,dp)
*
* $Id: rrfr.f,v 2.1 94/09/29 18:43:51 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 one,zero,ak,c,ci,cp,ct,d,dd,di,dn,dp,dt,g
      parameter ( zero=0.0d0,one=1.0d0 )
      dimension g(0:n),c(0:np),cp(0:np),d(0:np),dp(0:np)
*
*     rrfr : converts polynomial g to a rational fraction
*              g    = vector g(k), k=0,n
*              c,cp = numerator polynomials, convergents n and n-1
*              d,dp = denominator polynomials, convergents n and n-1
*
      do 10  i = 1, np
         c(i) = zero
         cp(i) = zero
         d(i) = zero
         dp(i) = zero
   10 continue
      c(0) = g(0)
      cp(0) = g(0)
      d(0) = one
      dp(0) = one
      ct = zero
      dt = one
      dn = g(0)
      do 40 k = 1, n
         dd = dn
         dn = zero
         nq = ( k + 1 ) / 2
c$dir scalar
         do 20 i = 0, nq
            dn = dn + g(k-i) * d(i)
   20    continue
         ak = - dn / dd
         if ( ak .eq. zero ) return
         do 30 i = 1, nq
            ci = c(i)
            c(i) = ci + ak * ct
            ct = cp(i)
            cp(i) = ci
            di = d(i)
            d(i) = di + ak * dt
            dt = dp(i)
            dp(i) = di
   30    continue
         ct = g(0)
         dt = one
   40 continue
      return
      end
      subroutine rrfvl (nt,ip,xi,xc,xd,xcp,xdp,dr,dpr,aerr,aepr,ierr)
*
* $Id: rrfvl.f,v 2.1 94/09/29 18:43:53 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 xi,xc,xd,xcp,xdp,dr,dpr,aerr,aepr,csr,cspr,dsr,dspr,
     x       dr1,dpr1
      dimension xc(0:nt),xd(0:nt),xcp(0:nt),xdp(0:nt)
*
*     rrfvl : evaluates real rational fraction
*
      call rpolpp (xc,xi,csr,cspr,nt,ip)
      call rpolp (xd,xi,dsr,dspr,nt)
      dr = csr / dsr
      dpr = ( cspr * dsr - csr * dspr ) / ( dsr * dsr )
*     error estimates
      if ( ierr .ne. 0 ) then
         nt1 = nt - 1
         call rpolpp (xcp,xi,csr,cspr,nt1,ip)
         call rpolp (xdp,xi,dsr,dspr,nt1)
         dr1 = csr / dsr
         dpr1 = ( cspr * dsr - csr * dspr ) / ( dsr * dsr )
         aerr = abs(dr-dr1)
         aepr = abs(dpr-dpr1)
      endif
      return
      end
      subroutine lwint (nmax,mcore,pt0)
*
* $Id: lwint.f,v 2.1 94/09/29 18:44:00 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      common /io/ iin,iout
      common /lw/ vmat,enrg,tc,tev,tp,tsave,pdiag,r4,space,tpev,b1,
     x            rmtx,cupest,ipvt
*
      nm2 = nmax * nmax
*
      vmat = pt0
      enrg = vmat + (nm2+nmax)/2
      tc = enrg + nmax
      tev = tc + nm2
      tp = tev + nm2
      tsave = tp + nm2
      pdiag = tsave + nm2
      r4 = pdiag + 4 * nmax
      space = r4 + nm2
      tpev = space + 3 * nm2
      b1 = tpev + nm2
      rmtx = b1 + nm2
      cupest = rmtx + nm2
      ipvt = wpadti(cupest+nmax)
      last = iadtwp(ipvt+nmax) + 1
      if ( last .gt. mcore ) then
         write (iout,1020) last,mcore
         stop
      endif
      pt0 = last + 1
      return
*
 1020 format (' lwint : storage overflow, last =',i10,' mcore =',i10)
      end
      subroutine lwprop (enrg,tev,tpev,tsave,pdiag,r4,keep,nstp,
     x                  nmax,etot,prntr,debug,rstart,rfinal,
     x                  space,b1,rmtx,idiscr,ipvt,rmat)
*
* $Id: lwprop.f,v 2.1 94/09/29 18:44:02 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 enrg,tev,tpev,tsave,pdiag,r4,etot,rstart,
     x                 r,rfinal,space,zero,one,rmtx,dr,rl,b1,rr,
     x                 rmat
      logical prntr,debug
      parameter ( zero=0.0d0,one=1.0d0 )
      dimension enrg(nmax),tev(nmax,nmax),tsave(nmax,nmax),
     x          pdiag(nmax,2),r4(nmax,nmax),
     x          space(nmax,nmax,3),tpev(nmax,nmax),b1(nmax,nmax),
     x          rmtx(nmax,nmax),ipvt(nmax),rmat(*),keep(nstp)
*
*     lwprop    LW Propagator Driver
*
*     keep(nstep) decides how many channels are kept at each step.
      nstep = 1
      n = keep(nstep)
      k = 0
      do 20 j = 1, nmax
         do 10 i = 1, j
            k = k + 1
            r4(i,j) =  - rstart * rmat(k)
            r4(j,i) =  - rstart * rmat(k)
 10      continue
 20   continue
      if (prntr) call printm ('r4mt',r4,n,n,nmax)
      read (idiscr, rec = 1) r,dr,tev,tpev
      rr = r
      call dcopy (n*nmax,r4,1,rmtx,1)
      call dgemm ('n','n',n,n,n,one,rmtx,nmax,tev,nmax,zero,tsave,nmax)
      call dgemm ('n','n',n,n,n,one,tpev,nmax,tsave,nmax,zero,r4,nmax)
      do 50 nstep = 2, nstp
         n = keep(nstep)
         read (idiscr, rec = nstep) rl,dr,enrg,tsave,tev,tpev
         if (debug) call printm ('tev ',tev,n,n,nmax)
*     propagation code
         call rpropl (enrg,-dr,etot,r4,pdiag,tsave,n,nmax,space,ipvt)
         if (debug) call printm ('r4mt',r4,n,n,nmax)
         rr = rl + dr
         if ( rr .ge. rfinal ) go to 60
         rl = rr
   50 continue
   60 rfinal = rr
      if (prntr) call printm ('rfin',r4,n,n,nmax)
*     apply boundary conditions

      call dgemm ('n','n',n,n,n,one,r4,nmax,tpev,nmax,zero,b1,nmax)
      call dgemm ('n','n',n,n,n,one,tev,nmax,b1,nmax,zero,r4,nmax)

      k = 0
      do 140 j = 1, n
         do 130 i = 1, j
            k = k + 1
            rmat(k) = - r4(i,j) / rfinal
  130    continue
  140 continue
      return
      end

      subroutine mprint ( name, a, n, m, nmax )
*
* $Id: mprint.f,v 2.1 94/09/29 18:44:04 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 a
      character*4 name
      dimension a(nmax,m)
      common /io/ iin,iout
*
*     mprint -- prints a real matrix a(n,m)
*               name      -- 4 character identification
*               a(nmax,m) -- array to be printed
*               n         -- number of rows to be printed
*               m         -- number of columns to be printed ( .le. 10 )
*               nmax      -- row dimension of a
*
      write (iout,1000) name, (a(1,i), i = 1, m )
      if ( n .gt. 1 ) then
         do 10 i = 2, n
            write (iout,1010) (a(i,j), j = 1, m )
   10    continue
      endif
      return
*
 1000 format ( '0', a4, 2x, 10g12.4 )
 1010 format ( 7x, 10g12.4 )
      end


      subroutine pot (r,nchan,v,nv,nchx,ismax,ion,lchl,cf,ethr)
*
* $Id: pot.f,v 2.1 94/09/29 18:44:06 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 r,v,cf,ethr,zero,one,two,rk,elx,vc
      parameter ( zero=0.0d0,one=1.0d0,two=2.0d0 )
      dimension v(nv),lchl(nchx),ethr(nchx),cf(nchx,nchx,ismax)
*
*     potential at radius r
*
      do 10 i = 1, nv
         v(i) = zero
   10 continue
      vc = - two * dble(ion) / r
      do 50 k = 1, ismax
         rk = one / r**(k+1)
         ij = 0
         do 40 j = 1, nchan
            do 30 i = 1, j
               ij = ij + 1
               v(ij) = v(ij) + cf(i,j,k) * rk
   30       continue
   40    continue
   50 continue
      rk = one / ( r * r )
      do 60 i = 1, nchan
         elx = dble( lchl(i) * ( lchl(i) + 1 ) )
         ij = (i*(i+1))/2
         v(ij) = v(ij) + elx * rk + ethr(i) + vc
   60 continue
      return
      end



      subroutine potp (r,nchan,v,nv,nchx,ismax,ion,lchl,cf)
*
* $Id: potp.f,v 2.1 94/09/29 18:44:07 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 r,v,cf,zero,two,rk,elx,vc
      parameter ( zero=0.0d0,two=2.0d0 )
      dimension v(nv),lchl(nchx),cf(nchx,nchx,ismax)
*
*     potential gradient at radius r
*
      do 10 i = 1, nv
         v(i) = zero
   10 continue
      vc =  two * dble (ion) / ( r * r )
      do 50 k = 1, ismax
         rk = dble (k + 1) / r**(k+2)
         ij = 0
         do 40 j = 1, nchan
            do 30 i = 1, j
               ij = ij + 1
               v(ij) = v(ij) - cf(i,j,k) * rk
   30       continue
   40    continue
   50 continue
      rk = two / r**3
      do 60 i = 1, nchan
         ij = (i*(i+1))/2
         elx = dble( lchl(i) * ( lchl(i) + 1 ) )
         v(ij) = v(ij) - elx * rk + vc
   60 continue
      return
      end
      subroutine printm (name,array,n,m,nmax)
*
* $Id: printm.f,v 2.1 94/09/29 18:44:09 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 array
      character*4 name,namec
      dimension array(*)
*
*     printm -- driver for mprint
*
      namec(1:4) = name(1:4)
      nh = 1
      mp = 6
      i = 0
      iz = ichar('0')
      do 10 nl = m, 1, -mp
         nc = min( nl, mp )
         call mprint (namec,array(nh),n,nc,nmax)
         i = i + 1
         namec(4:4) = char(iz+mod(i,10))
         nh = nh + nmax * nc
   10 continue
      end
      subroutine printv (name,n,nv,vmat)
*
* $Id: printv.f,v 2.1 94/09/29 18:44:11 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 vmat
      character*4 name,namec
      parameter ( ncol=9 )
      dimension vmat(nv)
      common /io/ iin,iout
c
      namec(1:4) = name(1:4)
      iz = ichar('0')
      do 10 ncf = 1, n, ncol
         ncl = min(ncf+ncol-1,n)
         write (iout,1010) namec
         do 20 ir = ncf, n
            nc = (ir*(ir-1))/2
            nc1 = nc + ncf
            nc2 = nc + min(ir,ncl)
            write (iout,1000) (vmat(j),j=nc1,nc2)
   20    continue
         write (iout,*)
         namec(4:4) = char(iz+mod(mod(ncf,ncol)+1,10))
   10 continue
      return
 1000 format (9f8.4)
 1010 format (2x,a4)
      end
      subroutine rpropl (en,du,e,r4,dchan,tot,n,nmax,space,ipvt)
*
* $Id: rpropl.f,v 2.1 94/09/29 18:44:12 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 r4,dchan,tot,space,zero,one,en,du,e,ek,
     x                 hek,s,sh,ch
      parameter ( zero=0.0d0,one=1.0d0 )
      common /io/ iin,iout
      dimension r4(nmax*nmax),tot(nmax,n),dchan(nmax,2),en(nmax),
     x          space(nmax,nmax,3),ipvt(nmax)
*
*     space(1,1,1) = z inverse ; space(1,1,2) = sector r(2)
*     space(1,1,3) = sector r(3)
*
      do 10 i = 1, n
         ek = e - en(i)
         if ( ek .ge. zero ) then
            ek = sqrt(ek)
            hek = du * ek
            s = - one
            sh = sin(hek)
            ch = cos(hek)
         else
            s = one
            ek = sqrt(-ek)
            hek = du * ek
            sh = sinh(hek)
            ch = cosh(hek)
         endif
         dchan(i,1) = one / ( s * sh * ek )
         dchan(i,2) = ch * dchan(i,1)
   10 continue
*     build sector r(2) and start building sector r(1)
      n2 = n * nmax
      do 45 ij = 1, n2
         j = (ij-1)/nmax + 1
         i = ij - (j-1) * nmax
         space(i,j,1) = r4(ij)
         space(i,j,2) = tot(i,j) * dchan(j,1)
         space(i,j,3) = tot(i,j) * dchan(j,2)
 45   continue
*     finish sector r(1) and build z inverse
      call dgemm ('n','t',n,n,n,one,space(1,1,3),nmax,tot,nmax,
     x            one,space,nmax)
*     save sector r(3) and invert z inverse
      call transp (space(1,1,2),space(1,1,3),n,nmax)
      call dgetrf (n,n,space,nmax,ipvt,ier)
      if ( ier .ne. 0 ) then
         write (iout,1000) ier
         stop
      endif
      call dgetrs('n',n,n,space,nmax,ipvt,space(1,1,2),nmax,ier)
      if ( ier .ne. 0 ) then
         write (iout,1010) ier
         stop
      endif
*     z*r(2) computed in space(1,1,2) ; r(3)*zr(2) + old r4
      call dgemm ('n','n',n,n,n,-one,space(1,1,3),nmax,space(1,1,2),
     x            nmax,zero,r4,nmax)
      ii = 1
      n1 = nmax + 1
      do 150 i = 1, n
         r4(ii) = dchan(i,2) + r4(ii)
         ii = ii + n1
  150 continue
      return
 1000 format (' rpropl : error return from dgetrf, info = ',i6)
 1010 format (' rpropl : error return from dgetrs, info = ',i6)
      end
      subroutine slwpot (nstart,vmat,enrg,tc,tp,tsave,n,nmax,
     x                  prntv,debug,rstart,rfinal,lwsect,lwkeep,
     x                  lwran,beta,cuplim,space,tev,tpev,ismax,ion,
     x                  lchl,cf,ethr,idiscr,cupest)
*
* $Id: slwpot.f,v 2.1 94/09/29 18:44:14 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 vmat,enrg,tc,tp,tsave,rstart,cf,ethr,rlim,stp,
     x                 rfinal,drstrt,stpmin,stpmax,beta,space,tev,tpev,
     x                 expt,zero,one,half,dr,rl,step,r,e1,couple,vprime,
     x                 or,rr,oe1,cupest,cuplim,lwsect
      logical prntv,debug
      parameter ( zero=0.0d0,one=1.0d0,half=0.5d0,expt=1.0d0/6.0d0 )
      dimension vmat(nmax,nmax),enrg(nmax),tc(nmax,nmax),tev(nmax,nmax),
     x          tp(nmax,nmax),tsave(nmax,nmax),space(nmax,nmax,3),
     x          tpev(nmax,nmax),lchl(nmax),cf(nmax,nmax,ismax),
     x          ethr(nmax),cupest(nmax),rlim(2),stp(2),
     x          lwsect(lwran),lwkeep(lwran)
      common /io/ iin,iout
*
*     slwpot -  diagonalises the potential in each sector
*               nstart = 1 at first entry. On exit,
*               nstart = number of steps taken (nstart can be used to
*               restart slwpot ready to propagate further)
*
*     test whether sectors are required
*
      if (nstart .eq. lwran) then
         write(iout,1060) nstart
         go to 520
      end if
*
*     calculate the first step length using the gradient of the
*     potential at rstart and the max using rfinal
*
      n2 = n * n
      nv = ( n2 + n ) / 2
      ik = 1
      if (nstart .gt. 1) ik = 2
      rlim(1) = rstart
      rlim(2) = rfinal
      do 10 k = ik, 2
         call potp (rlim(k),n,vmat,nv,nmax,ismax,ion,lchl,cf)
         call vdiag (vmat,enrg,tc,n,space,nv,nmax)
         e1 = zero
         do 20 i = 1,n
            e1 = e1 + enrg(i)
   20    continue
         vprime = (e1/n) ** 2
         stp(k) = (beta / vprime) ** expt
   10 continue
      stpmax = stp(2)
      if (nstart .eq. 1) then
         lwkeep(1) = n
         stpmin = stp(1)
         drstrt = stpmin
         if (prntv .or. debug)
     x       write (iout,1000) rstart,rfinal,stpmin,stpmax,beta
         call pot (rstart,n,vmat,nv,nmax,ismax,ion,lchl,cf,ethr)
         if (prntv) then
            write (iout,1010) rstart
            call printv ('vcup',n,nv,vmat)
         endif
         call vdiag (vmat,enrg,tc,n,space,nv,nmax)
         do 35 ik = 1, n2
            i = (ik-1)/n + 1
            k = ik - (i-1) * n
            tev(i,k) = tc(i,k)
            tpev(k,i) = tc(i,k)
            tp(k,i) = tc(i,k)
   35    continue
         if (debug) write (iout,1020)
         write (idiscr,rec=nstart) rstart,drstrt,tev,tpev
         or = rstart
      else
         read (idiscr,rec = nstart+1) or,rstart,drstrt,(enrg(i),i=1,n),
     x                                ((tp(i,k),i=1,n),k=1,n)
      endif
      if (rfinal .le. rstart) then
         write(iout,1070) nstart,rstart,rfinal
         go to 520
      endif
      olim = lwkeep(nstart)
      lwsect(nstart) = rstart
      rl = rstart
      dr = drstrt
      stpmin = drstrt
      e1 = zero
      do 50 i = 1, n
         e1 = e1 + enrg(i)
   50 continue
      e1 = e1 / n
      do 500 nstep = nstart+1, lwran
         r = rl + half * dr
         rr = rl + dr
         lwsect(nstep) = rr
         call pot (r,n,vmat,nv,nmax,ismax,ion,lchl,cf,ethr)
         if (prntv) call printv ('vcup',n,nv,vmat)
         call vdiag (vmat, enrg, tc, n, space, nv, nmax)
         if (debug) call printm ('tcc ',tc,n,n,nmax)
         do 95 ixy = 1, n2
            ix = (ixy-1)/n + 1
            iy = ixy - (ix-1) * n
            tev(ix,iy) = tc(ix,iy)
            tpev(iy,ix) = tc(ix,iy)
   95    continue
         oe1 = e1
         e1 = zero
         do 110 i = 1, n
            e1 = e1 +  enrg(i)
  110    continue
         e1 = e1 / n
*     construct adiabatic overlap coupling matrix
         call dgemm ('n','n',n,n,n,one,tp,nmax,tc,nmax,zero,tsave,nmax)
         if (debug) call printm ('tsve', tsave,n,n,nmax)
         call dcopy (n*nmax,tpev,1,tp,1)
         write (idiscr,rec = nstep) rl, dr, enrg, tsave , tev , tpev
*
*     estimate the coupling between each channel with higher channels
         lim = 0
         j2 = 0
         do 100  j = 2, n
            j2 = j2 + 2
            cupest(j) = zero
            do 105  k = 1, j - 1
               cupest(j) = cupest(j) + tsave(j,k) * tsave(j,k)  +
     x                     tsave(k,j) * tsave(k,j)
  105       continue
            cupest(j) = sqrt(cupest(j)) / (dble(j2) * dr)
            if (cupest(j) .gt. cuplim) lim = j
  100    continue
         if (lim .lt. olim) then
            lim = olim - 1
         elseif (lim .gt. olim) then
            do 115 i = nstep-1, 1, -1
               lwkeep(i) = max (lwkeep(i),lim)
  115       continue
         endif
         lwkeep(nstep) = lim
         olim = lim
         nprnt = min(n,2)
         if (debug) then
            couple = zero
            do 120 i = 1, n
               couple = couple + abs( tsave(i,i) )
  120       continue
            couple = one - couple / dble(n)
            write (iout,1040) nstep, rl, r, rr, dr, couple,
     x                        (enrg(i),i=1,nprnt)
         endif
*     next step size code
         vprime = ( ( oe1 - e1 ) / ( or - r ) )**2
         step = ( beta / vprime )**expt
         step = max( step, stpmin )
*     to force propogation to end at rfinal, uncomment next statement
*        step = min( step, rfinal-rr )
         or = r
         rl = rr
         if (rr .ge. rfinal) then
            nstart = nstep
            go to 510
         endif
         step = min( step, stpmax )
         dr = step
  500 continue
      nstep = lwran
      nstart = lwran
  510 write (idiscr,rec = nstep+1) or,rl,dr,enrg,tpev
      if (prntv .or. debug) write (iout,1050) nstep,step
  520 return
*
 1000 format (' R-Matrix Propagation : Light-Walker Method',/,
     x        ' Initial radius =',f10.5,5x,'; Final radius   =',f10.5,/,
     x        ' Min. step size =',e15.5,   '; Max. step size =',e15.5,/,
     x        ' Beta =',e15.5)
 1010 format (' Potential at r =',f10.5)
 1020 format ( '1Light-Walker Propagation',/,' nstep ',
     x         ' r(left)  r(cntr)  r(rght)     step coupling',
     x         ' Lowest potl evalues:..' )
 1040 format ( i6, 3f9.3, f9.4, f9.4, 1p2e11.3 )
 1050 format (' nstep =',i5,5x,'step =',f15.6)
 1060 format (' subroutine slwpot',/,
     x        ' nstart =',i4,/,
     x        ' nstart = lwran - no more sectors available',/,
     x        ' increase parameter lwran')
 1070 format (' subroutine slwpot',/,
     x        ' nstart =',i4,' rstart =', f9.4,/,
     x        ' rfinal =', f9.4,'no sectors required')
      end
      subroutine transp ( a, b, n, nmax )
*
* $Id: transp.f,v 2.1 94/09/29 18:44:18 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8   a, b
      dimension a(nmax,n), b(nmax,n)
*
*     transp -- copies the transpose of a(nmax,n) to b(nmax,n)
*
      n2 = n * n
      do 10 ij = 1, n2
         i = (ij-1)/n + 1
         j = ij - (i-1)*n
         b(j,i) = a(i,j)
   10 continue
      return
      end
      subroutine vdiag (v,e,t,n,work,nv,nmax)
*
* $Id: vdiag.f,v 2.1 94/09/29 18:44:20 vmb Exp Locker: vmb $
*
      integer n,nmax,info,iin,iout
      real*8 v,e,t,work
      common /io/ iin,iout
      dimension e(n), t(nmax,n), work(3*n), v(nv)
*
*     v(nv)     real symmetric matrix to be diagonalized
*     e(n)      eigenvalues of v
*     t(nmax,n) eigenvectors of v
*
      call DSPEV ( 'V', 'U', N, v, e, t, nmax, WORK, INFO )
      if ( info .ne. 0 ) then
         write (iout,1000) info
         stop
      endif
      return
*
 1000 format (' vdiag : error return by dspev, info =', i3)
      end
      subroutine asyms (ie,e,rmat,akmat,nchan,nschan,ispch,ichsp,
     x                  lchl,ion,lamax,cf,r,bbloch,
     x                  ethr,f,fp,et,nbasis,nampx,nranmx,rpflg,
     x                  prtflg,icol,apeps,degeny,iasy,lmax,
     x                  ewron,nopen,raf,raflw,asect,lwsect,lwkeep,
     x                  keep,nfinal,lwran,iapflg,diaflg,idiscr,
     x                  lwdb,beta,cuplim,deplim,ncrit,rfact,
     x                  gfact,maxtry,radius,rstats,rafix,x,ix)
*
* $Id: asyms.f,v 2.1 94/09/29 18:39:17 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 e,rmat,akmat,cf,bbloch,r,x,f,fp,et,ewron,apeps,degeny,
     x     rafin,asect,raflw,raf,lwsect,rfact,gfact,rafmax,ethr,
     x     r1,r2,beta,gmax,gmaxp,radius,rafix,rlwmax,cuplim,deplim
      logical fail,lwdb
      dimension lchl(nchan),cf(nchan,nchan,*),ix(*),
     x          nbasis(nranmx),rmat(*),akmat(*),x(*),
     x          prtflg(*),iapflg(4),rpflg(6),diaflg(20),
     x          f(nchan*nchan*2),fp(nchan*nchan*2),et(nchan),asect(*),
     x          raf(*),lwsect(*),radius(maxtry),ethr(nchan),
     x          rstats(maxtry),nschan(2),ispch(nchan),ichsp(nchan),
     x          lwkeep(lwran),keep(lwran),lwdb(2)
      common /io/ iin,iout
      common /pptr/ ampa,eigen,bigvec,wo
      common /kp/ bb,iwrk,bigvc
      common /lw/ vmat,enrg,tc,tev,tp,tsave,pdiag,r4,space,tpev,b1,
     x            rmtx,cupest,ipvt
      save gmaxp,nopenp
*
*   keep channels if less deeply bound than deplim
*   but keep at least ncrit channels
*   further keep at each propagation sector at least the number of
*   channels decided by the strength of the coupling in the LW
*   package
*
      do 40  i = nopen + 1,nchan
         if (et(i) .lt. -deplim) then
            nkeep = i - 1
            go to 45
         end if
 40   continue
      nkeep = nchan
 45   nkeep = min (nkeep + ncrit, nchan)
      do 50  i = 1,nfinal
         keep(i) = max(lwkeep(i), nkeep)
 50   continue
*
*  calculate gailitis coefficients for all channels as these are
*  independent of radius
*
      call gailco (nchan,iapflg,et,lchl,ion,cf,lamax,iasy,
     x             degeny,gmax,x,ix)
*
*  on the basis of the magnitude of the maximum coefficient decide
*  whether to use at least the final radius calculated for the previous
*  electron energy so long as a threshold has not been crossed
*
      if ( ie .gt. 1  ) then
         if ( nopen .eq. nopenp  .and. gmax .ge. gmaxp * gfact ) then
            raf(ie) = max (raf(ie),raf(ie - 1))
         end if
      end if
      rafin = raf(ie)
      rafmax = lwsect(nfinal)
      rlwmax = rafix
      if (nfinal .eq. lwran) rlwmax = lwsect(lwran)
*
*  redefine the trial radius to lie on a sector end
*
      call secend (rafin,nrange,nranmx,asect,isect,nfinal,lwsect,lwran)
*
*  try at successive radii until obtain converged wavefunctions
*
      do 60 itry = 1,maxtry
         ntry = itry
         raf(ie) = rafin
         nc = keep(isect)
*
*  calculate wavefunctions at rafin for the retained channels
*
         call gailit (f,fp,nc,nchan,iapflg,rafin,et,lchl,ion,
     x                iasy,apeps,lmax,x)
         radius(itry) = rafin
*
*  test the multi-channel wronskian for convergence
*
         call wronsk (nc,nopen,f,fp,diaflg(5),diaflg(6),ewron,fail)
         if ( diaflg(4) .eq. 1 ) then
            k1 = 1
            do 70 k = 1, 2
               write (iout,1110) rafin
               call wrtmat (f(k1),nc,nc,nc,icol)
               k1 = k1 + nc*nc
   70       continue
            k1 = 1
            do 80 k = 1, 2
               write (iout,1120) rafin
               call wrtmat (fp(k1),nc,nc,nc,icol)
               k1 = k1 + nc*nc
   80       continue
         end if
*
*  if the test fails and rafin is not at its allowed maximum try
*  extending rafin by a factor rfact
*
         if ( fail .and. rafin .lt. rlwmax ) then
            gmaxp = gmax
            rafin = min(rfact * rafin, rlwmax)
*
*  if necessary extend the range over which the propagator is
*  initiated otherwise redefine the new trial radius at a sector
*  end
*
            if (rafin .gt. rafmax) then
               nrange = nranmx
               r1 = rafmax
               r2 = rafin
               nf = nfinal
               call slwpot (nfinal,x(vmat),x(enrg),x(tc),x(tp),x(tsave),
     x                      nchan,nchan,lwdb(1),lwdb(2),r1,r2,lwsect,
     x                      lwkeep,
     x                      lwran,beta,cuplim,x(space),x(tev),x(tpev),
     x                      lamax,ion,lchl,cf,ethr,idiscr,x(cupest))
               rafin = lwsect(nfinal)
               if (nfinal .eq. lwran) rlwmax = lwsect(lwran)
               do 100 i = nf + 1,nfinal
                  keep(i) = max(lwkeep(i),nkeep)
 100           continue
               isect = nfinal
               rafmax = rafin
            else
               call secend (rafin,nrange,nranmx,asect,isect,nfinal,
     x                      lwsect,lwran)
            end if
         else
            go to 90
         end if
 60   continue
      write (iout,1130) ntry, rafin
      stop
 90   if (fail) then
         if (nfinal .lt. lwran) write (iout,1140) rafin
         if (nfinal .eq. lwran) write (iout,1145) rafin, lwran
      end if
      nopenp = nopen
      rstats(ntry) = rstats(ntry) + 1
      if (ie .eq. 1) gmaxp = gmax
      if (diaflg(12) .eq. 1) then
         write (iout,1150)ie,e,(radius(i),i=1,ntry)
      end if
*
*     propagate r-matrix to radius rafin
*     use BBM propagator as far as raflw and LW to rafin
*
      if (raflw .ne. r ) then
         if ( diaflg(3) .ne. 0 ) then
            write (iout,*)'rmat before splitting spins'
            call wrttmt (rmat,nchan,icol)
         end if
         call splitr(nchan,rmat,x(bigvec),ichsp)
         if ( diaflg(3) .ne. 0 ) then
            write (iout,*)'rmat after splitting spins'
            call wrttmt (rmat,nchan,icol)
         end if
         call rprop2 (rmat,e,nschan,nbasis,x(ampa),x(eigen),x(bigvec),
     x                nranmx,nampx,rpflg(4),1,nrange,asect,0)
         call backr (nchan,rmat,x(bigvec),ispch)
         if ( diaflg(3) .ne. 0 ) then
            r1 = min (raflw,rafin)
            write (iout,1160) r1
            call wrttmt (rmat,nchan,icol)
         endif
      end if
      if (raflw .lt. rafin) then
         call lwprop (x(enrg),x(tev),x(tpev),x(tsave),x(pdiag),x(r4),
     x                keep,nfinal,nchan,e,lwdb(1),lwdb(2),raflw,
     x                rafin,x(space),x(b1),x(rmtx),idiscr,ix(ipvt),
     x                rmat)
         if ( diaflg(3) .ne. 0 ) then
            write (iout,1170) rafin
            call wrttmt (rmat,keep(nfinal),icol)
         endif
      end if
      call kmat (nc,rafin,bbloch,x(bigvc),x(bb),ix(iwrk),nopen,
     x           f,fp,rmat,akmat)
      if ( prtflg(1) .gt. 1 ) then
         write (iout,1180)
         call wrttmt (akmat,nopen,icol)
      endif
      return
*
 1110 format (/' solutions at rafin=',f9.3)
 1120 format (/' derivatives at rafin',f9.3)
 1130 format (/' failed to calculate accurate wave functions'/
     x        /' tried ',i3,' times.  final radius was ',f9.3)
 1140 format (/,' *****  warning - wronskian test failed  at ',f9.3,/
     x          ' the calculation continues at this radius',/,
     x          ' this radius is limited by rafix set in subr cntrl',
     x          ' *********' )
 1145 format (/,' *****  warning - wronskian test failed  at ',f9.3,/
     x          ' the calculation continues at this radius',/,
     x          ' this radius is at the end of the last possible',
     x          ' LW sector,lwran =',i4,/,
     x          ' lwran is set in subr cntrl',' *********' )
 1150 format (/12x,'energy    radii at which functions calculated',/
     x        ' ',i5,2x,f12.6,5(2x,f12.5),/,23x,5(f12.5,2x),
     x        /,23x,5(f12.5,2x),/,23x,5(f12.5,2x))
 1160 format (/,' r-matrix after propagation by BBM at radius =',f10.5)
 1170 format (/,' r-matrix after propagation by LW at radius =',f10.5,/)
 1180 format (/' k-matrix :')
      end
      subroutine backr(nchan,rmat,sprmat,ichsp)
*
* $Id: backr.f,v 2.1 94/09/29 18:40:25 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 rmat,sprmat
      dimension rmat(*),sprmat(nchan,nchan),ichsp(nchan)
*
*     split rmat into spins by rearranging channels
*
      kk = 0
      do 10 i2 = 1,nchan
         ich2 = ichsp(i2)
         do 10 i1 = 1,i2
            kk = kk + 1
            ich1 = ichsp(i1)
            sprmat(ich1,ich2) = rmat(kk)
            sprmat(ich2,ich1) = rmat(kk)
 10   continue
      kk = 0
      do 20 i2 = 1,nchan
         do 20 i1 = 1,i2
            kk = kk + 1
            rmat(kk) = sprmat(i1,i2)
 20   continue
      return
      end
      double precision function but0 (nbut,fkn,uk,u)
*
* $Id: but0.f,v 2.1 94/09/29 18:39:20 vmb Exp Locker: vmb $
*
      implicit real*8 (a-h,o-z)
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0,pt4=0.4d0,
     x           sr=0.04d0,srm=-0.04d0,closep=0.3d0)
      logical pole
      dimension fkn(0:nbut),uk(0:nbut)
*
*     buttle correction
*
      s1 = one / dfloat(3)
      s2 = one / dfloat(45)
      s3 = two / dfloat(945)
      s4 = one / dfloat(4725)
      b = zero
*
*     u .gt. sr
*
      if ( u .gt. sr ) then
         fk = sqrt(u)
         pole = .false.
         do 10 n = 0, nbut
            if ( abs(fk-fkn(n)) .gt. closep ) then
               b = b + one / ( u - uk(n) )
            else
               pole = .true.
               d1 = fk - fkn(n)
            endif
   10    continue
         if (pole) then
            d2 = d1**2
            d = d1 * ( s1 + d2 * ( s2 + d2 * ( s3 + s4 * d2 ) ) )
            b = two * b + ( d + one / ( two * fk - d1 ) ) / fk
         else
            b = two * b + tan(fk) / fk
         endif
*
*     sum for u .lt. sr
*
      else
         do 20 n = 0, nbut
            b = b + one / ( u - uk(n) )
   20    continue
*
*    u .lt. sr and u .gt. srm=-sr
*
         if ( u .gt.  srm ) then
            b = two * b + one + s1 * u * ( one + pt4 * u )
*
*     u lt. -sr
*
         else
            fk = sqrt(-u)
            b = two * b + tanh(fk) / fk
         endif
      endif
      but0 = b
      return
      end
      subroutine butco (cf,cfbut1,lchl,nchan,lrang2)
*
* $Id: butco.f,v 2.1 94/09/29 18:40:08 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8  cf,cfbut1
*     dimension cf(3,lrang2),cfbut1(3,lrang2),lchl(nchan)
      dimension cf(3,nchan),cfbut1(3,lrang2),lchl(nchan)
*
*     butco selects the required buttle coefficients
*
      do 10 i= 1,nchan
         l = lchl(i) + 1
         do 20 j = 1,3
*DEBUG
*        print *,'j = ',j,'    i = ',i,'      l = ',l
            cf(j,i) = cfbut1(j,l)
 20      continue
 10   continue
      return
      end
      subroutine case(llo,lup,ldel,nlran,lrglx,ncases,mintsp,
     x                maxtsp,minlsp,maxlsp,lrgs,lrgl,lrgpi)
*
* $Id: case.f,v 2.1 94/09/29 18:40:55 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      common /io/ iin, iout
      dimension llo(nlran),lup(nlran),ldel(nlran),lrgs(ncases),
     x          lrgl(ncases),lrgpi(ncases)
*
*   construct the list of slpi cases
*
      icase = 0
      do 30 i = 1,nlran
         ld = max(ldel(i),1)
         do 40 l = llo(i), lup(i), ld
            if (l .le. lrglx) then
               do 50 ipi = 0,1
                  do 60 ilsp =  minlsp, maxlsp, 2
                     icase = icase + 1
                     lrgl(icase) = l
                     lrgpi(icase) = ipi
                     lrgs(icase) = ilsp
 60               continue
 50            continue
            else
               do 70 ipi = 0,1
                  do 80 itsp = mintsp, maxtsp, 2
                     icase = icase + 1
                     lrgl(icase) = l
                     lrgpi(icase) = ipi
                     lrgs(icase) = -itsp
 80               continue
 70            continue
            end if
 40      continue
 30   continue
*
      if (icase .gt. ncases ) then
          write(iout,1000)icase, ncases
          stop
      end if
*
      return
*
 1000 format(' **** error in case - array overflow',/
     x        'icase, ncases',2i4)
      end



      subroutine cbut (fkn,uk,nbutx)
*
* $Id: cbut.f,v 2.1 94/09/29 18:40:07 vmb Exp Locker: vmb $
*
      implicit integer  (a-z)
      real*8 fkn,uk,g,pi,half
      parameter  (half = 0.5d0,pi = 3.14159265d0)
      dimension fkn(0:nbutx),uk(0:nbutx)
*
*     cbut sets up the arrays needed by but0
*
      g = half * pi
      g = -g
      do 10 i = 0,nbutx
         g = g + pi
         fkn(i) = g
         uk(i) = g * g
 10   continue
      return
      end

      subroutine cntrl(x,ix,mcore)
*
* $Id: cntrl.f,v 2.1 94/09/29 18:39:14 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      character*80 title
      character*128 fout,filh,fkmat,fpha,ftmat,fcros,fom,ftop,
     x              fsmat,flwp,frmp
      real*8 zero,x,dfloat,bbloch,rmatr, apeps, rafin, sw1,
     x       quad,esc1,abvthr,belthr,delta,
     x       ewron,rfirst,rlast,degeny, e0,
     x       ebig, scalx, rscalx, flmin,dlamch,one,
     x       facmin,rafest,rafmax,raflw,raflw1,beta,
     x       lwsect, flion, cuplim, cuplm, mingap, ming,
     x       deplim, rfact,gfact,rafix,radwm,evhigh
      logical jj,lwdb,finis,adwm
      parameter ( one=1.0d0,facmin=1.0d+15 )
      parameter ( nbmx=100,lwran=1000,ncase=80,nlran=10 )
      dimension nschan(2),nspin(2),idscbb(2),
     x          ne1(10),esc1(20),ls(ncase),ll(ncase),
     x          lpi(ncase),llo(nlran),lup(nlran),ldel(nlran),
     x          x(mcore),ix(*),
     x          nbasis(nbmx),rpflg(6),iapflg(4),
     x          diaflg(20),prtflg(7),nmt(7),lwsect(lwran),
     x          lwkeep(lwran),keep(lwran),lwdb(2),frmp(2)
      common /io/ iin, iout
      common /lw/ vmat,enrg,tc,tev,tp,tsave,pdiag,r4,space,tpev,b1,
     x            rmtx,cupest,ipvt
      common /mach/ scalx,rscalx
      namelist /phzfil/fout,filh,fkmat,fpha,ftmat,fcros,fom,ftop,
     x                 fsmat,flwp,frmp,
     x                 iout,nfth,nmt,idiscr,idscbb
      namelist /phzin/ title,prtflg,diaflg,ne1,esc1,
     x                 ls,ll,lpi,jj,llo,lup,ldel,lrglx,lrglr,
     x                 abvthr,belthr,degeny,ewron,apeps,iasy,rafin,
     x                 rafix,rfact,gfact,maxtry,iapflg,
     x                 idbbw,rpflg,raflw,beta,cuplim,deplim,ncrit,
     x                 lwdb,icon,quad,ncol,ibuttl,lamax,radwm,
     x                 evhigh,nphys
*
*  phzin default values
*
      data prtflg/0,0,0,3,0,0,0/,diaflg/20*0/,ne1/10*0/,
     x     ll/80*-1/,jj/.false./,llo/10*0/,lup/10*-1/,ldel/10*1/,
     x     lrglx/12/,lrglr/1000/,
     x     abvthr/1.0d-3/,belthr/1.0d-3/,degeny/1.0d-3/,ewron/1.0d-4/,
     x     apeps/1.0d-4/,iasy/20/,
     x     rafix/1.0d4/,rfact/1.333d0/,gfact/0.5d0/,maxtry/33/,
     x     iapflg/4*0/,idbbw/0/,rpflg/6*0/,raflw/0.0d0/,
     x     beta/1.0d-6/,cuplim/1.0d-6/,deplim/1.0d-2/,ncrit/1/,
     x     lwdb/.false.,.false./,icon/1/,quad/0/,ncol/4/,
     x     lamax/20/,ibuttl/1/,radwm/2.0d0/,evhigh/1.0d30/,nphys/0/
*
*  phzfil default values
*
      data fout/'OUT'/,filh/'H'/,
     x     fkmat/'KMATB'/,fpha/'PHASEF'/,ftmat/'TMATB'/,
     x     fcros/'CROSSB'/,fom/'OMEGAB'/,ftop/'TOPB'/,
     x     fsmat/'S00'/,
     x     flwp/'lwpscr'/,frmp/'bbscr1','bbscr2'/
      data nfth/10/,idscbb/11,12/,idiscr/13/,
     x     nmt/21,22,23,24,25,26,30/
*
*  topin default values for NO top up, prtflg(6)=0 (see subroutine topdat)
*
      data lrgle/999/,enum2/0/,nl/0/
*
      data zero/0.0d0/
*
*   data criteria
*   abvthr, belthr > 10**-8, sqrt(accur) in coulfg
*   rafix, max value of rafin, < abort/sqrt(max e), in coulfg
*   maxtry defined by rmatr*rfact**maxtry .ge. rafix
*
      iout = 3
      read (iin,phzfil,end=18,err=18)
      go to 19
 18   write (6,phzfil)
      stop
 19   open (iout,file = fout,
     x      access = 'sequential',status = 'unknown',
     x      form = 'formatted')
      write (iout, 1000)
*
*  calculate the safe minimum for the machine
*  calculate the scaling factor for the gailitis coefficients allowing
*  the factor facmin above the minimum to prevent underflow
*  calculate the reciprocal
*
c      flmin = dlamch('S')
cz - ?
      flmin = 1.0d-95
      scalx = facmin * flmin
      rscalx = one / scalx
      write(iout,*)  'scalx,rscalx', scalx,rscalx
*
      do 10 i = 1,80
          title(i:i) = ' '
   10 continue
*
*  read first record on H file
*
      call rdhed(nfth,filh,nelc,nz,lrang2,ismax,ntarg,rmatr,
     x           bbloch)
      write (iout,2010) filh,nelc,nz,lrang2,ntarg,ismax,rmatr
      write (iout,2030)
      ion = nz - nelc
      zion = ion
*
*  for positive ions use z-scaling
*
      if (ion .gt. 0) then
         ion = 1
         flion = dfloat (zion)
         rmatr = flion * rmatr
      end if
*
*  set storage for target arrays
*
      etarg = 1
      ltarg = wpadti (etarg + ntarg+1)
      starg = ltarg + ntarg
      cfbut1 = iadtwp (starg + ntarg + 1)
      tlast = cfbut1 + 3 * lrang2
      itlast = wpadti(tlast + 1)
      if (tlast .gt. mcore ) then
         write (iout, 1010) tlast,mcore
*        stop
      end if
*     print *,' after target arrays:',tlast,mcore
*
      it = itlast
      mcor1 = mcore -tlast +1
      call rdtar (nfth,x(etarg),ix(ltarg),ix(starg),x(cfbut1),ntarg,
     x            lrang2,zion,nbutx,maxcas,itlast,ix(it),mcor1)
      if ( zion .gt. 0) then
         nxsn = (ntarg * ( ntarg - 1 )) / 2
      else
         nxsn = (ntarg * ( ntarg + 1 )) / 2
      end if
      iind = itlast
      jind = iind + nxsn
      itlast = jind + nxsn
      tlast = iadtwp (itlast + 1)
*     print *,' before xsnind:',tlast,mcore
      call xsnind(ntarg,nxsn,zion,ix(iind),ix(jind))
*
*  set storage for line strength and initialise
*
      sline = tlast
      tlast = sline + nxsn
      itlast = wpadti (tlast + 1)
*     print *,' after line strength:',tlast,itlast,mcore
      do 17 i = 1,nxsn
         x(sline + i - 1) = -one
 17   continue
*
*  set storage for buttle correction
*
      if (ibuttl .ne. 0 .and. nbutx .ne. 0) then
         fkn = tlast
         uk = fkn +  nbutx + 1
         tlast = uk + nbutx + 1
         itlast = wpadti(tlast + 1)
*        print *,'after buttle stor:',tlast,itlast,mcore
         if (tlast .gt. mcore ) then
            write (iout, 1010) tlast,mcore
*           stop
         end if
         call cbut (x(fkn),x(uk),nbutx)
      end if
      rafin = rmatr
      read (iin,phzin,end=20,err=20)
C
C  Change that sets nphys to ntarg if the default value of 0 has been used
C
      if (nphys.eq.0) nphys = ntarg
C
C  end of Change
C
      go to 21
 20   write (iout,phzin)
      stop
*
*  final radius reset to the rmatrix radius if it had been wrongly input
*  to a value less than rmatr.   It is reset to rafix if it has been input
*  at a value too large for the coulomb wave function package.
*  raflw is the radius after which the LW propagator is used.   This must
*  also be initialised to at least the rmatrix radius and less than rafix.
*
 21    rafin = max(rafin,rmatr)
       rafin = min(rafin,rafix)
       raflw = max(raflw,rmatr)
       raflw = min(raflw,rafix)
*
*  allow the number of multipoles retained in the potential to be reduced
*  by the phzin namelist value of lamax.  If lamax = 0  the channel
*  wavefunctions will be coulomb functions. If the input value of rafin
*  causes propagation the long range potential is included in the propagation.
*  If rafin is defaulted there will be no propagation beyond rmatr and the
*  long range potential will be neglected entirely.
*
      lamax = min(lamax,ismax)
      if (lamax .eq. 0) then
         iasy = 2
      else if (lamax .lt. ismax) then
         ismax = lamax
      end if
      if (ismax .eq. 0) then
         rafin = rmatr
         raflw = rmatr
      end if
*
*  ensure iasy is an even number gt 0
*
       iasy = (iasy/2) * 2
       iasy = max(2,iasy)
*
*  print basic input data or complete namelist phzin if diaflg(1) set
*
      write(iout,1020)title
      if ( diaflg(1) .gt. 0 ) then
         write (iout,phzin)
      else
         call prtdat(prtflg,rafin,raflw,jj,lrglx,lrglr,radwm,
     x               lamax,iasy,abvthr,belthr,degeny,ewron,
     x               beta,cuplim,deplim,ncrit)
      end if
*
*  check prtflg
*
      if (prtflg(4) .gt. 0  .and.  prtflg(5) .gt. 0) then
          write(iout,1150)
          stop
      end if
      if (prtflg(6) .gt. 0  .and.
     x    prtflg(4) .eq. 0  .and.  prtflg(5) .eq. 0) then
         write(iout,1160)
         stop
      end if
      if (prtflg(6) .gt. 0  .and. prtflg(7) .ne. 0 ) then
         write(iout,1165)
         stop
      end if
*
*  open files required for output
*
      if ( prtflg(1) .eq. 1  .or.  prtflg(1) .eq. 3) then
         open (nmt(1),file = fkmat,access = 'sequential',
     x         status = 'unknown',form = 'unformatted')
         write(iout,2011) fkmat
      end if
      if ( prtflg(2) .eq. 1  .or.  prtflg(2) .eq. 3) then
         open (nmt(2),file = fpha,access = 'sequential',
     x         status = 'unknown',form = 'formatted')
         write(iout,2011) fpha
      end if
      if ( prtflg(3) .eq. 1  .or.  prtflg(3) .eq. 3) then
         open (nmt(3),file = ftmat,access = 'sequential',
     x         status = 'unknown',form = 'unformatted')
         write(iout,2011) ftmat
      end if
      if ( prtflg(4) .eq. 1  .or.  prtflg(4) .eq. 3) then
         open (nmt(4),file = fcros,access = 'sequential',
     x         status = 'unknown',form = 'unformatted')
         write(iout,2011) fcros
      end if
      if ( prtflg(5) .eq. 1  .or.  prtflg(5) .eq. 3) then
         open (nmt(5),file = fom,access = 'sequential',
     x         status = 'unknown',form = 'unformatted')
         write(iout,2011) fom
      end if
      if ( prtflg(6) .eq. 1  .or.  prtflg(6) .eq. 3) then
         open (nmt(6),file = ftop,access = 'sequential',
     x         status = 'unknown',form = 'unformatted')
         write(iout,2011) ftop
      end if
      if ( prtflg(7) .eq. 1  .or.  prtflg(7) .eq. 3) then
         open (nmt(7),file = fsmat,access = 'sequential',
     x         status = 'unknown',form = 'unformatted')
         write(iout,2011) fsmat
      end if
*
*  direct access files are opened for smatrices in prsmat.
*
*  make thresholds closer than degeny truly degenerate
*  print the new energy levels if changed anf write headings to
*  output files for k-matrix and t-matrix
*
      call degtar (nelc,nz,ntarg,x(etarg),ix(ltarg),ix(starg),
     x             degeny,prtflg,nmt,title)
*
*  read ranges of L or J
*
      call rdlran (jj,llo,lup,ldel,nlran,lrglx,ncases,ntarg,
     x             ix(starg),mintsp,maxtsp,minlsp,maxlsp)
*
*  if any ranges of L have been specified construct pointers for
*  arrays of S L Pi
*  else copy arrays of S L Pi from namelist
*
      if (ncases .gt. 0) then
         lrgs = itlast
         lrgl = lrgs + ncases
         lrgpi = lrgl + ncases
         itlast = lrgpi + ncases
         tlast = iadtwp (itlast + 1)
*        print *,' in ncases block:',tlast,itlast,mcore
         if (tlast .gt. mcore ) then
            write (iout, 1010) tlast,mcore
*           stop
         end if
         call case(llo,lup,ldel,nlran,lrglx,ncases,mintsp,maxtsp,
     x             minlsp,maxlsp,ix(lrgs),ix(lrgl),ix(lrgpi))
      else
         lrgs = itlast
         lrgl = lrgs + ncase
         lrgpi = lrgl + ncase
         itlast = lrgpi + ncase
         tlast = iadtwp (itlast + 1)
         call rdcas (ls,ll,lpi,ix(lrgs),ix(lrgl),ix(lrgpi),ncase,
     x               ncases)
      end if
*
*  set storage for energy mesh
*
      nume1 = necnt(ne1)
      sce1 = tlast
      noch = wpadti (sce1 + nume1 + 1)
      itlast = noch + nume1
      tlast = iadtwp ( itlast + 1)
*     print *,' after energy mesh :',tlast,itlast,mcore
      scest = tlast
      norder = wpadti(scest + nume1 + 1)
*
*  store the energy mesh after pruning out energies too close to
*  threshold
*
      call evals (ntarg,nume1,ne1,esc1,x(etarg),abvthr,belthr,x(sce1),
     x            enum1,ebig,x(scest),ix(norder))
      if (enum1 .lt. nume1) then
         write (iout,2020) enum1
      end if
*
*  read in data and set up arrays required for outputting a file for the
*  top up program
*
      if (prtflg(6) .gt. 0) then
         call topdat(ntarg,x(etarg),abvthr,belthr,ebig,sce2,
     x               enum2,ncases,ix(lrgl),lrgle,nl,fl,lstart,lfin,
     x               nxsn,xsnl,tlast,x,ix)
         if (lstart .ge. lrgle) enum1 = 0
      end if
*
*  set storage for accumulation of cross sections or collision strengths
*
      if (prtflg(4) .ne. 0  .or.  prtflg(5) .ne. 0) then
         xslpi = tlast
         xsn = xslpi + nxsn
         len = nxsn * enum1
         tlast = xsn + len
*        print *,' after cross section accum.',tlast,itlast,mcore
         do 30 i = 1,len
            x(xsn + i -1) = zero
 30      continue
      end if
*
*  set storage for s-matrices for line broadening
*  and initialise
*
      if (prtflg(7) .ne. 0) then
         itlast = wpadti(tlast + 1)
         iet = itlast
         call prsmat(enum1,x(sce1),ntarg,x(etarg),ix(ltarg),
     x               ix(iet),nmt(7),ncases,sline,nelc,nz,itlast,
     x               x,ix)
         tlast = iadtwp(itlast +1)
         if (tlast .gt. mcore ) then
            write (iout, 1010) tlast,mcore
*           stop
         end if
      end if
      rafest = rafin
      raflw1 = raflw
*
*  set store for arrays concerned with successive trial radii
*
      radius = tlast
      rstats = wpadti(radius + maxtry + 1)
      itlast = rstats + maxtry
      tlast = iadtwp(itlast + 1)
*print tlast
*     print *,'tlast = ',tlast
*
*  loop over SLPI cases
*
*  if multipole coefficients have been read decide to calculate
*  line strengths if LS coupling is used
*
      finis = jj
      if (ismax .eq. 0) finis = .true.
*
      do 50 icase = 1,ncases
         call rdhslp(nfth,ncases,icase,ix(lrgl),ix(lrgs),ix(lrgpi),
     x               lrgl2,nspn2,npty2,nchan,nstat,iflg,ix,mcore)
         if (iflg .lt. 0) go to 50
*
         write (iout, 2000) nspn2,lrgl2,npty2
         if (lrgl2 .lt. lrgle) then
            enum = enum1
            sce = sce1
         else
            enum = enum2
            sce = sce2
         end if
*  if asymptotic distorted wave method is to be used set starting
*  radius as radwm else use the radius from the H-file
*
         if (lrgl2 .gt. lrglr) then
            rfirst = radwm
            adwm = .true.
         else
            rfirst = rmatr
            adwm = .false.
         end if
*
*  set store ready to read H file
*
         lwmat = nstat * nchan
         leig = nstat
         eig = tlast
         wmat = eig + leig
         lchl = wpadti (wmat + lwmat + 1)
         ichl = lchl + nchan
         nltarg = ichl + nchan
         ncn2 = nchan * nchan
         lcf = ismax * ncn2
         cf = iadtwp (nltarg + ntarg + 1)
         cc = cf + lcf
         last = cc + nchan
         ilast = wpadti (last + 1)
         if (last .gt. mcore) then
            write (iout, 1010) last,mcore
*           stop
         end if
*
*  read channel information for current SLPI case and print if
*  diaflg(9) is set
*
         mcor1 = mcore - last + 1
         call rdch (nfth,ix(nltarg),ix(lchl),ix(ichl),x(cf),x(eig),
     x              x(wmat),ntarg,nstat,nchan,ncn2,ismax,rmatr,zion,
     x              diaflg(9),x,ix,mcor1,evhigh,nphys)
*
*  calculate line strengths
*
      if (.not. finis) then
         call lstren(ntarg,ix(nltarg),ix(ltarg),ix(starg),nxsn,
     x               x(sline),x(cf),nchan,ismax,ix(lchl),lrgl2,finis)
      end if
*
*  set storage for case dependent arrays
*
         raf = last
         ethr = raf + enum
         et = ethr + nchan
         f = et + nchan
         fp = f + 2 * ncn2
         rmat = fp + 2 * ncn2
         last = rmat + (ncn2 + nchan)/2 + 1
*
*  store relevant buttle coefficients for this case or omit for
*  debug purposes
*
         write (iout,*) ibuttl,ibuttl,ibuttl
         if (ibuttl .ne. 0) then
            cfbut = last
            last = cfbut + 3 * nchan + 1
*           write (iout,1763)
1763        format(' before butco',/)
            call butco (x(cfbut),x(cfbut1),ix(lchl),nchan,lrang2)
*           write (iout,1764)
1764        format(' after butco',/)
         end if
*
*   calculate channel parameters
*   ethr = threshold energy for each channel
*   cc   = l(l+1) for each channel
*   noch = number of open channels at each energy
*
         e0 = x(etarg)
*        write (iout,1765)
1765     format(' before cparm1',/)
         call cparm1 (nchan,ntarg,x(etarg),ix(lchl),ix(ichl),x(ethr),e0,
     x                enum,x(sce),ix(noch),x(cc))
*        write (iout,1766)
1766     format(' after cparm1',/)
*
*  calculate estimated final radius at each energy raf(ie) and its
*  maximum value rafmax
*
         call rafmx (nchan,x(ethr),enum,x(sce),x(cc),rafest,rafmax,
     x               x(raf),rmatr,ion,lamax)
*        print *,'after rafmx'
*
         lmax = maxl(nchan,ix(lchl))
*
*  set store for gailitis wave function package
*
         pt0 = last
         call galptr (nchan,iasy,lmax,pt0,mcore)
*        print *,'after galptr'
         last = pt0
         if (last .gt. mcore) then
            write (iout, 1010) last,mcore
*           stop
         end if
*
*  index the channels according to target spins to enable the rmatrix
*  to be divided into spin dependent parts which can be propagated
*  separately - used for the BBM propagator to save space and time
*  not implemented for the LW propagator which is highly vectorised
*  over channels.
*
            ispch = wpadti (last + 1)
            ichsp = ispch + nchan
            last = iadtwp (ichsp + nchan + 1)
*
            call split (nchan,ix(ichl),ntarg,ix(starg),ix(ispch),
     x                  ix(ichsp),nschan,nspin,ns)
*           print *,'after split'
*
*  set store for propagators and define the sector ends
*
            asect = last
            mcor1 = mcore - last + 1
            pt0 = last
            ilast = wpadti (last +1)
*
*  the radius at which LW propagator starts is raflw
*  no store is allotted for BBM if raflw = rmatrix radius
*  otherwise raflw is redefined to be at the end of the last BBM sector
*  note that raflw can be greater than rafmax so that BBM can be
*  initialised further than the first estimate of the propagation radius
*
            raflw = raflw1
C            raflw = min(raflw,rafmax)
*            print *,'raflw before prppnt = ',raflw
*            print *,'last before prppnt = ',last
           call prppnt(rfirst,raflw,nschan,ebig,ion,nranmx,idscbb,idbbw,
     x                  nampx,x(last),ix(ilast),nleg,ixmax,delta,pt0,
     x                  mcor1)
*            print *,'after prppnt',last,mcore
            last = pt0
            if (last .gt. mcore) then
               write (iout, 1010) last,mcore
*              stop
            end if
            nran = nranmx
            rafin = rafmax
            if ( nran .gt. nbmx ) then
               write (iout,1050) nran, nbmx
               stop
            endif
*
         if (nran .gt. 0) then
*
*  open files for BBM propagator
*
         do 55 i = 1,ns
            mch = nschan(i)
            mhd = 2 * mch * ixmax
            mamp = mch * mhd
            irecl = mamp + mamp + mhd
            irecl = irecl * 8
            open (idscbb(i),file = frmp(i),access = 'direct',
     x         recl = irecl,status = 'unknown')
 55      continue
         end if
*
            pt0 = last
*
*  set store for LW propagator
*
            call lwint(nchan,mcore,pt0)
            last = pt0
         if (last .gt. mcore) then
            write (iout, 1010) last,mcore
*           stop
         end if
*DEB
*        print *,'after store for LW:',last,mcore
*
*  open file for light-walker propagator
*
         irecl = 2 + nchan + 3*nchan*nchan
         irecl = irecl * 8
         open (idiscr,file = flwp,access = 'direct',
     x         recl = irecl,status = 'unknown')
*
      if (diaflg(8) .eq. 1) then
         write (iout,1060) rafin,iasy,apeps,degeny,iapflg
      end if
*
         ipr = 0
         rlast = raflw
         ming = mingap(x(ethr),nchan,degeny)
         cuplm = cuplim * ming
         if (nran .ne. 0) then
            do 60 i = 1, nran
               nbasis(i) = 0
   60       continue
            call spropi (x(ethr),nchan,nschan,ix(ispch),nran,nbasis,
     x                x(asect),bbloch,rpflg,idscbb,ipr,ismax,ion,
     x                ix(lchl),x(cf),nampx,x)
            if (diaflg(7) .eq. 1) then
               write (iout,1070) ipr,rfirst,raflw,rpflg
               write (iout,1090) nran,nchan,nleg,idscbb,ixmax,ismax
            end if
         end if
*
        if (raflw .lt. rafin) then
           if (diaflg(7) .eq. 1) then
              write (iout,1061)raflw,rafin,beta
           end if
            nstart = 1
            call slwpot (nstart,x(vmat),x(enrg),x(tc),x(tp),x(tsave),
     x                  nchan,nchan,lwdb(1),lwdb(2),raflw,rafin,
     x                  lwsect,lwkeep,lwran,beta,cuplm,x(space),
     x                  x(tev),x(tpev),ismax,ion,ix(lchl),x(cf),
     x                  x(ethr),idiscr,x(cupest))
            nfinal = nstart
         else
            nfinal = 1
            lwsect(1) = raflw
            lwkeep(1) = nchan
         end if
        if (diaflg(12) .eq. 1) then
           if (raflw .gt. rfirst) write (iout,1110) raflw
           if ( nfinal .gt. 1)   write (iout,1120) lwsect(nfinal)
           write (iout,1140)
        end if
*
*  set store for kmat - use common/kp/
*
         pt0 = last
         call kptr (nchan,pt0,mcore)
         last = pt0
         kmat = last
         last = kmat + ( nchan * nchan + nchan ) / 2
         if (last .gt. mcore) then
            write (iout, 1010) last,mcore
*           stop
         end if
*
*  set store for eigpha and etable
*
         c = 0
         if (prtflg(2) .ne. 0) then
            wrk = last
            phz = wrk + 3 * nchan
            phzv = phz + enum
            c = phzv + nchan
            last = c + ncn2
            if (last .gt. mcore) then
               write (iout, 1010) last,mcore
*              stop
            end if
            do 40 i = 1,enum
               x(phz + i -1) = zero
 40         continue
         end if
*
*  set store for tmat
*
         if ( prtflg(3) .ne. 0 .or. prtflg(4) .ne. 0
     x        .or. prtflg(5) .ne. 0 .or. prtflg(7) .ne. 0
     x        .or. diaflg(11) .ne. 0) then
            tr = last
            ti = tr + ncn2
            a = ti + ncn2
            ind = wpadti(a + ncn2 + 1)
            last = iadtwp(ind + nchan + 1) + 1
            if (c .eq. 0) then
               c=last
               last = c + ncn2
            end if
         endif
         if ( last .gt. mcore ) then
            write (iout,1010) last,mcore
*           stop
         endif
*
*  set store for crostn and totcol - calculate statistical weights
*
         if ( prtflg(4) .ne. 0 .or. prtflg(5) .ne. 0 .or.
     x        diaflg(11) .ne. 0 ) then
            sk = last
            sw2 = sk + ntarg
            ji = wpadti(sw2 + ntarg + 1)
            last = iadtwp(ji + enum + 1)
            call statw (lrgl2,nspn2,ntarg,sw1,x(sw2),ix(ltarg),
     x                  ix(starg))
         end if
*
         if ( last .gt. mcore ) then
            write (iout,1010) last,mcore
*           stop
         endif
*DEB
*        print *,' after store for crostn:',last,mcore
         if ( diaflg(15) .ne. 0 ) write (iout,1040) last,mcore
*        print *,'before phased:',last,mcore
         call phased (prtflg,nmt,nspn2,lrgl2,npty2,nl,x(fl),enum,
     x                x(sce),ix(noch),x(cfbut),x(rmat),
     x                ix(lchl),nstat,nchan,nschan,ix(ispch),ix(ichsp),
     x                x(cf),x(eig),x(wmat),e0,ibuttl,x(etarg),
     x                ix(ichl),x(cc),ion,zion,lamax,rfirst,adwm,
     x                diaflg,x(kmat),
     x                bbloch,x(wrk),ix(ind),ncol,x(phz),x(phzv),
     x                x(f),x(fp),x(et),nbasis,nampx,nranmx,rpflg,
     x                apeps,degeny,iasy,lmax,ewron,icon,quad,title,
     x                x(ethr),raflw,x(raf),lwdb,idiscr,x(a),
     x                x(c),x(xslpi),x(xsn),x(xsnl),ntarg,
     x                ix(nltarg),ix(ji),x(tr),
     x                x(ti),nxsn,ix(iind),
     x                ix(jind),icase,ncases,sw1,x(sw2),x(sk),
     x                x(asect),lwsect,lwkeep,keep,nfinal,lwran,iapflg,
     x                nbutx,x(fkn),x(uk),beta,cuplm,deplim,ncrit,rfact,
     x                gfact,maxtry,x(radius),ix(rstats),rafix,x,ix)
         close (idiscr,status = 'DELETE')
         close (idscbb(1))
         close (idscbb(2))
 50   continue
*
      call wtcas (ix(lrgs),ix(lrgl),ix(lrgpi),ncases)
*
*   print line strengths
*
      if (.not. jj .and. ismax .gt. 0) call wrslin(ntarg,xsn,x(sline))
*
*   print total x sections/collision strengths if required
*
      if (prtflg(4) .gt. 1  .or.  prtflg(5) .gt. 1 )then
         ixsn = 2
         if (prtflg(4) .gt. 1) then
            write(iout,1170)
         else
            write(iout,1180)
         end if
         if (enum1 .gt. 0) call pxsntb (enum1,nxsn,
     x       x(xsn),x(sce1),ixsn,ix(iind),ix(jind),ncol,nphys)
         if (enum2 .gt. 0) call pxsnl (nl,x(fl),enum2,
     x       nxsn,x(xsnl),x(sce2),ixsn,ix(iind),ix(jind),ncol,nphys)
      end if
      return
 1000 format (//, ' program farm', 10x,
     x         ' CPC version  (1994)', //)
 1010 format (' storage overflow ( cntrl ) : last =',i10,' icore =',
     x          i10/'  increase value of icore in main program')
 1020 format (//1x,a80//)
 1040 format (' core monitor ( cntrl ) core used =',i10,' icore =',i10)
 1050 format (' propagator subranges, nran =',i4,' exceeds maximum ',i4)
 1060 format (/,5x,'expansion radius,             rafin  =',d16.8,/,
     x        5x,'terms in series,              iasy   =',i16,/,
     x        5x,'convergence criterion,        apeps  =',d16.8,/,
     x        5x,'channel degeneracy criterion, degeny =',d16.8,/,
     x        5x,'print flags,                  iapflg =',4i4,/)
 1061 format (' r-matrix propagation : light-walker method',//,
     x        ' parameters for r-matrix propagation',/,
     x        ' initial radius=',f15.5,10x,'final radius=',f15.5,/,
     x        ' beta=',d15.5,/)
 1070 format (' r-matrix propagator :',/,5x,'ipr   =',i3,' rfirst =',
     x        d16.8,' rlast =',d16.8,/,5x,'rpflg  :',6i3)
 1090 format (' nran   =',i4,' nchan  =',i4,' nleg   =',i4,/,
     x        ' idscbb  =',2i4,' ixmax  =',i4,' ismax  =',i4)
 1110 format (/' Burke Baluja Morgan propagator initialised to r =',
     x           f12.6,' scaled Bohr')
 1120 format (/' Light Walker propagator initialised to r =',f12.6,
     x           ' scaled Bohr')
 1140 format (/' Light Walker propagator will be used to extend the',
     x         '  radius')
 1150 format (' **** error ****  cross sections AND collision ',
     x        'strengths cannot be calculated together')
 1160 format (' **** error **** top up output is required',/,
     x        ' prtflg(4) OR prtflg(5) must be non zero')
 1165 format (' **** error **** top up output is required',/,
     x        ' prtflg(7) must be zero')
 1170 format (//,15x,'TOTAL CROSS-SECTIONS in units of a0**2'/,
     x           15x,'**************************************',/)
 1180 format (//,15x,'TOTAL COLLISION STRENGTHS'/,
     x           15x,'*************************',/)
 2000 format (/,15x,'****************',/,15x,'s l p =',3i3,/
     x           15x,'****************', /)
 2010 format (/,' Internal region data read from file ',a128,/
     x        '   nelc =',i4,' nz=',i4,' lrang2=',i4,/
     x        '   number of target states =',i4,
     x        ' highest multipole =',i2/
     x        '   rmatrix radius =',d16.8,' Bohr',//)
 2011 format (' opened output file ',a128)
 2020 format (' number of discrete energies away from thresholds = ',i4)
 2030 format (' ***** scaled units are used by the code:',/
     x        ' for positively charged ions energies are scaled',/
     x        ' by a factor 1/z**2 and radii by z',/
     x        ' where z is the charge on the ion (nz-nelc)',/,
     x        ' for neutrals the scaling factor is 1 ',/
     x        ' ***** data read from namelists is in scaled units:',/
     x        ' energies in scaled rydbergs and radii in scaled Bohr'//)
       end
      subroutine cparm1 (nchan,ntarg,etarg,lchl,ichl,eth,e0,
     x                   enum,sce,noch,cc)
*
* $Id: cparm1.f,v 2.1 94/09/29 18:40:01 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 dfloat,atry,etarg,eth,e0,sce,cc,eps,zero
      parameter ( atry=2.0d0,zero=0.0d0 )
      common /io/ iin,iout
      dimension etarg(ntarg),lchl(nchan),ichl(nchan),eth(nchan),
     x          sce(enum),noch(enum),cc(nchan)
*
*     cparm : channel parameters eth, lch
*
      e0 = etarg(1)
c      write (iout,1000) e0
      do 30 i = 1, nchan
         ii = ichl(i)
         eth(i) = atry * ( etarg(ii) - e0 )
c         write (iout,1010) i,lchl(i),eth(i)
         cc(i) = dfloat (lchl(i)*(lchl(i) + 1))
 30   continue
      do 60 ie = 1,enum
         nopen = 0
         do 70 i = 1,nchan
            eps = sce(ie) - eth(i)
            if (eps .gt. zero) nopen = nopen + 1
  70     continue
         noch(ie)=nopen
  60  continue
      return
*
c 1000 format (/,' cparm1 :',/,' e0 =',d16.8,)
c 1010 format (' i =',i6,' lch =',i6,' eth =',d16.8)
      end
      subroutine crostn (e,ntarg,nltarg,et,sw1,sw2,a,xslpi,xsn,
     x                   idbug,nxsn,sk,zion,ji,iind,jind,tr,ti,
     x                   nopen,ien,nl,nume2,xsnl,il)
*
* $Id: crostn.f,v 2.1 94/09/29 18:40:18 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 e,et,a,xslpi,xsn,xsnl,tr,ti,an2,zero,one,half,c1,c2,
     x       sw1,sw2,cros,pi,cx,cxp,cxq,dx,sk,flion2
      parameter ( an2=0.28002836d+00,zero=0.d0,one=1.0d0,half=0.5d0,
     x            pi=3.1415926535897932d+00 )
      dimension nltarg(ntarg),et(*),tr(*),ti(*),ji(*),a(*),sw2(ntarg),
     x          xslpi(nxsn),xsn(nxsn),xsnl(nxsn,nume2,nl),
     x          sk(ntarg),iind(nxsn),jind(nxsn)
*kbion
      real*8 x1s(1000),x2s(1000),x2p(1000)
*kbion
      common /io/ iin,iout
*
*     crostn : total cross sections from the t-matrix
*              an2 is the conversion factor from a0**2 to angstroms**2
*
*  calculate energy scaling factor ready to calculate cross-sections
*
      ion = zion
      if (ion .eq. 0) then
         flion2 = one
      else
         flion2 = dfloat (ion * ion)
      end if
      nop2 = nopen * nopen
      do 12 i = 1, nop2
         a(i) = tr(i)**2 + ti(i)**2
   12 continue
      j = 1
      do 14 i = 1, ntarg
         if ( et(j) .gt. zero ) then
            sk(i) = et(j) * flion2
            nlt = i
            j = j + nltarg(i)
         endif
   14 continue
      ix = 0
      c1 = half * pi * sw1
      nf = 0
      do 60 i = 1, nlt
         ni = nf + 1
         nf = nf + nltarg(i)
         mf = 0
         jup = i
         if (ion .gt. 0) jup = i - 1
         do 40 j = 1,jup
            ix = ix + 1
            mi = mf + 1
            mf = mf + nltarg(j)
            c2 = sw2(j) * sk(j)
            cros = zero
            do 30 nc = ni, nf
               do 20 mc = mi, mf
                  nmc = nopen * ( mc - 1 ) + nc
                  cros = cros + a(nmc)
   20          continue
   30       continue
            xslpi(ix) =  c1 * cros / c2
   40    continue
         ji(ien) = ix
   60 continue
*     debug print of cross-sections for partial wave
      if ( idbug .gt. 0 .and. ix .gt. 0 ) then
         write (iout,1000)e
         ix = 0
*kbion
         do 117 i=1,ntarg
          x1s(i) = 0.0
          x2s(i) = 0.0
          x2p(i) = 0.0
117      continue
*kbion
         do 90 i = 1, nlt
            jup = i
            if (ion .gt. 0) jup = i - 1
            do 80 j = 1, jup
               ix = ix + 1
               cxp = xslpi(ix)
               cx =  cxp / pi
               cxq = an2 * cxp
               dx = sw2(j) * sk(j) * cx
*kbion
               if(j.eq.1) x1s(i) = cxp
               if(j.eq.2) x2s(i) = cxp
               if(j.eq.3) x2p(i) = cxp
*kbion
               if (cxp .ne. zero) then
               write (iout,1010) j,i,cx,cxp,cxq,dx
               end if
   80       continue
   90    continue
*kbion
         write(81,8181) e
         write(82,8181) e
         write(83,8181) e
         write(81,8182) (x1s(i),i=1,ntarg)
         write(82,8182) (x2s(i),i=1,ntarg)
         write(83,8182) (x2p(i),i=1,ntarg)
8181     format(1pd18.10)
8182     format(1p5d16.8)
*kbion
      endif
      if (il .le. 0) then
         do 100 i = 1,ix
            xsn(i) = xsn(i) + xslpi(i)
 100     continue
      else
         do 101 i = 1,ix
            xsnl(i,ien,il) = xsnl(i,ien,il) + xslpi(i)
 101     continue
      end if
      return
*
 1000 format (/' cross sections at energy ',d12.5,'  (scld ryd)',
     x        /'   i   j      pi*a0**2',9x,
     x        'a0**2',11x,'ang**2',4x,'collision strength',/)
 1010 format (2i4,4d16.6)
c 1030 format (1x,6(i3,'-',i3,5x))
c 1040 format (1x,6d12.5)
      end

      subroutine degtar (nelc,nz,ntarg,etarg,ltarg,starg,degeny,
     x                   prtflg,nmt,title)
*
*  make thresholds closer than degeny truly degenerate
*  print the new energy levels if changed and write headings to
*  output files for k-matrix, t-matrix, cross-section, collision strength
*  and top-up.
*
      implicit integer (a-z)
      real*8 etarg,zero,atry,degeny,etargr,d0
      character*80 title
      character*20 name
      logical change
      parameter ( zero=0.0d0, atry=2.0d0)
      common /io/ iin,iout
      dimension etarg(ntarg),ltarg(ntarg),starg(ntarg),prtflg(*),
     x          nmt(*),name(20)
      data name(1),name(3),name(4),name(5),name(6),name(7)
     x  /' K-MATRIX',' T-MATRIX',' X-SECTION ',' OMEGA ',
     x   ' X-SECT-FOR-TOPUP',' OMEGA-FOR-TOPUP '/
*
      change = .false.
      do 1 i = 2,ntarg
*    prevent vectorization of this loop
      write(iout,*)
         if (abs(etarg(i)-etarg(i-1))*atry .le. degeny) then
            etarg(i) = etarg(i-1)
            change = .true.
         end if
 1    continue
      if (change) then
         write (iout,1000)
         d0 = etarg(1)
         do 10  i=1,ntarg
             etargr = atry * (etarg(i) - d0)
             write (iout,1010) i,ltarg(i),starg(i),etargr
 10      continue
      end if
cz
cz      nutar=88
cz      open(nutar,file='farm_tar',form='UNFORMATTED')
cz
cz      if (prtflg(1) .eq. 1  .or.  prtflg(1) .eq. 3) then
cz         write(nutar) name(1),title                         ! nmt(1)
cz         write(nutar) nelc,nz,ntarg,(etarg(i),i=1,ntarg),
cz     x                (ltarg(i),i=1,ntarg),(starg(i),i=1,ntarg)
cz      end if
cz      if (prtflg(3) .eq. 1  .or.  prtflg(3) .eq. 3) then
cz         write(nutar) name(3),title                             ! nmt(3)
cz         write(nutar) nelc,nz,ntarg,(etarg(i),i=1,ntarg),
cz     x                (ltarg(i),i=1,ntarg),(starg(i),i=1,ntarg)
cz      end if
cz      if (prtflg(4) .eq. 1  .or.  prtflg(4) .eq. 3) then
cz         write(nutar) name(4),title                            ! nmt(4)
cz         write(nutar) nelc,nz,ntarg,(etarg(i),i=1,ntarg),
cz     x                (ltarg(i),i=1,ntarg),(starg(i),i=1,ntarg)
cz      end if
cz      if (prtflg(5) .eq. 1  .or.  prtflg(5) .eq. 3) then
cz         write(nutar) name(5),title
cz         write(nutar) nelc,nz,ntarg,(etarg(i),i=1,ntarg),       ! nmt(5)
cz     x                (ltarg(i),i=1,ntarg),(starg(i),i=1,ntarg)
cz      end if
cz      if (prtflg(6) .eq. 1  .or.  prtflg(6) .eq. 3) then
cz         if (prtflg(4) .gt. 0) then
cz            write(nutar) name(6),title                           ! nmt(6)
cz         else if (prtflg(5) .gt. 0) then
cz            write(nutar) name(7),title
cz         end if
cz         write(nutar) nelc,nz,ntarg,(etarg(i),i=1,ntarg),
cz     x                (ltarg(i),i=1,ntarg),(starg(i),i=1,ntarg)
cz      end if
cz

*
 1000 format(20x, 'target states modified by degeny'/
     x       20x, '********************************',//
     x       10x,'index',5x,'total l',3x,'(2*s+1)',8x,'energy'/
     x       43x,'scaled ryd')
 1010 format(3x,3i10,7x,f12.6)
      return
      end


      subroutine eigpha (iph,phz,phzv,nchan,ien,nc,kmat,wrk)
*
* $Id: eigpha.f,v 2.1 94/09/29 18:39:39 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 phz,phzv,kmat,wrk,vc,zero
      parameter ( zero=0.d0 )
      common /io/ iin,iout
      dimension kmat(*),phz(*),phzv(nchan),wrk(3*nchan),vc(1)
*
*     eigenphases of k-matrix
*
      call DSPEV ( 'N', 'U', nc, kmat, phzv, vc, nchan, wrk, info )
      if ( info .ne. 0 ) then
         write (iout,1000) info
         stop
      endif
      do 80 no = 1, nc
         phzv(no) = atan( phzv(no) )
   80 continue
      if ( iph .ge. 2 ) then
         write (iout,1010) (phzv(no),no = 1,nc)
         istart = 1
c         do 85 no=1,nc
c            write (iout,1020) (vc(i),i=istart,istart+nc-1)
c            istart=istart+nc
c 85      continue
      end if
      phz(ien) = zero
      do 90 no = 1,nc
         phz(ien) = phz(ien) + phzv(no)
 90   continue
      return
*
 1000 format (' eigpha : error return by dspev, ierr =', i3)
 1010 format (' eigenphases :',/,(6f12.6))
 1020 format (' eigenvectors :',/,(6f12.6))
      end
      subroutine etable (icon,quad,title,nume,phz,en,nspn2,lrgl2,
     x                   npty2,nfts)
*
* $Id: etable.f,v 2.1 94/09/29 18:39:41 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 quad,quad1,phz,en,pi,one,e,dl,dlx,dlm,dlp,dlx1,
     x       dlx2,dlx3
      character*80 title
      parameter ( pi=3.1415926535897932d+00,one=1.0d0)
      common /io/ iin,iout
      dimension phz(nume),en(nume)
*
*     etable : smoothing of eigenphase sum:
*              icon = 0  no smoothing
*                   = 1  add pis to minimize energy differences
*
*     second-order divided differences, assume first three points smooth
*
      if ( icon .eq. 1 ) then
         do 40 ie = 1, nume
            e = en(ie)
            dl = phz(ie)
            if ( ie .eq. 1 ) then
               phz(ie) = dl + quad * pi
               quad1 = quad
            else
                  dlx = dl + quad1 * pi
                  dlm = dlx - pi
                  dlp = dlx + pi
                  dlx1 = abs(dlx-phz(ie-1))
                  dlx2 = abs(dlm-phz(ie-1))
                  dlx3 = abs(dlp-phz(ie-1))
                  if ( dlx2 .lt. dlx1 ) then
                     if ( dlx3 .ge. dlx2 ) then
                        phz(ie) = dlm
                        quad1 = quad1 - one
                     else
                        phz(ie) = dlp
                        quad1 = quad1 + one
                     endif
                  elseif ( dlx3 .lt. dlx1 ) then
                     phz(ie) = dlp
                     quad1 = quad1 + one
                  else
                     phz(ie) = dlx
                  endif
               endif
   40    continue
      endif
*
*     print table
*
      write (nfts,1000) title
         do 360 ie = 1, nume
            e = en(ie)
               write (nfts,1010) nspn2,lrgl2,npty2,e,phz(ie)
  360    continue
      return
*
 1000 format (1h1,1x,a72,/,30x,'table of eigenphases',//,
     x       1x,'  2s+1   L   Pi      energy      phase sum')
 1010 format (1x,3i5,2d14.5)
      end

      subroutine evals (ntarg,nume,ne,esc,etarg,abvthr,belthr,sce,
     x                  enum,ebig,scest,norder)
*
* $Id: evals.f,v 2.1 94/09/29 18:39:42 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 esc,etarg,abvthr,belthr,sce,e,zero,etest,einc,estart,ebig,
     x        e1,atry,scest,fix,desmat
      parameter ( zero=0.0d0, atry=2.0d0, fix=1.0d6, desmat=1.0d-8 )
      common /io/ iin,iout
      dimension ne(*),esc(*),etarg(ntarg),sce(nume),scest(nume),
     x          norder(nume)
*
      ebig = zero
      enum = 0
      e1 = etarg(1)
      if ( ne(1) .eq. 0 ) go to 61
      neig = 10
      iec = 1
      do 30 ieg = 1, neig
         neigs = ne(ieg)
         if (neigs .eq. 0) go to 40
         estart = esc(iec)
         if (abs(neigs) .eq. 1 ) then
            einc = zero
            neigs = abs(neigs)
         else if (neigs .lt. -1) then
            neigs = abs(neigs)
            einc = (esc(iec+1) - estart)/(neigs-1)
         else
            einc = esc(iec+1)
         endif
         iec = iec + 2
         e = estart
         if ( e .lt. zero ) then
            write (iout,1010)
            stop
         endif
         do 20 ie = 1, neigs
            do 10 i = 1, ntarg
               etest = e - (etarg(i) - e1) * atry
               if ( etest .ge. zero ) then
                  if ( etest .le. abvthr ) go to 21
               else
                  if ( - etest .le. belthr ) go to 21
               endif
   10       continue
            enum = enum + 1
            scest(enum) = e
 21         e = e + einc
 20      continue
 30   continue
*
*   since the meshes could overlap etc. reorder the energies in
*   ascending order.  Duplicating energies are set to a very large
*   value, fix, ready to be pruned out
*
 40   call order(scest,enum,norder,fix,desmat)
      do 60  i = 1,enum
         sce(i) = scest(norder(i))
 60   continue
      maxnum = 0
      do 70 i = enum,1,-1
         if (sce(i) .lt. fix) then
            maxnum = i
            go to 75
         end if
 70   continue
 75   enum = maxnum
      ebig = sce(enum)
 61   ebig = max (ebig,(etarg(ntarg)-e1)*atry)
*
*   ebig may need modifying for bound state calculations
      return
 1000 format (' evals : input error, ne(1) =',i6,' nume =',i6)
 1010 format (' evals : input error, e =',d16.8,' negative')
      end
      subroutine kmat (nchan,ra,bsto,aa,bb,iwrk,nopen,f,fp,rmat,
     x                 akmat)
*
* $Id: kmat.f,v 2.1 94/09/29 18:40:30 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 aa,akmat,bb,bsto,df,f,fp,ra,rmat,kscal,temp,
     x       scalx,rscalx
      dimension rmat(*),akmat(*),f(nchan,nchan,2),fp(nchan,nchan,2),
     x          aa(nchan,nchan),bb(nchan,nchan),iwrk(nchan)
*     x         ,kscal(nchan)
      common /mach/ scalx,rscalx
*
      k = 2
      do 20 j = 1, nchan
         if ( j .gt. nopen ) k = 1
         do 10 i = 1, nchan
            aa(i,j)   = f(i,j,k)
   10    continue
 20   continue
*
      do 50 j = 1, nchan
         kj = 2
         if ( j .gt. nopen ) kj = 1
         do 40 k = 1, nchan
            df = ra * fp(k,j,kj) - bsto * f(k,j,kj)
            do 30 i = 1, nchan
               kix = max(i,k)
               kin = min(i,k)
               kk = ( kix * ( kix - 1 ) ) / 2 + kin
               aa(i,j) = aa(i,j) - rmat(kk) * df
   30       continue
   40      continue
* mod by ts 7-jan-94 --- scale closed-open k-matrix
         if ( j .gt. nopen ) then
            temp = 1./max(aa(j,j),scalx)
*            kscal(j) = temp
            do 45 i = 1, nchan
               aa(i,j) = aa(i,j) * temp
 45         continue
         endif
* end mod ---
   50 continue
      nc = nchan
      do 80 j = 1, nopen
         do 70 i = 1, nchan
            bb(i,j) = - f(i,j,1)
   70    continue
   80 continue
      do 110 j = 1, nopen
         do 100 k = 1, nchan
            df = ra * fp(k,j,1) - bsto * f(k,j,1)
            do 90 i = 1, nchan
               kix = max(i,k)
               kin = min(i,k)
               kk = ( kix * ( kix - 1 ) ) / 2 + kin
               bb(i,j) = bb(i,j) + rmat( kk) * df
   90       continue
  100    continue
  110 continue
         if ( nc .eq. 1 ) then
            bb(1,1) = bb(1,1) / aa(1,1)
         else
            call sleqn (aa,nchan,nc,bb,nopen,iwrk)
         endif
         k = 0
         do 130 j = 1, nopen
            do 120 i = 1, j
               k = k + 1
               akmat(k) = bb(i,j)
 120        continue
* mod by ts 7-jan-94 --- rescale closed-open k-matrix
* rescaling is not necessary here as the array bb is not used
* again
*            do 125 i = nopen+1, nchan
*               bb(i,j) = bb(i,j) * kscal(i)
* 125        continue
* end mod ---
  130     continue
*
      return
      end

      subroutine kptr (nchan,pt0,mcor)
*
* $Id: kptr.f,v 2.1 94/09/29 18:39:44 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      common /io/ iin,iout
      common /kp/ bb,iwrk,bigvc
*
      bigvc = pt0
      nchsq = nchan * nchan
      bb = bigvc + nchsq
      iwrk = wpadti(bb +  nchsq) +1
      last = iadtwp(iwrk + nchan + 1)
      if ( last .gt. mcor ) then
         write (iout,1000) last,mcor
         stop
      endif
      pt0 = last
      return
 1000 format (' kptr : last =',i10,' core =',i10)
      end
      subroutine llims(ncases,lrgl,lstart,lfin)
*
* $Id: llims.f,v 2.1 94/09/29 18:41:05 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      dimension lrgl(ncases)
*
      lstart = 999
      lfin = 0
      do 10 il = 1,ncases
         l = lrgl(il)
         lstart = min(lstart, l)
         lfin = max(lfin, l)
 10   continue
      return
      end
      subroutine lstren(ntarg,nltarg,ltarg,starg,nxsn,sline,cf,
     x                 nchan,lamax,lchl,lrgl2,finis)
      implicit integer(a - z)
      real*8 sline,cf,dfloat,wabs1,c,zero,quartr
      logical finis
      dimension nltarg(ntarg),ltarg(ntarg),starg(ntarg),sline(nxsn),
     x          cf(nchan,nchan,lamax),lchl(nchan)
      data zero/0.0d0/,quartr/0.25d0/
*  calculate line strength in LS coupling

      nc2=nltarg(1)
      do 10 it=2,ntarg
         if(nltarg(it).gt.0) then
         nc1=nc2+1
         nc2=nc2+nltarg(it)
         mc2=0
         do 20 jt=1,it-1
            if(nltarg(jt).gt.0) then
            mc1=mc2+1
            mc2=mc2+nltarg(jt)
            k=((it-1)*(it-2))/2+jt
            if(sline(k).lt.zero)then
            do 30 nc=nc1,nc2
               do 40 mc=mc1,mc2
                  if(cf(nc,mc,1).ne.zero) then
                  c=cf(nc,mc,1)/
     x              wabs1(ltarg(it),ltarg(jt),lchl(nc),lchl(mc),lrgl2)
                  c=c*c
                  sline(k)=.25*starg(it)*c/dfloat(max(lchl(nc),
     x                     lchl(mc)))
                  goto 20
                  end if
   40             continue
   30          continue
            sline(k) = zero
            end if
            end if
   20    continue
         end if
   10 continue
      finis = .true.
      k = 0
      do 50 it = 2,ntarg
         si = starg(it)
         do 60 jt = 1,it - 1
            k = k + 1
            if (starg(jt) .ne. si) sline(k) = zero
            if (sline(k) .lt. zero) finis = .false.
 60      continue
 50   continue
      return
      end

      subroutine lvals (ncases,lrgl,lrgle,nl,fl)
*
* $Id: lvals.f,v 2.1 94/09/29 18:41:03 vmb Exp Locker: vmb $
*
      real*8 fl,f
      common /io/ iin, iout
      dimension lrgl(ncases),fl(*)
      nl = 0
      do 10 n=1,ncases
         l = lrgl(n)
         if (l .ge. lrgle) then
            f = dble(l)
            do 20 i = 1,nl
               if (f .eq. fl(i)) go to 10
 20         continue
            nl = nl + 1
            fl(nl) = f
         end if
 10   continue
      write (iout,1000) nl
      return
 1000 format (' number of l values in top up scheme, nl =',i3)
      end

      function maxl (nchan,l2p)
*
* $Id: maxl.f,v 2.1 94/09/29 18:39:46 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      dimension l2p(nchan)
      lmax = -1
      do 10 i = 1, nchan
         lmax = max(lmax,l2p(i))
   10 continue
      maxl = lmax
      return
      end
      function mingap (eth,nchan,degeny)
*
* $Id: mingap.f,v 2.1 94/09/29 18:39:35 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 eth,degeny,ming,gap,mingap
      dimension eth(nchan)
*
*     calculate the minimum gap between target thresholds
*
      ming = eth(nchan) - eth(1)
      do 10   i = 2,nchan
         gap = eth(i) - eth(i-1)
         if (gap .gt. degeny) ming = min(ming, gap)
 10   continue
      mingap = ming
      return
      end
      function necnt (ne)
*
* $Id: necnt.f,v 2.1 94/09/29 18:39:47 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      common /io/ iin,iout
      dimension ne(10)
*
*     necnt = number of impact energies
*
      ienum = 0
      do 80 ieg = 1, 10
         ienum = ienum + abs(ne(ieg))
   80 continue
      necnt = ienum
      write (iout,1000) necnt
      return
 1000 format (/' necnt : number of energies =',i6)
      end
      subroutine order(en,nchan,norder,fix,desmat)
c
c --- define norder, which points to values of en in ascending order.
c
c --- input parameters ...
c     (en(m),m=1,nchan) = energies in arbitrary order;
c --- output parameter ...
c     (norder(m),m=1,nchan) = the position of the m-th biggest energy
c                             in the en array.
c
      implicit integer (a-z)
      real*8 en,fix,desmat,x
      dimension en(nchan),norder(nchan)
c
      norder(1)=1
      if(nchan.eq.1) return
      do 9 m=2,nchan
         j=m
         x=en(j)
         j1=j-1
         do 7 i=1,j1
            if (j.ge.m) then
               i1=norder(i)
               if (abs(x - en(i1)) .le. desmat) then
                  en(j) = fix
                  x = fix
               end if
               if (x.ge.en(i1)) go to 7
            end if
            j=j-1
            norder(j+1)=norder(j)
    7    continue
         norder(j)=m
    9 continue
      return
      end
      subroutine parcol(e,sw1,a,nchan,ichl,lchl,tr,ti,nopen)
*
* $Id: parcol.f,v 2.1 94/09/29 18:40:28 vmb Exp Locker: vmb $
*
*
*   partial collision strengths
*
      implicit integer (a-z)
      real*8 e,sw1,a,tr,ti,half
      parameter ( half=0.5d0)
      common /io/ iin,iout
      dimension ichl(nchan),lchl(nchan),a(*),tr(*),ti(*)
*
      nop2 = nopen * nopen
      do 10 i = 1, nop2
         a(i) = tr(i)**2 + ti(i)**2
   10 continue
*
         write (iout,1000)
         do 20 i = 1,nopen
         do 20 j = i,nopen
            k = nopen * (i-1) + j
            a(k) = a(k) * half * sw1
            write (iout,1020) e,ichl(i),lchl(i),ichl(j),
     x                        lchl(j),a(k)
 20      continue
      return
*
 1000 format (/,' Partial collision strength :',/
     x        10x,'energy',10x,'TI LI    TJ LJ         OMEGA')
 1020 format (5x,d14.6,5x,2i3,3x,2i3,5x,d14.6)
      end
      subroutine phased (prtflg,nmt,nspn2,lrgl2,npty2,nl,fl,
     x                   nume,sce,noch,cfbut,rmat,lchl,nstat,nchan,
     x                   nschan,ispch,ichsp,cf,eig,wmat,e0,ibuttl,
     x                   etargs,ichl,cc,ion,zion,ismax,rmatr,adwm,
     x                   diaflg,kmat,bbloch,
     x                   wrk,ind,ncol,phz,phzv,f,fp,et,nbasis,nampx,
     x                   nranmx,rpflg,apeps,degeny,iasy,lmax,ewron,
     x                   icon,quad,title,ethr,raflw,raf,lwdb,idiscr,
     x                   a,c,xslpi,xsn,xsnl,
     x                   ntarg,nltarg,ji,tr,ti,nxsn,
     x                   iind,jind,icases,ncases,
     x                   sw1,sw2,sk,asect,lwsect,lwkeep,keep,nfinal,
     x                   lwran,iapflg,nbutx,fkn,uk,beta,cuplim,deplim,
     x                   ncrit,rfact,gfact,maxtry,radius,rstats,rafix,
     x                   x,ix)
*
* $Id: phased.f,v 2.1 94/09/29 18:39:51 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      character*80 title
      real*8 sce,cfbut,rmat,cf,eig,wmat,e0,etargs,rmatr,e,
     x       etot,ee,half,ryev,bbloch,wrk,phz,phzv,f,fp,et,apeps,
     x       raflw,raf,
     x       degeny,ewron,quad,kmat,x,ethr,tr,ti,a,c,xslpi,xsn,
     x       xsnl,sw1,sw2,sk,fkn,uk,cc,radius,rafix,
     x       asect,lwsect,beta,rfact,gfact,
     x       cuplim,deplim,fl
      logical lwdb,adwm,firste
      parameter ( half=0.5d0,ryev=13.605d0)
      common /io/ iin,iout
      common /ie/ ee,ierp
      dimension et(*),sce(nume),cfbut(*),rmat(*),lchl(*),
     x          cf(*),eig(*),wmat(*),noch(*),wrk(*),phz(*),phzv(*),
     x          nbasis(*),etargs(*),ichl(*),kmat(*),f(*),fp(*),rpflg(6),
     x          x(*),ix(*),ind(*),ethr(*),tr(*),ti(*),a(*),c(*),
     x          xslpi(nxsn),xsn(nxsn,nume),xsnl(nxsn,nume,nl),
     x          iind(nxsn),jind(nxsn),fl(nl),nltarg(*),ji(*),sw2(ntarg),
     x          sk(*),iapflg(4),fkn(nbutx),uk(nbutx),cc(*),asect(*),
     x          lwsect(lwran),lwkeep(lwran),keep(lwran),lwdb(2),raf(*),
     x          diaflg(*),prtflg(7),nmt(7),nschan(2),
     x          ispch(nchan),ichsp(nchan),radius(maxtry),rstats(maxtry)
*
      ierp = iout
*DEB
*     print *,' started phased with ierp = ',ierp
      il = 0
      do 5 i = 1 ,nl
         if (lrgl2 .eq. nint(fl(i))) il = i
 5    continue
*
*  loop over scattering energies
*
      firste = .true.
      do 100 ie = 1, nume
         e = sce(ie)
         ee = ryev * e
         nopen = noch(ie)
         do 10 i = 1, nchan
            et(i) = e - ethr(i)
 10      continue
*
*  if there are no open channels or in the case of ions fewer than 2
*  open states skip to the next energy
*
         if (nopen .eq. 0) go to 100
         if (ion .gt. 0 .and. prtflg(1) .eq. 0) then
             do 15 i = 1,ntarg
                if (nltarg(i) .gt. 0) then
                    k = i
                    go to 16
                end if
   15        continue
             stop
   16        if (nltarg(k) .ge. nopen) go to 100
          end if
*
*  calculate r-matrix at the boundary
*
*  use surface amplitudes from H-file or asymptotic distorted wave
*  method
*
         if (.not. adwm) then
*           print *,' before rmatd'
            call rmatd (nstat,nchan,eig,wmat,e0,cfbut,ibuttl,ichl,
     X                  etargs,rmatr,e,rmat,fkn,uk)
*           print *,' after rmatd'
         else
            call rmatf1(nchan,lchl,et,lmax,ion,iasy,rmatr,rmat,x)
         end if
*
         if ( diaflg(14) .gt. 0 ) then
            etot = e0 + half * e
            write (iout,1010) ie,e,ee,etot,nopen,(et(i),i = 1,nchan)
         endif
         if ( diaflg(2) .gt. 0 ) then
            write (iout,1020) rmatr,e
            call wrttmt (rmat,nchan,ncol)
         endif
*
*  calculate the k-matrix
*
         call asyms (ie,e,rmat,kmat,nchan,nschan,ispch,ichsp,
     x               lchl,ion,ismax,cf,rmatr,bbloch,ethr,
     x               f,fp,et,nbasis,nampx,nranmx,rpflg,prtflg,
     x               ncol,apeps,degeny,iasy,lmax,ewron,nopen,
     x               raf,raflw,asect,lwsect,lwkeep,keep,nfinal,lwran,
     x               iapflg,diaflg,idiscr,lwdb,beta,cuplim,deplim,
     x               ncrit,rfact,gfact,maxtry,radius,rstats,rafix,x,ix)
*
         if (prtflg(1) .eq. 1  .or. prtflg(1) .eq. 3) then
            if (firste) then
               numen = nume-ie+1
               write (nmt(1)) nspn2,lrgl2,npty2,numen,nchan,
     x                       (ichl(ich),lchl(ich),ich=1,nchan)
            end if
            nop2 = (nopen * (nopen + 1))/2
            write (nmt(1)) e,nopen,nop2,(kmat(i),i=1,nop2)
         end if
*
*  calcuate the eigenphase sum
*
         if (prtflg(2) .ge. 1) then
            call dcopy ((nopen*nopen+nopen)/2,kmat,1,c,1)
            call eigpha (prtflg(2),phz,phzv,nchan,ie,nopen,c,wrk)
         end if
*
*  calculate the t-matrix
*
         if ( prtflg(3) .ge. 1 .or. prtflg(4) .ge. 1
     x   .or. prtflg(5) .ge. 1 .or. prtflg(7) .ge. 1
     x   .or. diaflg(11) .ge. 1) then
            call tmat (ie,e,kmat,tr,ti,prtflg(3),nopen,a,c,ind,
     x                 nspn2,lrgl2,npty2,nchan,ichl,lchl,firste,
     x                 nume,nmt(3),ncol)
         end if
*
*  calculate cross sections
*
         if ( prtflg(4) .ge. 1 ) call crostn
     x                 (e,ntarg,nltarg,et,sw1,sw2,a,xslpi,xsn(1,ie),
     x                  diaflg(10),nxsn,sk,zion,ji,iind,
     x                  jind,tr,ti,nopen,ie,nl,nume,xsnl,il)
*
*  calculate partial collision strengths
*
         if (diaflg(11) .gt. 0 ) call parcol
     x                 (e,sw1,a,nchan,ichl,lchl,tr,ti,nopen)
*
*  calculate total collision strengths
*
         if (prtflg(5) .ge. 1) then
           nlt = ichl(nopen)
           call totcol (ntarg,nltarg,nlt,sw1,a,xslpi,xsn(1,ie),
     x                  diaflg(10),
     x                  nxsn,ion,ji,iind,jind,tr,ti,nopen,ie,
     x                  nl,nume,xsnl,il)
         end if
         if (prtflg(7) .ge. 1) then
            call smat1(nmt(7),ie,tr,ti,prtflg(7),nopen,nume,
     x                 nspn2,lrgl2,npty2,nltarg,ntarg,lchl,ncol,
     x                 icases,ncases,sce,x,ix)
         end if
      firste = .false.
  100 continue
      if (prtflg(2) .eq. 1 .or. prtflg(2) .eq. 3) then
         call etable (icon,quad,title,nume,phz,sce,nspn2,lrgl2,
     x                npty2,nmt(2))
      end if
      if (prtflg(4) .eq. 1  .or. prtflg(4) .eq. 3) then
         if (il .eq. 0) then
            call wrxsn(nmt(4),nxsn,nume,sce,xsn)
         else
            call wrxsnl(nmt(6),nxsn,il,nume,fl,sce,xsnl)
         end if
      end if
      if (prtflg(5) .eq. 1  .or. prtflg(5) .eq. 3) then
         if (il .eq. 0) then
            call wrxsn(nmt(5),nxsn,nume,sce,xsn)
         else
            call wrxsnl(nmt(6),nxsn,il,nume,fl,sce,xsnl)
         end if
      end if
      if (icases .eq. ncases .and. diaflg(13) .eq. 1) then
         if (maxtry .gt. 10) then
            sum = 0
            do 50 i = 10,maxtry
               sum = sum + rstats(i)
 50         continue
            rstats(10) = sum
         end if
         maxtri = min(10,maxtry)
         write (iout,1100)(rstats(i),i=1,maxtri)
      end if
      return
*
 1010 format (/,i4,' incident energy e =',f10.6,' scaled ryd ',
     x        f10.6,' scaled Bohr ',
     x        /,' etot =',f14.6,' scaled a.u.,  nopen =',i3,/,
     x        ' channel energies (scaled ryd) :',/,(2x,6f12.5))
 1020 format (' r-matrix at  r  =',f9.5,
     x        ' scaled Bohr  e =',f9.5,' scaled ryd :'/)
 1100 format (///1x,'trial radius statistics'/
     x           1x,'***********************'/
     x'    once twice   x3    x4    x5    x6    x7    x8    x9  over9',
     x        /, 1x,10i6 )
      end
      subroutine prsmat(nume,sce,ntarg,etarg,ltarg,iet,nmt7,
     x                  ncases,sline1,nelc1,nz1,ilast,x,ix)
      implicit integer (a-z)
      real*8 sce,e,etarg,e0,zero,atry,x
      character*1 num
      character*3 skk
      parameter ( zero=0.0d0, atry=2.0d0  )
*
      common /io/ iin,iout
      common /smpt/ ietc,nslp,ispst,illst,iprst,lln,llx,srmat,
     x              simat,sline,nelc,nz,ldim,nastm
      dimension etarg(ntarg),ltarg(ntarg),sce(nume),
     x          iet(ntarg),x(*),ix(*),num(0:9)
*
      data num/'0','1','2','3','4','5','6','7','8','9'/
*
      sline = sline1
      nelc = nelc1
      nz = nz1
      e = sce(nume)
      nastm = ntarg
      e0 = etarg(1)
      do 10 it = 1,ntarg
         do 20 ie = 1,nume
            if (sce(ie) .ge. (etarg(it) - e0)*atry) then
               iet(it)=ie
               go to 30
            end if
 20      continue
         nastm=it-1
         goto 40
 30      lat = ltarg(it)
         irecl=8*2*(((lat+1)*(lat+2))/2)*(nume-iet(it)+1)
         nunit=nmt7+it
         skk = 'S'//num(it/10)//num(it-10*(it/10))
         open(nunit,file=skk,status='unknown',access='direct',
     x        recl=irecl)
 10   continue
 40   ietc = ilast
      nslp = ietc + nastm
      ispst = nslp + nastm
      ndim = nastm * ncases
      illst = ispst + ndim
      iprst = illst + ndim
      lln = iprst + ndim
      llx = lln + ndim
      srmat = iadtwp(llx + ndim + 1)
      latm = ltarg(1)
      do 50 it = 2,nastm
         latm = max(latm,ltarg(it))
 50   continue
      ldim = ((latm+1)*(latm+2))/2
      sdim = ldim * nume * nastm
      simat = srmat + sdim
      last = simat + sdim
      ilast = wpadti(last + 1)
      do 60 it = 1,nastm
         ix(nslp + it -1) = 0
 60   continue
      return
      end

      subroutine prtdat(prtflg,rafin,raflw,jj,lrglx,lrglr,radwm,
     x                  lamax,iasy,abvthr,belthr,degeny,ewron,
     x                  beta,cuplim,deplim,ncrit)
*
* $Id: prtdat.f,v 2.1 94/09/29 18:39:37 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 rafin, raflw, radwm, abvthr, belthr, degeny, ewron,
     x     beta,cuplim,deplim
      logical jj
      dimension prtflg(7),nmt(7)
      common /io/ iin,iout
*
*     print input data
*
      write (iout,1000)
      if (prtflg(1) .ge. 1)
     x    write(iout,1001)
      if (prtflg(2) .ge. 1)
     x    write(iout,1002)
      if (prtflg(3) .ge. 1)
     x    write(iout,1003)
      if (prtflg(4) .ge. 1)
     x    write(iout,1004)
      if (prtflg(5) .ge. 1)
     x    write(iout,1005)
      if (prtflg(6) .ge. 1)
     x    write(iout,1006)
      if (prtflg(7) .ge. 1)
     x    write(iout,1007)
      if (jj) then
         write(iout,1020)
      else
         write(iout,1021)
      end if
      write(iout,1022) rafin, raflw
      write(iout,1023)lrglx
      write(iout,1024)lrglr,radwm
      write(iout,1025) lamax,iasy
      write(iout,1026) abvthr, belthr, degeny
      write(iout,1027) ewron
      write(iout,1028) beta, cuplim, deplim, ncrit
      return
*
 1000 format(/' namelist phzin data :')
 1001 format(' K matrices will be calculated')
 1002 format(' eigenphase sums will be calculated')
 1003 format(' T matrices will be calculated')
 1004 format(' X sections will be calculated')
 1005 format(' total collision strengths will be calculated')
 1006 format(' top up data will be calculated')
 1007 format(' S matrices will be calculated')
 1020 format(' JJ coupling is assumed.  L=2*J')
 1021 format(' LS coupling is assumed')
 1022 format(' estimated final radius, rafin =',f10.3,' scld Bohr',/,
     x       ' radius to which BBM propagator should be used,',
     x       ' raflw =',f10.3,' scld Bohr')
 1023 format(' Exchange is neglected for total L .gt. ',i3)
 1024 format(' Asymptotic Distorted Wave Method used at total L .gt.'
     x         ,i5,/,
     x       '    using an initial radius of',f6.3,
     x       'scld Bohr at which the R-matrix is calculated')
 1025 format(' no. of multipoles retained in potential, lamax =',i3,/,
     x       ' no. of terms in Gailitis expansion, iasy =',i3)
 1026 format(' scattering energies closer to a threshold',/
     x       ' than abvthr =',d16.8,' scld Ryd',/
     x       ' and  belthr =',d16.8,' scld Ryd are omitted.',/
     x       ' target levels closer than degeny =',d16.8,' scld Ryd',/
     x       ' are taken as degenerate in the Gailitis expansion.')
 1027 format(' Wronskian convergence parameter, ewron =',d16.8)
 1028 format(' Light-Walker sector size parameter beta =',d16.8,/
     x       ' closed channels with coupling less than cuplim =',d16.8,/
     x       ' are dropped if the channel energy .lt. -deplim =',d16.8,
     x       ' scld ryd',/
     x       ' except that ncrit =',i2,' of such channels are retained')
*
      end
      subroutine pxsnl (nl,fl,nume,nxsn,xsnl,en,
     x                  ixsn,iind,jind,ncol,nphys)
*
* $Id: pxsnl.f,v 2.1 94/09/29 18:40:14 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 xsnl,en,fl
      dimension xsnl(nxsn,nume,nl),en(nume),iind(nxsn),jind(nxsn),
     x          fl(nl)
      common /io/ iin,iout
*
*     pxsnl : print cross section table for each L
*
      do 10 il = 1,nl
         write (iout,1000)il,fl(il)
         call pxsntb (nume,nxsn,xsnl(1,1,il),en,
     x                ixsn,iind,jind,ncol,nphys)
 10   continue
      return
*
 1000 format (//' ',i4,' L=',f4.0, '(summed over spin and parity)')
      end
      subroutine pxsntb (nume,nxsn,xsn,en,ixsn,iind,jind,ncol,nphys)
*
* $Id: pxsntb.f,v 2.1 94/09/29 18:40:13 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 xsn,en
      dimension xsn(nxsn,nume),en(nume),iind(nxsn),jind(nxsn)
      common /io/ iin,iout
*
*     pxsntb : print cross section table
*
      if ( ixsn .ne. 0 ) then
         write(iout,*) ' nphys, ncol, nume:',nphys, ncol, nume
         do 20 i = 1, nphys*(nphys+1)/2, ncol
            ip = min(i+ncol-1,nphys*(nphys+1)/2)
            write (iout,1010) (jind(ij),iind(ij),ij = i,ip)
            do 10 ie = 1, nume
               write (iout,1020) ie,en(ie),(xsn(ntg,ie),ntg = i,ip)
   10       continue
   20    continue
      endif
      return
*
 1010 format (//,4x,'i',3x,'e(scld ryd )',6(i3,'-',i3,5x))
 1020 format (1x,i4,2x,d12.5,6d12.5)
      end
      subroutine rafmx (nchan,eth,enum,sce,cc,rafest,rafmax,raf,
     x                  rmatr,ion,lamax)
*
* $Id: rafmx.f,v 2.1 94/09/29 18:40:19 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 eth,sce,cc,eps,nu2,rafest,rafmax,raf,rafi,r,rmatr,
     x     zero,one,three,rion,rion2,atry,dfloat
      parameter ( atry=2.0d0,zero=0.0d0,one=1.0d0,three=3.0d0 )
      common /io/ iin,iout
      dimension eth(nchan),cc(nchan),sce(enum),raf(enum)
*
*     estimate the max value of rafin for the energy range
*
      rafmax = rafest
*     if there is no long range potential take rafin as default at
*      all energies
*     for positive ions
*      for open channels take rafin as 3 times the point of inflection
*      for closed channels midway between inner and outer
*     for neutrals
*      for open channels take rafin as .5*point of inflection
*      for closed take default
*     for negative ions take rafin as default
      if (ion .gt. 0 .and. lamax .gt. 0) then
         rion = dfloat (ion)
         rion2 = rion * rion
         do 40 ie = 1,enum
            rafi = rafest
            r = rafest
            do 50 i = 1,nchan
               eps = (sce(ie) - eth(i))/rion2
               if (eps .ge. zero) then
                  r = three / rion *
     x                  (sqrt(one + eps*cc(i)) - one) / eps
               else
                  nu2 = one / (-eps)
                  if ( nu2 .ge. cc(i) ) r = nu2 / rion
               end if
               rafi = max (rafi,r)
 50         continue
            rafmax = max (rafmax,rafi)
            raf(ie) = rafi
 40      continue
      else if (ion .eq. 0 .and. lamax .gt. 0) then
         do 60 ie = 1,enum
            rafi = rafest
            r = rafest
            do 70 i = 1,nchan
               eps = sce(ie) - eth(i)
               if (eps .gt. zero) then
                  r = sqrt(cc(i)/eps)
               end if
               rafi = max(rafi,r)
 70         continue
            rafmax = max(rafmax,rafi)
            raf(ie) = rafi
 60      continue
      else
         rafmax = rafest
         do 80 ie = 1,enum
            raf(ie) = rafest
 80      continue
      end if
c     write (iout,1000)rafmax
c1000 format (' estimated rafmax =',f10.4)
      return
*
      end
      subroutine rdcas (ls,ll,lpi,lrgs,lrgl,lrgpi,ncase,ncases)
*
* $Id: rdcas.f,v 2.1 94/09/29 18:39:59 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      common /io/ iin,iout
      dimension ls(ncase),ll(ncase),lpi(ncase),lrgs(ncase),
     x          lrgl(ncase),lrgpi(ncase)
*
*     copy list of cases
*
      ncases = 0
      do 10  i=1,ncase
         l = ll(i)
         if (l .lt. 0) go to 20
         lrgs(i) = ls(i)
         lrgl(i) = l
         lrgpi(i) = lpi(i)
         ncases = i
 10   continue
 20   return
      end
      subroutine rdlran(jj,llo,lup,ldel,nlran,lrglx,ncases,
     x                  ntarg,starg,mintsp,maxtsp,minlsp,maxlsp)
      implicit integer (a-z)
      logical jj
      common /io/ iin,iout
      dimension starg(ntarg),tsval(10),llo(nlran),lup(nlran),
     x          ldel(nlran)
*
*   calculate the range of target spin and large spin and
*   modify the ranges of large l to avoid overlap of the ranges
*   calculate the number of slpi cases
*
*   for j-j coupling large spin is taken as zero and l = 2 * j
*
      if (.not. jj) then
         ns = 1
         s1 = starg(1)
         tsval(1) = s1
         mintsp = s1
         maxtsp = s1
         do 10 i = 2,ntarg
            tsp = starg(i)
            do 20 j = 1,ns
               if (tsp .eq. tsval(j)) go to 10
 20         continue
            ns = ns + 1
            tsval(ns) = tsp
            mintsp = min (mintsp, tsp)
            maxtsp = max (maxtsp, tsp)
 10      continue
         ntspin = (maxtsp - mintsp + 2) / 2
         minlsp = mintsp - 1
         maxlsp = maxtsp + 1
         if (minlsp .eq. 0) then
            minlsp = 2
         end if
         nlspin = (maxlsp - minlsp + 2) / 2
      else
         mintsp = 0
         maxtsp = 0
         minlsp = 0
         maxlsp = 0
         ntspin = 1
         nlspin = 1
      end if
      ncases = 0
      nl = 0
      ltop = -1
      do 30 i = 1,nlran
         ld = max(ldel(i),1)
         if (llo(i) .le. ltop ) llo(i) = ltop + ld
         do 40 l = llo(i), lup(i), ld
            if (l .le. lrglx) then
               ncases = ncases + nlspin
            else
               ncases = ncases + ntspin
            end if
            ltop = l
 40      continue
 30   continue
      ncases = ncases * 2
      return
      end
      subroutine rmatd (nstat,nchan,eig,wmat,e0,cfbut,ibuttl,ichl,
     x                  etargs,rmatr,e,rmat,fkn,uk)
*
* $Id: rmatd.f,v 2.1 94/09/29 18:40:31 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 half,two,zero,cfbut,e,e0,eb,ediff,eig,etargs,etot,fkn,
     x       ra2,rmat,rmatr,sum,sumb,uk,wmat,but0,dfloat
      parameter ( zero=0.d0,half=0.5d0,two=2.0d0 )
      dimension rmat(*),eig(nstat),wmat(nstat,nchan),cfbut(3,nchan),
     x          etargs(*),ichl(nchan),fkn(*),uk(*)
*
*     rmatd : nchan-channel r-matrix for energy etot
*
      ra2 = rmatr * rmatr
      etot = e0 + half * e
      kk = 0
      do 70 i2 = 1, nchan
         do 60 i1 = 1, i2
            sum = zero
            do 50 i = 1, nstat
               ediff = eig(i) - etot
               sum = sum + wmat(i,i1) * wmat(i,i2) / ediff
   50       continue
            kk = kk + 1
            rmat(kk) = sum
C            print*,'kk, rmat(kk) =',kk,rmat(kk)
   60    continue
   70 continue
      if ( ibuttl .ne. 0 ) then
         do 80 i1 = 1, nchan
            eb = two * ( etot - etargs(ichl(i1)) )
            if ( cfbut(3,i1) .ge. dfloat(-10000) ) then
               sumb = eb * ( eb * cfbut(3,i1) + cfbut(2,i1) ) +
     x                cfbut(1,i1)
            else
               nbut = - int( cfbut(3,i1) ) / 10000
               sumb = cfbut(1,i1) * but0(nbut,fkn,uk,cfbut(2,i1) +
     x                ra2*eb)
            endif
                 kk = ( i1 * ( i1 + 1 ) ) / 2
            rmat(kk) = rmat(kk) + sumb
   80    continue
      endif
      return
*
      end
      subroutine rmatf(nchan,l2p,et,ion,rmatr,nfc,fc,fcp,gc,gcp,
     x                 iasy,iasx,xg,xc,xd,xcp,xdp,rmat)
      implicit integer (a-z)
      real*8 et,rmatr,rmat,fc,fcp,gc,gcp,dfloat,zero,one,skj,etaj,
     x       rho,flj,xg,xc,xd,xcp,xdp,zz,etj,u,uder
      parameter ( zero=0.d0, one=1.d0)
      common /io/ iin,iout
      dimension rmat(*),et(nchan),l2p(nchan),
     x          fc(nfc),gc(nfc),fcp(nfc),gcp(nfc),
     x          xg(0:iasy),xc(0:iasx),xd(0:iasx),xcp(0:iasx),
     x          xdp(0:iasx)
*
      zz = dfloat(ion)
      kk = 0
      do 70 i2 = 1, nchan
         do 60 i1 = 1, i2
            kk = kk + 1
            rmat(kk) = zero
   60    continue
   70 continue
      do 20 j = 1,nchan
         jj = (j * (j+1))/2
         lj = l2p(j)
         etj = et(j)
         if (etj .ge. zero) then
            skj = sqrt(abs(etj))
            etaj = -zz/skj
            rho = rmatr * skj
            flj = dfloat(lj)
            call coulfg(rho,etaj,flj,flj,fc,gc,fcp,gcp,1,0,ifail)
            if (ifail .ne. 0) then
*              write (iout,1000) ifail,lj,rho,etaj
*              write (iout,*) 'j,lj,skj,et(j)',j,lj,skj,et(j)
*   f = exp(- rho + abs(etaj) * log(rho+rho)) / sqrt(skj)
*   fp = - (skj - abs(etaj) / rmatr) * f
               rmat(jj) = -one / (skj - abs(etaj) )
            else
               rmat(jj) = fc(lj+1)/(rho * fcp(lj+1))
            end if
         else if ( ion .ne. 0 ) then
            call coul (lj,zz,etj,rmatr,u,uder,iasy,xg,xc,xd,xcp,xdp,
     x                 0,0)
c            skj = sqrt(abs(etj))
c            etaj = -zz/skj
c            rmat(jj) = -one / (skj - abs(etaj) )
             rmat(jj) = u/(rmatr * uder)
         else
            call decay (etj,lj,rmatr,u,uder)
            rmat(jj) = u/(rmatr * uder)
         endif
 20   continue
      return
 1000 format (' coulfg : ifail =',i4,' lj =',i4,' rho =',d16.8,
     x        ' etaj =',d16.8)
      end
      subroutine rmatf1(nchan,lchl,et,lmax,ion,iasy,rmatr,rmat,x)
      implicit integer (a-z)
      real*8 et,rmatr,rmat,x
      dimension rmat(*),et(nchan),lchl(nchan),x(*)
      common /glpt/ ga1,ga2,fc,fcp,gc,gcp,xg,xc,xd,xcp,xdp,fac1,fac2,
     x              fac3,fac4,iac1,deg,sma,smb
*
      nfc = lmax + 10
      iasx = iasy/2
      iasy1 = iasy-1
      call rmatf(nchan,lchl,et,ion,rmatr,nfc,x(fc),x(fcp),x(gc),x(gcp),
     x           iasy1,iasx,x(xg),x(xc),x(xd),x(xcp),x(xdp),rmat)
      return
      end
      subroutine secend(rafin,nrange,nranmx,asect,isect,nfinal,
     x                  lwsect,lwran)
*
* $Id: secend.f,v 2.1 94/09/29 18:39:49 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 rafin,asect,lwsect
      dimension asect(*),lwsect(lwran)
         do 21 i = 1,nranmx + 1
            if (rafin .le. asect(i)) then
               rafin = asect(i)
               nrange = i-1
               isect = 1
               go to 23
            end if
  21     continue
         nrange = nranmx
         do 22 i = 1,nfinal
            if (rafin .le. lwsect(i)) then
               rafin = lwsect(i)
               isect = i
               go to 23
            end if
 22      continue
 23      return
         end
      subroutine sleqn (a,lda,n,b,m,ipvt)
*
* $Id: sleqn.f,v 2.1 94/09/29 18:39:53 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 a,b,one
      parameter ( one=1.0d0 )
      common /io/ iin,iout
      dimension a(lda,n),b(lda,m),ipvt(n)
*
      call dgetrf (n,n,a,lda,ipvt,ier)
      if ( ier .ne. 0 ) then
         write (iout,1000) ier
         stop
      endif
      call dgetrs('n',n,m,a,lda,ipvt,b,lda,ier)
      if ( ier .ne. 0 ) then
         write (iout,1010) ier
         stop
      endif
      return
 1000 format (' sleqn : error return from dgetrf, info = ',i6)
 1010 format (' sleqn : error return from dgetrs, info = ',i6)
      end
      subroutine smat (ie,iet,trmat,timat,srmat,simat,ipr,nopen,
     x                 nume,nastm,ldim,nspn2,lrgl2,npty2,nltarg,ncol)
      implicit integer (a-z)
      real*8 trmat,timat,srmat,simat,zero,one,two
      parameter ( zero=0.0d0,one=1.d0,two=2.d0 )
      dimension trmat(*),timat(*),iet(*),nltarg(*),
     x          srmat(ldim,nastm,nume),simat(ldim,nastm,nume)
      common /io/ iin,iout
*
*     tmat : t-matrix computed from k-matrix
*            s = ( 1 + i * k ) / ( 1 - i * k ) ;  s = del + t
*            t = 2 * i * k / ( 1 - i * k )
*
      nf = 0
      do 60 it = 1, nastm
         if (ie.lt.iet(it)) return
         if (nltarg(it).eq.0) goto 60
         ni = nf + 1
         nf = nf + nltarg(it)
            k = 0
            do 30 nc = ni, nf
               k = k + 1
               nmc = nopen * ( nc - 1 ) + nc
               srmat(k,it,ie) = trmat(nmc) + one
               simat(k,it,ie) = timat(nmc)
               do 20 mc = nc+1, nf
                  k = k + 1
                  nmc = nopen * ( nc - 1 ) + mc
                  srmat(k,it,ie) = trmat(nmc)
                  simat(k,it,ie) = timat(nmc)
   20          continue
   30       continue
   60 continue

*     if ipr = 2,3  print s-matrix
      if ( ipr .gt. 1 ) then
         do 70 it = 1,nastm
            nt = nltarg(it)
            if (nt .gt. 0) then
               write (iout,1020)
               call wrttmt (srmat(1,it,ie),nt,ncol)
               write (iout,1030)
               call wrttmt (simat(1,it,ie),nt,ncol)
            end if
 70      continue
      endif
      return
*
 1020 format (/,' real part of s-matrix :')
 1030 format (/,' imag part of s-matrix :')
      end
      subroutine smat1(nmt7,ie,tr,ti,ipr,nopen,nume,nspn2,lrgl2,
     x                 npty2,nltarg,ntarg,lchl,ncol,icase,ncases,
     x                 sce,x,ix)
      implicit integer (a-z)
      real*8 tr,ti,sce,x
      common /smpt/ iet,nslp,ispst,illst,iprst,lln,llx,srmat,
     x              simat,sline,nelc,nz,ldim,nastm
      dimension tr(*),ti(*),nltarg(ntarg),x(*),ix(*),lchl(nopen),
     x          sce(nume)
*
      call smat (ie,ix(iet),tr,ti,x(srmat),x(simat),ipr,
     x           nopen,nume,nastm,ldim,nspn2,lrgl2,npty2,nltarg,
     x           ncol)
      if ((ipr .eq. 1  .or. ipr .eq. 3)  .and. ie .eq. nume) then
         call wrsxx(nmt7,ix(iet),x(srmat),x(simat),nltarg,lchl,
     x              ix(nslp),ix(lln),ix(llx),ix(ispst),ix(illst),
     x              ix(iprst),ldim,nastm,ncases,nume,nspn2,lrgl2,
     x              npty2)
         if (icase .eq. ncases ) then
             call wrs00(nmt7,nelc,nz,ntarg,ix(nslp),ix(ispst),
     x                 ix(illst),ix(iprst),ix(iet),ix(lln),ix(llx),
     x                 sce,x(sline),nastm,ncases,nume)
         end if
      end if
      return
      end
      subroutine split (nchan,ichl,ntarg,starg,ispch,ichsp,nschan,
     x                  nspin,ns)
*
* $Id: split.f,v 2.1 94/09/29 18:40:21 vmb Exp Locker: vmb $
*
      implicit integer (  a - z )
      dimension ichl(nchan),starg(ntarg),ispch(nchan),ichsp(nchan),
     x          nschan(2),nspin(2)
*
*     split the channels into two possible spins
*     array ispch relates the new channel order to the old
*
      sl = starg (ichl(1))
      nch1 = 1
      ispch(1) = 1
      ichsp(1) = 1
      do 10  ich = 2,nchan
         if (sl .eq. starg(ichl(ich))) nch1 = nch1 + 1
 10   continue
      nschan(1) = nch1
      nspin(1)  = sl
      if (nch1 .eq. nchan) then
         ns = 1
         nschan(2) = 0
         do 20  ich = 2,nchan
            ispch(ich) = ich
            ichsp(ich) = ich
 20      continue
      else
         nch2 = nch1
         nch1 = 1
         ns = 1
         do 30  ich = 2,nchan
            sn = starg(ichl(ich))
            if (sn .ne. sl) ns = -ns
            if (ns .gt. 0) then
               nch1 = nch1 + 1
               ispch(nch1) = ich
               ichsp(ich) = nch1
            else
               nch2 = nch2 + 1
               ispch(nch2) = ich
               ichsp(ich) = nch2
            end if
            sl = sn
 30      continue
      ns = 2
      nspin(2) = starg(ichl(nch1 + 1))
      nschan(2) = nchan - nch1
      end if
      return
      end
      subroutine splitr(nchan,rmat,sprmat,ispch)
*
* $Id: splitr.f,v 2.1 94/09/29 18:40:23 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 rmat,sprmat
      dimension rmat(*),sprmat(nchan,nchan),ispch(nchan)
*
*     split rmat into spins by rearranging channels
*
      kk = 0
      do 10 i2 = 1,nchan
         ich2 = ispch(i2)
         do 10 i1 = 1,i2
            kk = kk + 1
            ich1 = ispch(i1)
            sprmat(ich1,ich2) = rmat(kk)
            sprmat(ich2,ich1) = rmat(kk)
 10   continue
      kk = 0
      do 20 i2 = 1,nchan
         do 20 i1 = 1,i2
            kk = kk + 1
            rmat(kk) = sprmat(i1,i2)
 20   continue
      return
      end
      subroutine statw (lrgl,nspn,ntarg,sw1,sw2,ltarg,starg)
      implicit integer (a-z)
      real*8 sw1,sw2,dfloat
      dimension sw2(ntarg),ltarg(ntarg),starg(ntarg)
*
*  nspn < 0 corresponds to no-exchange in inner region
*  nspn = 0 corresponds to Breit Pauli in inner region lrgl=2 * J
*           ltarg(i)=2*jtarg(i)
*
      if (nspn .gt. 0) then
         sw1 = ( 2 * lrgl + 1 ) * nspn
      else if (nspn .lt. 0) then
         sw1 = -2 * ( 2 * lrgl + 1 ) * nspn
      else
         sw1 = lrgl + 1
      end if
      if (nspn .ne. 0) then
         do 10 i = 1, ntarg
            lf = 2 * ltarg(i) + 1
            sw2(i) = dfloat( lf*starg(i) )
   10    continue
      else
         do 20 i = 1,ntarg
            sw2(i) = dfloat(ltarg(i) + 1)
   20    continue
      end if
      return
      end
      subroutine tmat (ie,e,akmat,trmat,timat,ipr,nopen,a,b,ind,
     x                 nspn2,lrgl2,npty2,nchan,ichl,lchl,firste,
     x                 nume,nmt3,ncol)
*
* $Id: tmat.f,v 2.1 94/09/29 18:40:16 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 e,akmat,trmat,timat,a,b,x,zero,one,two
      logical firste
      parameter ( zero=0.0d0,one=1.d0,two=2.d0 )
      dimension akmat(*),trmat(*),timat(*),a(*),b(*),ind(*),
     x          ichl(nchan),lchl(nchan)
      common /io/ iin,iout
*
*     tmat : t-matrix computed from k-matrix
*            s = ( 1 + i * k ) / ( 1 - i * k ) ;  s = del + t
*            t = 2 * i * k / ( 1 - i * k )
*
      nop2 = nopen * nopen
      if ( nopen .eq. 1 ) then
         x = akmat(1)
         a(1) = x * x
         b(1) = one / ( one + x * x )
         timat(1) = x * b(1)
         trmat(1) = a(1) * b(1)
      else
         ix = 1
         do 4 i = 1, nopen
            do 2 j = 1, i
               b((i-1)*nopen+j) = akmat(ix)
               b((j-1)*nopen+i) = akmat(ix)
               ix = ix + 1
    2       continue
    4    continue
         call dcopy (nop2,b,1,timat,1)
         call dgemm ('n','n',nopen,nopen,nopen,one,b,nopen,b,nopen,
     x                zero,a,nopen)
         call dcopy (nop2,a,1,trmat,1)
         do 10 i = 1, nopen
            ii = nopen * ( i - 1 ) + i
            a(ii) = a(ii) + one
   10    continue
         call dgetrf (nopen,nopen,a,nopen,ind,ier)
         if ( ier .ne. 0 ) then
            write (iout,1000) ier
            stop
         endif
         call dgetrs('n',nopen,nopen,a,nopen,ind,trmat,nopen,ier)
         call dgetrs('n',nopen,nopen,a,nopen,ind,timat,nopen,ier)
         if ( ier .ne. 0 ) then
            write (iout,1010) ier
            stop
         endif
      endif
      do 30 i = 1, nop2
         trmat(i) = - two * trmat(i)
         timat(i) = two * timat(i)
   30 continue
*     if ipr = 1,3 write t-matrix to nmt3 ; = 2,3  print t-matrix
      if ( ipr .gt. 1 ) then
         write (iout,1020)
         call wrtmat (trmat,nopen,nopen,nopen,ncol)
         write (iout,1030)
         call wrtmat (timat,nopen,nopen,nopen,ncol)
      endif
      if ( ipr .eq. 1 .or. ipr .eq. 3  ) then
         if (firste) then
            numen = nume-ie+1
            write (nmt3) nspn2,lrgl2,npty2,numen,nchan,
     x                  (ichl(ich),lchl(ich),ich=1,nchan)
         end if
*  write the lower triangle of the t-matrix to file
         k = 0
         l = -nopen
         do 40 i = 1,nopen
            l = l+nopen
            m = l
            do 50 j = 1,i
               k = k+1
               m = m+1
               a(k) = trmat(m)
               b(k) = timat(m)
 50         continue
 40      continue
         write (nmt3) e,nopen,k,(a(i),b(i),i=1,k)
*         write (nmt3) e,nopen,nop2,(trmat(i),timat(i),i=1,nop2)
      endif
      return
*
 1000 format (' tmat : error return from dgetrf, info = ',i6)
 1010 format (' tmat : error return from dgetrs, info = ',i6)
 1020 format (/,' real part of t-matrix :')
 1030 format (/,' imag part of t-matrix :')
      end
      subroutine topdat(ntarg,etarg,abvthr,belthr,ebig,sce2,
     x                  enum2,ncases,lrgl,lrgle1,nl,fl,lstart,lfin,
     x                  nxsn,xsnl,tlast,x,ix)
*
* $Id: topdat.f,v 2.1 94/09/29 18:41:01 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 esc2,x,etarg,abvthr,belthr,ebig,zero
      common /io/ iin, iout
      dimension ne2(10),esc2(20),x(*),ix(*),etarg(ntarg),
     x          lrgl(ncases)
      namelist /topin/ lrgle,ne2,esc2
      data zero/0.0d0/,ne2/10*0/
      read (iin,topin,end=10,err=10)
      go to 11
 10   write (iout,topin)
      stop
 11   lrgle1 = lrgle
      write (iout,1000) lrgle
      nume2 = necnt(ne2)
      sce2 = tlast
      tlast = sce2 + nume2 + 1
      scest = tlast
      norder = wpadti(scest + nume2 + 1)
*
*  store the energy mesh after pruning out energies too close to
*  threshold
*
      call evals (ntarg,nume2,ne2,esc2,etarg,abvthr,belthr,
     x            x(sce2),enum2,ebig,x(scest),ix(norder))
      if (enum2 .lt. nume2) then
         write (iout,2010) enum2
      end if
      call lvals (ncases,lrgl,lrgle,nl,x(tlast))
      fl =tlast
      tlast =fl + nl
*
      call llims(ncases,lrgl,lstart,lfin)
*
      xsnl = tlast
      len = nl * enum2 * nxsn
      tlast = xsnl + len
*
*  zeroise the xsnl array
*
      do 20 i = 1,len
         x(xsnl + i - 1) = zero
 20   continue
      return
 1000 format (' start top up procedure at L=lrgle = ',i4)
 2010 format (' coarse energy mesh for top-up',/,
     x        ' number of discrete energies away from thresholds = ',i4)
      end
      subroutine totcol (ntarg,nltarg,nlt,sw1,a,omslpi,omega,
     x                   idbug,nxsn,ion,ji,iind,jind,tr,ti,nopen,ien,
     x                   nl,nume2,oml,il)
*
* $Id: totcol.f,v 2.1 94/09/29 18:40:26 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 a,omslpi,omega,oml,tr,ti,zero,half,c1,sw1,cros
      parameter ( zero=0.d0,half=0.5d0 )
      dimension nltarg(ntarg),tr(*),ti(*),ji(*),a(*),
     x          omslpi(nxsn),omega(nxsn),iind(nxsn),jind(nxsn),
     x          oml(nxsn,nume2,nl)
      common /io/ iin,iout
*
      nop2 = nopen * nopen
      do 12 i = 1, nop2
         a(i) = tr(i)**2 + ti(i)**2
   12 continue
      ix = 0
      c1 = half * sw1
      nf = 0
      do 60 i = 1, nlt
         ni = nf + 1
         nf = nf + nltarg(i)
         mf = 0
         jup = i
         if (ion .gt. 0) jup = i-1
         do 40 j = 1, jup
            ix = ix + 1
            mi = mf + 1
            mf = mf + nltarg(j)
            cros = zero
            do 30 nc = ni, nf
               do 20 mc = mi, mf
                  nmc = nopen * ( nc - 1 ) + mc
                  cros = cros + a(nmc)
   20          continue
   30       continue
            omslpi(ix) = c1 * cros
   40    continue
         ji(ien) = ix
   60 continue
*     debug print of collision strength
      if ( idbug .gt. 0 .and. ix .gt. 0 ) then
         write (iout,1000)
         nloop = ix/6 + 1
         n2 = 0
         do 90 ni = 1, nloop
            n1 = n2 + 1
            n2 = min(n1 + 5,ix)
            write (iout,1010) (jind(ij),iind(ij),ij = n1,n2)
            write (iout,1020) (omslpi(ij),ij = n1,n2)
   90    continue
      endif
      if (il .le. 0) then
         do 100 i = 1,ix
            omega(i) = omega(i) + omslpi(i)
 100     continue
      else
            do 101 i = 1,ix
               oml(i,ien,il) = oml(i,ien,il) + omslpi(i)
 101        continue
      end if
      return
*
 1000 format (/,10x,'collision strengths'/10x,'*******************')
 1010 format (1x,6(i3,'-',i3,5x))
 1020 format (1x,6d12.5)
      end

      double precision function wabs1(a,b,c,d,f)
*
*  calculates absolute value of racah coefficient
*  w(a,b,c,d,1,f)
*
      implicit integer(a-z)
      real*8  dfloat,half
      common /io/ iin, iout
      data half/0.5d0/
*
      x=0
      y=0
      if (f .le. min(b+d,a+c)
     x    .and. f .ge. max(abs(b-d),abs(a-c))) then
         if (b .eq. a-1) then
            if (d .eq. c-1) then
               x=(b+d+f+2)*(b+d+f+3)*(b+d-f+1)*(b+d-f+2)
               y=(b+1)*(2*b+1)*(2*b+3)*(d+1)*(2*d+1)*(2*d+3)
            else if (d .eq. c) then
               x=(b+d+f+2)*(-b+d+f)*(b-d+f+1)*(b+d-f+1)
               y=(b+1)*(2*b+1)*(2*b+3)*d*(d+1)*(2*d+1)
            else if (d .eq. c+1) then
               x=(-b+d+f-1)*(-b+d+f)*(b-d+f+1)*(b-d+f+2)
               y=(b+1)*(2*b+1)*(2*b+3)*d*(2*d-1)*(2*d+1)
            end if
         else if (b .eq. a) then
            if (d .eq. c-1) then
               x=(b+d+f+2)*(-b+d+f+1)*(b-d+f)*(b+d-f+1)
               y=b*(b+1)*(2*b+1)*(d+1)*(2*d+1)*(2*d+3)
            else if (d .eq. c) then
               x=b*(b+1)+d*(d+1)-f*(f+1)
               y=b*(b+1)*(2*b+1)*d*(d+1)*(2*d+1)
               x=x*x
            else if (d .eq. c+1) then
               x=(b+d+f+1)*(-b+d+f)*(b-d+f+1)*(b+d-f)
               y=b*(b+1)*(2*b+1)*d*(2*d-1)*(2*d+1)
            end if
         else if (b .eq. a+1) then
            if (d .eq. c-1) then
               x=(-b+d+f+1)*(-b+d+f+2)*(b-d+f-1)*(b-d+f)
               y=b*(2*b-1)*(2*b+1)*(d+1)*(2*d+1)*(2*d+3)
            else if (d .eq. c) then
               x=(b+d+f+1)*(-b+d+f+1)*(b-d+f)*(b+d-f)
               y=b*(2*b-1)*(2*b+1)*d*(d+1)*(2*d+1)
            else if (d .eq. c+1) then
               x=(b+d+f)*(b+d+f+1)*(b+d-f-1)*(b+d-f)
               y=b*(2*b-1)*(2*b+1)*d*(2*d-1)*(2*d+1)
            end if
         end if
      end if
      if (x .ne. 0  .and. y .ne. 0) then
         wabs1=half*sqrt(dfloat(x)/dfloat(y))
      else
         write(iout,2000)a,b,c,d,f
         stop
      end if
 2000 format(' stop in wabs1 while calculating',/,
     x      ' W(',i2,',',i2,',',i2,',',i2,'; 1,',i2,')')
*
      return
      end
      subroutine wronsk (nchan,nopen,f,fp,idiag1,idiag2,eps,fail)
*
* $Id: wronsk.f,v 2.1 94/09/29 18:39:54 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      logical fail
      real*8 f,fp,eps,zero,one,sum,tsum
      parameter ( zero=0.0d0,one=1.0d0 )
      common /io/ iin,iout
      dimension f(nchan,nchan,2),fp(nchan,nchan,2)
*
      fail = .false.
      if ( idiag2 .ne. 0 ) then
         write (iout,1010)
         do 20 i = 1, nchan
            write (iout,1000) (f(i,j,1),j = 1,nchan)
   20    continue
         write (iout,1020)
         do 30 i = 1, nchan
            write (iout,1000) (fp(i,j,1),j = 1,nchan)
   30    continue
         write (iout,1030)
         do 40 i = 1, nchan
            write (iout,1000) (f(i,j,2),j = 1,nchan)
   40    continue
         write (iout,1040)
         do  50  i = 1, nchan
            write (iout,1000) ( fp(i,j,2),j = 1,nchan)
   50    continue
      endif
*
*     check multichannel wronskian relations
*
      do 90 i0 = 1, 2
         if ( idiag2 .ne. 0 ) write (iout,1050) i0
         do 90 j1 = 1, nchan
            do 90 j2 = 1, nchan
               sum = zero
               tsum = zero
               do 70 i = 1, nchan
                  sum = sum + fp(i,j1,1) * f(i,j2,i0) - f(i,j1,1) *
     x                  fp(i,j2,i0)
   70          continue
               if ( i0 .eq. 2 .and. j1 .eq. j2 .and. j1 .le. nopen )
     x              tsum = one
               if ( abs(sum-tsum) .gt. eps ) then
                  fail = .true.
                  if (idiag1 .ne. 0) then
                     write (iout,1070) i0,j1,j2,sum
                  end if
               endif
               if ( idiag2 .ne. 0 ) write (iout,1060) j1,j2,sum
   90 continue
      return
*
 1000 format (1x,4d30.15)
 1010 format ('0regular functions')
 1020 format ('0regular derivatives')
 1030 format ('0irregular functions')
 1040 format ('0irregular derivatives')
 1050 format (/,' multichannel wronskian relations :',/,' i0 =',i3,/)
 1060 format (' i1 =',i3,5x,' i2 =',i3,5x,'wronskian =',d30.17)
 1070 format (' i0 =',i3,4x,' i1 =',i3,4x,' i2 =',i3,4x,'wronskian =',
     x        d26.17)
      end
      subroutine wrs00(nfts,nelc,nz,ntarg,nslp,ispst,illst,iprst,iet,
     x                 lln,llx,sce,sline,nastm,ncases,nume)
*
      implicit integer (a-z)
      real*8 sce,sline,zero
      parameter(zero = 0.0d0)
      dimension iet(*),nslp(*),sce(*),sline(*),
     x          lln(nastm,ncases),llx(nastm,ncases),
     x          ispst(nastm,ncases),illst(nastm,ncases),
     x          iprst(nastm,ncases)
*
      write(nfts) nz,nelc,ntarg,nastm
      write(3,*) nz,nelc,ntarg,nastm
      write(nfts) nume,(sce(ie),ie=1,nume),(zero,ie=1,nume)
      write(3,*) nume,(sce(ie),ie=1,nume),(zero,ie=1,nume)
      write(nfts) (sline(k),k=1,((ntarg-1)*ntarg)/2)
      write(3,*) (sline(k),k=1,((ntarg-1)*ntarg)/2)
      do 10 it=1,nastm
         write(nfts) iet(it),nslp(it),(ispst(it,n),illst(it,n),
     x               iprst(it,n),lln(it,n),llx(it,n),n=1,nslp(it))
         write(3,*) iet(it),nslp(it),(ispst(it,n),illst(it,n),
     x               iprst(it,n),lln(it,n),llx(it,n),n=1,nslp(it))
   10 continue
*
      return
      end




      subroutine wrslin(ntarg,nxsn,sline)
      implicit integer(a - z)
      real*8 sline,zero
      common /io/ iin, iout
      dimension sline(nxsn)
      data zero/0.0d0/
*  print sline
      write(iout,2000)
      do 10 it=2,ntarg
         do 20 jt=1,it-1
            k=((it-1)*(it-2))/2+jt
            if(sline(k).gt.zero) write(iout,2010)jt,it,sline(k)
   20    continue
   10 continue
 2000 format(//,6x,' Line Strengths',/,6x,' **************',/
     x       ' lower  upper'/' state  state',/ )
 2010 format(2x,i3,4x,i3,4x,f10.4)
      return
      end

      subroutine wrsxx(nfts,iet,srmat,simat,nltarg,lchl,nslp,lln,llx,
     x                  ispst,illst,iprst,ldim,nastm,ncases,nume,
     x                  nspn2,lrgl2,npty2)
*
      implicit integer (a-z)
      real*8 srmat,simat
      dimension iet(*),srmat(ldim,nastm,nume),simat(ldim,nastm,nume),
     x          nltarg(*),nslp(*),lln(nastm,ncases),llx(nastm,ncases),
     x          ispst(nastm,ncases),illst(nastm,ncases),
     x          iprst(nastm,ncases),lchl(*)
*
      nf = 0
      do 10 it = 1, nastm
         if (nltarg(it).eq.0) goto 10
         ni = nf + 1
         nf = nf + nltarg(it)
         nslp(it) = nslp(it) + 1
         n = nslp(it)
         ln = lchl(ni)
         lx = lchl(nf)
         lln(it,n) = ln
         llx(it,n) = lx
         ispst(it,n) = nspn2
         illst(it,n) = lrgl2
         iprst(it,n) = npty2
         km = ((2 + lx - ln) * (4 + lx - ln)) / 8
         nunit = nfts + it
         write(nunit,rec=n) ((srmat(k,it,ie),simat(k,it,ie),k=1,km),
     x                      ie=iet(it),nume)
         write(3,*)' nunit,n',nunit,n,((srmat(k,it,ie),simat(k,it,ie),
     x   k=1,km),ie=iet(it),nume)
   10 continue
*
      return
      end




      subroutine wrtmat (a,n,m,l,ncol)
*
* $Id: wrtmat.f,v 2.1 94/09/29 18:39:56 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 a
      common /io/ iin,iout
      dimension a(l,*)
*
*     wrtmat : prints an n*m matrix stored in an l*l array
*
      ntim = m / ncol
      nf = 0
      if ( ntim .gt. 0 ) then
         do 20 i = 1, ntim
            ni = nf + 1
            nf = nf + ncol
            do 10 j = 1, n
               write (iout,1000) (a(j,k),k = ni,nf)
   10       continue
            write (iout,1010)
   20    continue
      endif
      ni = nf + 1
      if ( ni .le. m ) then
         do 40 j = 1, n
            write (iout,1000) (a(j,k),k = ni,m)
   40    continue
      endif
      return
*
 1000 format (7d13.5)
 1010 format (/)
      end
      subroutine wrttmt (a,n,ncol)
*
* $Id: wrttmt.f,v 2.1 94/09/29 18:39:58 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 a
      common /io/ iin,iout
      dimension a(*)
*
*     wrtmat : prints an n*n matrix stored in a lower triangular array
*
      nf = 0
      do 40 k = 1, n
         ntim = k / ncol
         if ( ntim .gt. 0 ) then
            do 20 i = 1, ntim
               ni = nf + 1
               nf = nf + ncol
               write (iout,1000) (a(kk),kk = ni,nf)
   20       continue
         endif
         ires = mod(k,ncol)
         if ( ires .ge. 1 ) then
            ni = nf + 1
            nf = nf + ires
            write (iout,1000) (a(kk),kk = ni,nf)
         endif
         write (iout,1010)
   40 continue
      return
*
 1000 format (7d13.5)
 1010 format (/)
      end
      subroutine wrxsn(nmt,nxsn,nume1,sce1,xsn)
*
* $Id: wrxsn.f,v 2.1 94/09/29 18:40:58 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 sce1,xsn
*
*  write file of total cross sections or total collision strengths
*
*  nume1 = number of energies in fine energy mesh
*  l1max = highest L for which cross sections will have been written
*
      dimension sce1(nume1),xsn(nxsn,nume1)
*
cz
      rewind(nmt)

      write(nmt)nxsn,nume1
      write(nmt)(sce1(ie1),ie1=1,nume1)
      write(nmt)((xsn(ix,ie1),ix=1,nxsn),ie1=1,nume1)

cz
cz      backspace(nmt)
cz      backspace(nmt)
cz      backspace(nmt)
cz
*
      return
      end

      subroutine wrxsnl(nmt,nxsn,nl,nume2,fl,sce2,xsnl)
*
* $Id: wrxsnl.f,v 2.1 94/09/29 18:41:00 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 fl,sce2,xsnl
*
*  write the top up file for  cross-sections or collision strengths
*
*  nl = L dimension of xsnl = max no. of L values attempted at the coarse
*         energy mesh
*  nume2 = number of energies in coarse energy mesh
*
      dimension fl(nl),sce2(nume2),xsnl(nxsn,nume2,nl)
*
cz
      rewind(nmt)
cz
      write(nmt) nxsn,nume2,nl
      write(nmt) (sce2(ie2),ie2 = 1,nume2)
      write(nmt) (fl(il),il = 1,nl)
      write(nmt) (((xsnl(ix,ie2,il), ix = 1,nxsn),ie2 = 1,nume2),
     x           il = 1,nl)
cz
cz      backspace (nmt)
cz      backspace (nmt)
cz      backspace (nmt)
cz      backspace (nmt)
*
      return
      end
      subroutine wtcas (lrgs,lrgl,lrgpi,ncases)
*
* $Id: wtcas.f,v 2.1 94/09/29 18:40:05 vmb Exp Locker: vmb $
*
*
*     print list of cases whose data has not been found
*
      dimension lrgs(ncases),lrgl(ncases),lrgpi(ncases)
*
      common /io/ iin,iout
*
      do 10  i=1,ncases
         if (lrgl(i) .ge. 0) then
            write (iout,1010) lrgs(i),lrgl(i),lrgpi(i)
         end if
 10   continue
      return
 1010 format ('NO DATA FOUND FOR 2S + 1 =',i3,'  L =',i4,'  PI =',i2)
*
      end
      subroutine xsnind(ntarg,nxsn,ion,iind,jind)
*
* $Id: xsnind.f,v 2.1 94/09/29 18:40:57 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      dimension iind(nxsn),jind(nxsn)
*
*  calculate initial and final target index for each x-section
*
      ij = 0
      do 11 i = 1,ntarg
         jup = i
         if (ion .gt. 0) jup = i-1
         do 11 j = 1,jup
            ij = ij + 1
            iind(ij) = i
            jind(ij) = j
 11   continue
      return
      end
      subroutine ampltd (ra1,ra2,hmlt,ampa1,ampa2,eigen,eigo,vecto,
     x                   wrk,lbug,nchan,nhsize,nhd,nleg)
*
* $Id: ampltd.f,v 2.1 94/09/29 18:44:27 vmb Exp Locker: vmb $
*
      implicit real*8 (a-h,o-z)
      parameter ( zero=0.0d0,one=1.0d0 )
      dimension ampa1(nchan,nhd),ampa2(nchan,nhd),eigen(nhd),wrk(3*nhd),
     x          hmlt(nhsize),eigo(nhd),vecto(nhd,nhd)
      common /io/ iin,iout
*
*      ampltd : diagonalizes hamiltonian on subrange, stores eigenvalues
*               uses eigenvectors to compute surface amplitudes
*
*               ra1,ra2 = radii of end points of current subrange
*               hmlt    = hamiltonian ( lower triangle )
*               lbug    = debug switch
*               eigo, vecto, wrk : work space
*               ampa1,ampa2 = surface amplitudes
*               eigen       = eigenvalues of hamiltonian
*
      if ( lbug .gt. 0 ) write (iout,1000)
      ra12 = ra2 - ra1
      call DSPEV ( 'V', 'U', nhd, hmlt, eigo, vecto, nhd, wrk, ierr )
      if ( ierr .ne. 0 ) then
         write (iout,1040) ierr
         stop
      endif
      if ( lbug .gt. 0 ) then
         write (iout,1010) (eigo(neig),neig = 1,nhd)
         if ( lbug .ge. 2 ) then
            do 10 k = 1, nhd
               write (iout,1020) k
               write (iout,1010) (vecto(nvec,k),nvec = 1,nhd)
 10         continue
         endif
      endif
      do 50 k = 1, nhd
         eigen(k) = eigo(k)
         do 20 i = 1, nchan
            ampa1(i,k) = zero
            ampa2(i,k) = zero
   20    continue
         ik = 0
         sign = - one
*      sum over index on legendre polynomials
         do 40 j = 1, nleg
            rootf = sqrt( ( dfloat(2*j -1) ) / ra12 )
            sign = - sign
            roots = sign * rootf
            do 30 i = 1, nchan
               ik = ik + 1
               ampa1(i,k) = ampa1(i,k) + vecto(ik,k) * roots
               ampa2(i,k) = ampa2(i,k) + vecto(ik,k) * rootf
   30       continue
   40    continue
*      output of surface amplitudes
         if ( lbug .ge. 2 ) then
            write (iout,1030) ra1
            write (iout,1010) (ampa1(i,k),i = 1,nchan)
            write (iout,1030) ra2
            write (iout,1010) (ampa2(i,k),i = 1,nchan)
         endif
   50 continue
      return
*
 1000 format (1h1,/,52x,' subroutine ampltd',/,53x,17(' '),
     x        /,5x,'eigenvalues in rydbergs',/)
 1010 format (1x,8f14.7)
 1020 format (/,5x,'eigenvector',i5,/)
 1030 format (/,5x,'amplitude at =',f10.5,3x,'a.u.',/)
 1040 format (' ampltd : error return from dspev, ierr =',i4)
      end
      subroutine jnvert (ndim,nf,b,c,d,a,dd,nbb,nff)
*
* $Id: jnvert.f,v 2.1 94/09/29 18:44:30 vmb Exp Locker: vmb $
*
      implicit real*8 (a-h,o-z)
      parameter ( one=1.0d0 )
      dimension b(nbb),a(nff),d(nff),c(ndim,nf),dd(ndim)
*
*     jnvert : evaluates d = a - c(transpose) * b(inverse) * c
*              cholesky type factorization,
*              nesbet, j. comp. phys 8 (1971) 483
*              lower triangles of a, b, d stored row-wise
*
*     replace b with triangular factorization, nesbet eqns 12, 13
*
      kij = 0
      do 70 i = 1, ndim
         mm = i - 1
         ii = kij
         if ( i .eq. 1 ) go to 40
         jj = 0
         do 30 j = 1, mm
            q = b(kij+j)
            if ( j .eq. 1 ) go to 20
            nn = j - 1
            do 10 nu = 1, nn
               q = q - b(ii+nu) * dd(nu) * b(jj+nu)
   10       continue
   20       jj = jj + j
            b(kij+j) = q / ( b(jj) * dd(j) )
   30    continue
   40    kij = kij + i
         q = b(kij)
         if ( i .eq. 1 ) go to 60
         do 50 nu = 1, mm
            q = q - b(ii+nu) * dd(nu) * b(ii+nu)
   50    continue
   60    dd(i) = sign(one,q)
         b(kij) = sqrt(abs(q))
   70 continue
*
*     replace c by matrix defined by nesbet eqn 14
*
      do 90 ich = 1, nf
         c(1,ich) = c(1,ich) / b(1)
         if ( ndim .eq. 1 ) go to 100
         kij = 1
         do 90 i = 2, ndim
            q = c(i,ich)
            i1 = i - 1
            do 80 k = 1, i1
               q = q - c(k,ich) * b(kij+k)
   80       continue
            kij = kij + i
            c(i,ich) =  q / b(kij)
   90 continue
*
*     set up final matrix in d using nesbet eqn 11
*
  100 kka = 0
      do 130 ich = 1, nf
         do 120 jch = 1, ich
            q = a(kka+jch)
            do 110 k = 1, ndim
               q = q - c(k,ich) * dd(k) * c(k,jch)
  110       continue
            d(kka+jch) = q
  120    continue
         kka = kka + ich
  130 continue
      return
      end
      subroutine legndr (pl,xi,nleg,ixmax)
*
* $Id: legndr.f,v 2.1 94/09/29 18:44:32 vmb Exp Locker: vmb $
*
      implicit real*8 (a-h,o-z)
      parameter ( half=0.5d0,one=1.0d0,onep5=1.5d0,two=2.d0 )
      dimension pl(nleg,ixmax),xi(ixmax)
*
*     legndr : normalized legendre polynomials at quadrature points
*              nleg  = maximum number of legendre polynomials in basis
*              ixmax = number of abscissae
*              xi    = quadrature abscissae
*              pl    = polynomials
*
      nleg1 = nleg-1
      do 20 i = 1, ixmax
         x = xi(i)
         pl(1,i) = sqrt(half)
         pl(2,i) = sqrt(onep5) * x
*
*      generate polynomials to order nleg-1 by recursion
*
         do 10 k = 2, nleg1
            k1 = k - 1
            pl(k+1,i) = ( sqrt(two*k1+one) * x * pl(k,i) - k1 * pl(k1,i)
     x                  / sqrt(two*k1-one) ) * sqrt(two*k+one) /
     x                  dfloat(k)
   10    continue
   20 continue
      return
      end
      subroutine mesh (asect,nrange,ixmax,xi,npts,r)
*
* $Id: mesh.f,v 2.1 94/09/29 18:44:33 vmb Exp Locker: vmb $
*
      implicit real*8 (a-h,o-z)
      parameter ( half=0.5d0,two=2.d0 )
      dimension xi(ixmax),r(npts),asect(nrange+1)
*
*     mesh : mesh at which potentials are required
*            asect   = radius at end points of each subrange
*            nrange  = number of subranges
*            ixmax   = number of abscissae in quadrature scheme
*            xi      = abscissae
*            npts    = number of points in mesh
*            r(npts) = radii
*
      ra2 = asect(1)
      i1 = - ixmax
*
      do 20 nrang = 1, nrange
         ra1 = ra2
         ra2 = asect (nrang + 1)
         dr2 = half * (ra2 - ra1)
         rmid = ra1 + dr2
         i1 = i1 + ixmax
         i2 = i1 + 2 * ixmax + 1
         do 10 ix = 1, ixmax
            i1 = i1 + 1
            i2 = i2 - 1
            r(i1) = rmid - dr2 * xi(ix)
            r(i2) = rmid + dr2 * xi(ix)
   10    continue
   20 continue
      return
      end
      subroutine potl (nchan,ispch,nbasis,r,v,nmx,lamax,ion,lchl,cf)
*
* $Id: potl.f,v 2.1 94/09/29 18:44:35 vmb Exp Locker: vmb $
*
      implicit real*8 (a-h,o-z)
      parameter ( zero=0.0d0,one=1.0d0,two=2.0d0 )
      dimension r(nbasis),v(nchan,nchan,nbasis),cf(nmx,nmx,lamax),
     x          lchl(nmx),ispch(nmx)
*
*     potl : passes potential to r-matrix propagator routine rprop.
*            potentials are expanded in inverse powers of radial
*            distance ; expansion coefficients in cf
*
      do 40 ichan = 1, nchan
         i = ispch(ichan)
         el2 = dfloat(lchl(i) * ( lchl(i) + 1 ))
         do 30 jchan = 1, nchan
            j = ispch(jchan)
            do 20 ir = 1, nbasis
               vp = zero
               rr = one / r(ir)
               if ( i .eq. j ) vp = - two * dfloat(ion) * rr + el2 * rr
     x                              * rr
               do 10 k = 1, lamax
                  vp = vp + cf(i,j,k) * rr ** (k+1)
   10          continue
               v(ichan,jchan,ir) = vp
   20       continue
   30    continue
   40 continue
      return
      end
      subroutine prop (etotr,ns,ampa11,ampa12,eigen1,nch1,nhd1,
     x                 ampa21,ampa22,eigen2,nch2,nhd2,cr)
*
* $Id: prop.f,v 2.1 94/09/29 18:44:41 vmb Exp Locker: vmb $
*
      implicit real*8 (a-h,o-z)
      parameter ( zero=0.0d0,one=1.d0 )
      dimension  cr(*),ampa11(nch1,nhd1),ampa12(nch1,nhd1),eigen1(nhd1),
     x           ampa21(nch2,nhd2),ampa22(nch2,nhd2),eigen2(nhd2)
*
*      prop : propagates r-matrix across one sub-range
*
*             etotr       = impact energy
*             ampa1,ampa2 = surface amplitudes for current subrange
*             eigen       = eigenvalues of subrange hamiltonian
*             cr          = on entry the lower triangle of r-multiplied
*                           r-matrix, on exit the propagated matrix
*
      nchan = nch1 + nch2
      nchan1 = nchan * ( nchan + 1 ) / 2
      nchan2 = nchan1 + 1
      nchan3 = nchan1 + nchan * nchan
      nchan4 = nchan3 + 1
      nchan5 = nchan3 + nchan1
      nchan6 = nchan5 + 1
*     curly r-matrices stored in cr, r-multiplied r-matrix already in
*     first nchan1 locations, zero rest
*
      do 10 k = nchan2, nchan5
         cr(k) = zero
   10 continue
      do 40 kl = 1, nhd1
         edenom = one / ( eigen1(kl) - etotr )
         k = 0
         do 20 i = 1, nch1
            fac1 = ampa11(i,kl) * edenom
            fac2 = ampa12(i,kl) * edenom
            do 20 j = 1, i
               k = k + 1
               cr(k) = cr(k) + fac1 * ampa11(j,kl)
               cr(k+nchan3) = cr(k+nchan3) + fac2 * ampa12(j,kl)
   20    continue
*
*     non-symmetric matrix r12
*
         k = nchan1 - nch2
         do 30 j = 1, nch1
            fac1 = edenom * ampa12(j,kl)
            k = k + nch2
            do 30 i = 1, nch1
               k = k + 1
               cr(k) = cr(k) + ampa11(i,kl) * fac1
   30    continue
   40 continue
*
      if (ns .eq. 2) then
         ks1 = (nch1 * (nch1 + 1))/2
         ks12 = nch1 * nchan
         do 70 kl = 1, nhd2
            edenom = one / ( eigen2(kl) - etotr )
            k = ks1
            do 50 i = 1, nch2
               fac1 = ampa21(i,kl) * edenom
               fac2 = ampa22(i,kl) * edenom
               k = k + nch1
               do 50 j = 1, i
                  k = k + 1
                  cr(k) = cr(k) + fac1 * ampa21(j,kl)
                  cr(k+nchan3) = cr(k+nchan3) + fac2 * ampa22(j,kl)
   50       continue
*
*     non-symmetric matrix r12
*
            k = nchan1 + ks12
            do 60 j = 1, nch2
               fac1 = edenom * ampa22(j,kl)
               k = k + nch1
               do 60 i = 1, nch2
                  k = k + 1
                  cr(k) = cr(k) + ampa21(i,kl) * fac1
   60       continue
   70    continue
      end if
*
*     lower triangle of new r-multiplied r-matrix
*
      call jnvert (nchan,nchan,cr(1),cr(nchan2),cr(1),cr(nchan4),
     x             cr(nchan6),nchan2,nchan2)
      return
      end
      subroutine prppnt(rmatr,rlast,nschan,ebig,
     x                  ion,nranmx,idisc,idbbw,nampx,x,ix,
     x                  mleg,ixmx,delt,pt0,mcor1)
*
* $Id: prppnt.f,v 2.1 94/09/29 18:44:43 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 x,delta,delt,ebig,rmatr,rlast
      parameter ( nleg=10,delta=6.0d0 )
      parameter ( ixmax=5 )
      common /pptr/ ampa,eigen,bigvc,work
      dimension  x(*),ix(*),nschan(2),idisc(2)
      common /io/ iin, iout
*
*     print *,'start of prppnt'
      last = pt0
*     print *,' last = ',last
*
*   **** parameters nleg and ixmax should be set to the same values
*   **** in spropi
*   output parameters to cntrl for printing and opening files
*
      mleg = nleg
      ixmx = ixmax
      delt = delta
*
*   calculate nranmx and sector ends and modify rlast to lie on the
*   end of the last sector
*
      asect = 1
*     print *,' before sect in prppnt'
      call sect (ebig,ion,rmatr,rlast,delta,nranmx,x(asect))
*     print *,' after sect'
      lasect =  nranmx + 1
      last = pt0 + lasect
*
*  assign space only if there is to be propagation using the BBM method
*
      if (nranmx .gt. 0) then
         nch1 = nschan(1)
         nch2 = nschan(2)
         nmax = max (nch1,nch2)
         nhd = nmax * nleg
         nbigvc = ( nhd * ( nhd + 1 ) ) / 2
         nwork = nhd*nhd + 4 * nhd
         lpr = lasect + 2 + nbigvc + nwork + nhd
         nch1s = nch1 * nch1
         nch2s = nch2 * nch2
         nmx2 = nch1s + nch2s
         nmx = nch1 + nch2
         nampx = 4 * ixmax * nmx2 * nranmx
         neigen = 2 * ixmax * nmx * nranmx
         lreqd = lpr + nampx + neigen
*
*  calculation cannot be done in core use disc files
*
         if ( lreqd .gt. mcor1 .or. idbbw .gt. 0) then
            if ( idisc(1) .eq. 0 .or. idisc(2) .eq. 0 ) then
               write (iout,3010) lreqd,mcor1
               stop
            endif
            nampx = 4 * ixmax * nmx2
            neigen = 2 * ixmax * nmx
            lreqd = lpr + nampx + neigen
*
*  stop if not enough core even when using disc files
*
            if ( lreqd .gt. mcor1 ) then
               write (iout,3020) lreqd,mcor1
               stop
            endif
         endif
         eigen = last
         ampa = eigen + neigen
         bigvc = ampa + nampx
         work = bigvc + nbigvc
         last = work + nwork + 1
      end if
      pt0 = last
*
 3010 format (' storage overflow (prppnt) : needed =',i10,
     x        ' available =',i10,/,' disk required ')
 3020 format (' storage overflow (prppnt) : needed =',i10,
     x        ' available =',i10)
       return
       end
      subroutine rprop2 (rmat,etotr,nschan,nbasis,ampa,eigen,bigvec,
     x                   nrange,nampx,lbug,nran1,nran2,asect,irev)
*
* $Id: rprop2.f,v 2.1 94/09/29 18:44:44 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 rmat,etotr,ampa,eigen,bigvec,eps,a1,a2,ra2,dfloat,asect
      parameter ( eps=1.d-10 )
      common /io/ iin,iout
      common /sp/ idw,ifa,nsampr(2),inpr,nshd(2),nsamp(2),isdisc(2),
     x            mlegp
      dimension rmat(*),nbasis(*),ampa(*),eigen(*),bigvec(*),
     x          asect(*),nschan(2),nramp(2),nreig(2),isnhd(2),
     x          isnamp(2)
*    propagates rmat where the channels are rearranged according to two
*    possible spins.
*
      nch1 = nschan(1)
      nch2 = nschan(2)
      ns = 1
      if (nch2 .ne. 0) ns = 2
      mleg = mlegp
      nchan = nch1 + nch2
      if ( lbug .ne. 0 ) write (iout,1060)
      if ( nrange .eq. 0 ) return
      if (irev .eq. 0) then
         nrang1 = min (nran1,nran2)
         nrang2 = max (nran1,nran2)
         jsign = 1
         a1 = asect (nrang1)
         a2 = asect (nrang2+1)
      else
         nrang1 = max (nran1,nran2)
         nrang2 = min (nran1,nran2)
         jsign = -1
         a1 = dfloat (jsign) * asect (nrang1 +1)
         a2 = dfloat (jsign) * asect (nrang2)
      end if
      if ( idw .eq. 0 ) then
         do 10 is = 1,ns
            nramp(is) = nsamp(is) * ( nrang1 - 1 ) + 1
            nreig(is) = nshd(is) * ( nrang1 - 1 ) + 1
 10      continue
         nramp(2) = nramp(2) + nsampr(1)*2
         nreig(2) = nreig(2) + nshd(1)*nrange
      else
         nramp(1) = 1
         nreig(1) = 1
         nramp(2) = 1 + nsampr(1) * 2
         nreig(2) = 1 + nshd(1)
      endif
*
*     store lower triangle of r-multiplied r-matrix
*
      kmax = ( nchan * ( nchan + 1 ) ) / 2
      do 170 k = 1, kmax
         bigvec(k) = a1 * rmat(k)
  170 continue
*      nhd = mleg * nchan
*      ihd = nhd
*
*     propagate across subranges
*
      do 290 jr = nrang1, nrang2, jsign
         i1 = 1
         j1 = 1
         do 20 is = 1,ns
            isnhd(is) = nbasis(jr) * nschan(is)
            isnamp(is) = 2 *isnhd(is) * nschan(is)
            i2 = i1 + isnhd(is) -1
            j2 = j1 + isnamp(is) -1
            if ( idw .ne. 0 ) then
               read (isdisc(is),rec=jr) (ampa(j),j = j1,j2),
     x                   (eigen(i),i = i1,i2)
               i1 = i2 + 1
               j1 = j2 + 1
            endif
*
            ieg =nreig(is) - 1
            do 240 neig = 1,isnhd(is)
               ieg = ieg + 1
               if ( abs(eigen(ieg)-etotr) .lt. eps ) then
                  write (iout,1050) etotr,eigen(ieg)
                  stop
               endif
  240       continue
 20      continue
*
         if ( irev .eq. 0 ) then
            call prop (etotr,ns,ampa(nramp(1)),ampa(nramp(1)+nsampr(1)),
     x              eigen(nreig(1)),nch1,isnhd(1),ampa(nramp(2)),
     x              ampa(nramp(2)+nsampr(2)),eigen(nreig(2)),nch2,
     x              isnhd(2),bigvec)
            ra2 = asect (jr+1)
         else
            call prop (etotr,ns,ampa(nramp(1)+nsampr(1)),ampa(nramp(1)),
     x           eigen(nreig(1)),nch1,isnhd(1),ampa(nramp(2)+nsampr(2)),
     x           ampa(nramp(2)),eigen(nreig(2)),nch2,isnhd(2),bigvec)
            ra2 = asect (jr)
         end if
         if ( lbug .ne. 0 ) then
            write (iout,1070) ra2
            k = 0
            do 250 j = 1, nchan
               write (iout,1080) (bigvec(k+i),i = 1,j)
               k = k + j
  250       continue
         endif
         if ( idw .eq. 0 ) then
            do 30 is = 1,ns
               nramp(is) = nramp(is) + nsamp(is) * jsign
               nreig(is) = nreig(is) + nshd(is) * jsign
 30         continue
         endif
  290 continue
      do 320 k = 1, kmax
         rmat(k) = bigvec(k) / a2
  320 continue
      return
*
 1050 format (/,' energy',f14.7,' too close to eigenvalue',f14.7)
 1060 format (1h1,/,52x,' subroutine rprop2',/,52x,/)
 1070 format (/,5x,'new  r-matrix at',f10.5,3x,'a.u.',/)
 1080 format (8f14.7)
      end
      subroutine sect(ebig,ion,rmatr,rlast,delt,nrange,asect)
*
* $Id: sect.f,v 2.1 94/09/29 18:44:28 vmb Exp Locker: vmb $
*
      implicit integer (a - z)
      real*8 ebig,rmatr,rlast,delt,asect,delg,dela,k,a,ka,
     x     iok,x,one,oneh,dfloat
*
*   calculate end points of sectors for propagation and store
*   in asect.   each sector spans between two and three loops
*   of the wavefunction.
*
      dimension asect(*)
      parameter (one=1.0d0, oneh=1.5d0)
*
      k = sqrt(ebig)
      a = rmatr
      asect(1) = rmatr
      if (rlast .eq. rmatr) then
         nrange = 0
         return
      end if
      if (ion .gt. 0) then
         n = 1
 10      ka = k * a
         iok = dfloat (ion) / (k * ka)
         x = delt / ka
 20      delg = ka * (x + iok * log (one + x)) / delt
         if (delg .gt. oneh) then
            x = x / delg
            go to 20
         else
            n = n + 1
            a = (one + x) * a
            asect (n) = a
         end if
         if (a .lt. rlast) go to 10
         nrange = n-1
      else
         dela = delt / k
         nrange = int( (rlast-rmatr) / dela ) + 1
         do 30 i = 2,nrange + 1
            asect(i) = asect(i - 1) + dela
 30      continue
      end if
      rlast = asect(nrange + 1)
      return
      end
      subroutine setmtr (ra1,ra2,hmlt,wts,pl,v,ethr,bbloch,lbug,
     x                   nchan,ispch,ixmax,mleg,nhsize,nleg)
*
* $Id: setmtr.f,v 2.1 94/09/29 18:44:46 vmb Exp Locker: vmb $
*
      implicit real*8 (a-h,o-z)
      parameter ( zero=0.0d0,one=1.0d0,two=2.d0 )
      dimension ethr(*),ispch(nchan),v(nchan,nchan,*),wts(ixmax),
     x          pl(mleg,ixmax),hmlt(nhsize)
      common /io/ iin,iout
*
*     setmtr : hamiltonian matrix on a subrange
*              ra1,ra2  = end points of subrange
*              xi,wts   = quadrature abscissae and weights
*              pl       = normalized legendre polynomials at abscissae
*              v        = potential matrix evaluated at the abscissae
*              ethr     = threshold energies
*              bbloch   = parameter in bloch operator
*              lbug     = debug switch
*              hmlt     = lower triangle of hamiltonian matrix
*
      if ( lbug .eq. 1 ) write (iout,1000)
      nhs = 0
      ixmax2 = 2 * ixmax + 1
      r12 = one / ( ra2 - ra1 )
      r1b = zero
      r2b = zero
      if ( bbloch .ne. zero ) then
         r1b = r12 * bbloch / ra1
         r2b = r12 * bbloch / ra2
      endif
      rsq = two * r12 * r12
      do 50 nn = 1, nleg
        do 40 i = 1, nchan
            do 30 mm = 1, nn
               jmn = mod(mm+nn,2) + 1
               sgn = (-1) ** jmn
*     matrix elements not involving the potential ( closed form )
               rootmn = sqrt( ( two * dfloat(nn) - one ) * ( two *
     x                  dfloat(mm) - one ) )
               diag = - rootmn * ( r2b + sgn * r1b )
               if ( nn .eq. mm ) diag = diag + ethr(ispch(i))
               if ( jmn .eq. 1 ) diag = diag + rsq * rootmn *
c     x                                  dfloat(nn*(nn-1))
     x                                  dfloat(mm*(mm-1))
               jm = nchan
               if ( mm .eq. nn ) jm = i
               do 20 j = 1, jm
                  sum1 = zero
                  if ( i .eq. j ) sum1 = diag
*     matrix elements of potential using gauss-legendre quadrature.
                  do 10 ix = 1, ixmax
                     sum1 = sum1 + wts(ix) * pl(mm,ix) * pl(nn,ix) *
     x                      ( v(i,j,ixmax2-ix) - sgn * v(i,j,ix) )
   10             continue
*     store result in hmlt
                  nhs = nhs + 1
                  hmlt(nhs) = sum1
   20          continue
   30       continue
   40    continue
   50 continue
      if ( lbug .ne. 0 ) then
         write (iout,1010) nhsize
         write (iout,1020)
         write (iout,1030) (hmlt(nv),nv = 1,nhsize)
      endif
      return
*
 1000 format ('1',/,52x,' subroutine setmtr',/)
 1010 format (/,5x,' hamiltonian array size=',i5,//)
 1020 format (/,5x,' hamiltonian matrix elements',//)
 1030 format (1x,8f14.7)
      end
      subroutine sprop  (ethr,nmx,nchan,is,ispch,nrange,nbasis,asect,
     x                   bbloch,lbug,idisk,ampa,eigen,bigvec,
     x                   work,ipr,lamax,ion,lchl,cf,mleg,ixmax,
     x                   pl,xi,wts)
*
* $Id: sprop.f,v 2.1 94/09/29 18:44:48 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 ethr,asect,bbloch,ampa,eigen,bigvec,work,cf,
     x       pl,ra1,ra2,xi,wts
      dimension ethr(nmx),nbasis(nrange),ampa(*),
     x          eigen(*),xi(ixmax),wts(ixmax),pl(mleg,ixmax),
     x          bigvec(*),lbug(3),work(*),lchl(nmx),cf(nmx,nmx,lamax),
     x          asect(*),ispch(*)
      common /io/ iin,iout
      common /sp/ idw,ifa,nsampr(2),inpr,nshd(2),nsamp(2),isdisc(2),
     x            mlegp
*
      isdisc(is)  = idisk
      mlegp = mleg
      ifa = 0
      iprr = abs(ipr)
      if ( iprr .eq. 2 ) ifa = sign(1,ipr)
*     potl data
      mhd = 2 * nchan * ixmax
      nshd(is) = mhd
      mamp = nchan * mhd
      nsamp(is) = mamp
      nampr = mamp * nrange
      nsampr(is) = nampr
      if (idisk .ne. 0) then
         nampr = mamp
         nsampr(is) = nampr
         idw = 1
      endif
      ra2 = asect(1)
      nramp = 1
      nreig = 1
*     loop over subranges
      do 100 nrang = 1, nrange
         ra1 = ra2
         ra2 = asect(nrang+1)
         i1 = nramp
         i2 = nramp + nampr
         nleg = nbasis(nrang)
         nhd = nchan * nleg
         nhsize = nhd * ( nhd + 1 ) / 2
         namp = 2 * nchan * nhd
*
*     potential matrix including centrifugal term. mesh is stored
*     in eigen, potentials in ampa. potential is set up for entire range
*     (a,b) on the first pass unless a scratch disk used to store
*     amplitudes. potential then set up separately for each subrange.
*
         if ( idw .ne. 0 ) then
            npts = 2 * ixmax
            call mesh (asect(nrang),1,ixmax,xi,npts,eigen)
            if ( lbug(1) .ne. 0 ) write (iout,1040)
     x                            (eigen(i),i = 1,npts)
            call potl (nchan,ispch,npts,eigen,ampa,nmx,lamax,ion,
     x                 lchl,cf)
         else
            if ( nrang .eq. 1 ) then
               npts = 2 * ixmax * nrange
***** nmesh should only be called for first spin if keep it in eigen2
               call mesh (asect,nrange,ixmax,xi,npts,eigen)
               if ( lbug(1) .ne. 0 ) write (iout,1040)
     x                               (eigen(i),i = 1,npts)
               call potl (nchan,ispch,npts,eigen,ampa,nmx,lamax,ion,
     x                    lchl,cf)
            endif
         endif
*
*     matrix elements of hamiltonian on subrange (ra1,ra2) in bigvec
*     diagonalize subrange hamiltonian
*
         call setmtr (ra1,ra2,bigvec,wts,pl,ampa(nramp),ethr,bbloch,
     x                lbug(2),nchan,ispch,ixmax,mleg,nhsize,nleg)
         call ampltd (ra1,ra2,bigvec,ampa(i1),ampa(i2),eigen(nreig),
     x                work,work(nhd+1),work(nhd*nhd+nhd+1),
     x                lbug(3),nchan,nhsize,nhd,nleg)
         if ( idw .eq. 0 ) then
            nramp = nramp + mamp
            nreig = nreig + mhd
         else
            write (idisk,rec=nrang) (ampa(i),i = 1,namp),
     x            (eigen(i),i = 1,nhd)
         endif
  100 continue
      return
*
* 1030 format (/,' no propagation required for a=',f14.7,3x,'b=',f14.7,
*     x        3x,'ebig=',f14.7,/)
 1040 format (' potential matrix is evaluated at the following radii',
     x        /,(12f10.4))
      end
      subroutine spropi (ethr,nmx,nschan,ispch,nrange,nbasis,asect,
     x                   bbloch,lbug,idisk,ipr,lamax,ion,lchl,
     x                   cf,nampx,x)
*
* $Id: spropi.f,v 2.1 94/09/29 18:44:50 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 ethr,bbloch,cf,x,asect,xi,wts,pl
      parameter ( mleg=10,ixmax=5 )
      dimension ethr(nmx),nbasis(nrange),lbug(*),lchl(nmx),
     x          nschan(2),ispch(nmx),idisk(2),
     x          cf(nmx,nmx,lamax),x(*),asect(nrange+1),
     x          xi(ixmax),wts(ixmax),pl(mleg,ixmax)
      common /pptr/ ampa,eigen,bigvec,work
      common /sp/ idw,ifa,nsampr(2),inpr,nshd(2),nsamp(2),isdisc(2),
     x            mlegp
*
      common /io/ iin,iout
*
*     gauss-legendre quadrature abscissae and weights
*
      data xi(1)/0.973906528517172d0/,wts(1)/0.066671344308688d0/
      data xi(2)/0.865063366688985d0/,wts(2)/0.149451349150581d0/
      data xi(3)/0.679409568299024d0/,wts(3)/0.219086362515982d0/
      data xi(4)/0.433395394129247d0/,wts(4)/0.269266719309996d0/
      data xi(5)/0.148874338981631d0/,wts(5)/0.295524224714753d0/
*
      ifa = 0
      inpr = 0
      iprr = abs(ipr)
      if ( iprr .eq. 2 ) ifa = sign(1,ipr)
      nleg = 0
      do 20 ir = 1, nrange
         if ( nbasis(ir) .le. 0 ) nbasis(ir) = mleg
         if ( nbasis(ir) .gt. mleg ) nbasis(ir) = mleg
         if ( nbasis(ir) .gt. nleg ) nleg = nbasis(ir)
   20 continue
      if ( lbug(1) .ne. 0 ) then
         if ( inpr .eq. 0 ) then
            write (iout,1010) nleg
         else
            write (iout,1010) (nbasis(ir),ir = 1,nrange)
         endif
      endif
*    evaluate normalized legendre polynomials pl(x) at quadrature points
      call legndr (pl,xi,mleg,ixmax)
      nch1 = nschan(1)
      nch2 = nschan(2)
      nmax = max (nch1,nch2)
      mindim = 2 * nmax * nmax * max(nleg,ixmax)
      nampr = 2 * (nch1*nch1 + nch2*nch2) * ixmax  * nrange
*      check storage space for surface amplitudes and if disk available
      if ( idw .eq. 1  .or.  2*nampr .gt. nampx ) then
         if ( idisk(1) .eq. 0 .or. idisk(2) .eq. 0
     x       .or. mindim .gt. nampx ) then
            write (iout,1020) nmx,nleg,nrange
            stop
         endif
         idw = 1
      endif
      idisc = 0
      if (idw .eq. 1) idisc = idisk(1)
      is = 1
      call sprop (ethr,nmx,nch1,is,ispch,nrange,nbasis,asect,bbloch,
     x            lbug,idisc,
     x            x(ampa),x(eigen),x(bigvec),x(work),
     x            ipr,lamax,ion,lchl,cf,mleg,ixmax,pl,xi,wts)
      if (nch2 .gt. 0) then
         if (idw .eq. 0) then
            nhd1x = 2 * ixmax * nch1 * nrange
            ampas2 = ampa + 2*nch1*nhd1x
            eigen2 = eigen + nhd1x
         else
            idisc = idisk(2)
            ampas2 = ampa
            eigen2 = eigen
         end if
         is = 2
         call sprop (ethr,nmx,nch2,is,ispch(nch1+1),nrange,nbasis,asect,
     x               bbloch,lbug,idisc,x(ampas2),x(eigen2),
     x               x(bigvec),x(work),ipr,lamax,ion,lchl,cf,
     x               mleg,ixmax,pl,xi,wts)
      end if
      return
 1010 format (' number of basis functions included in each subrange is',
     x        /,(20i5),/)
 1020 format (' insufficient space in ampa and eigen for nchan=',i3,
     x        3x,'nleg=',i2,3x,'nrange=',i3,/)
      end
      subroutine compc(lrgl2,nspn2,npty2,lrgl,lrgs,lrgpi,ncases,iflg)
*
* $Id: compc.f,v 2.1 94/09/29 18:45:07 vmb Exp Locker: vmb $
*
      dimension lrgl(ncases),lrgs(ncases),lrgpi(ncases)
*
*    determine whether current block of data corresponds to a required case
*
      do 10  i=1,ncases
         if (lrgl(i) .eq. lrgl2 .and. lrgs(i) .eq. nspn2 .and.
     &       lrgpi(i) .eq. npty2) then
            lrgl(i) = -1
            iflg = 0
            go to 20
         end if
 10   continue
      iflg = -1
 20   return
      end
      subroutine fincas(maxcas,lrgs,lrgl,lrgpi,ncases,ix,mcore)
*
* $Id: fincas.f,v 2.1 94/09/29 18:44:59 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
c
c     dummy routine for use with serial H file
c
      dimension ix(mcore)
      return
      end
      subroutine rdch (nfth,nltarg,lchl,ichl,cf,eig,wmat,ntarg,nstat,
     &                nchan,n2,lamax,rmatr,ion,iprint,x,ix,mcor1,
     &                evhigh,nphys)
*
* $Id: rdch.f,v 2.1 94/09/29 18:45:00 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 cf,eig,wmat,rmatr,fac,flion,x,one,two,dfloat
      common /io/ iin,iout
      dimension nltarg(ntarg),lchl(nchan),ichl(nchan),
     x          cf(n2,lamax),eig(nstat),wmat(nstat,nchan),
     x          ix(mcor1),x(mcor1)
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0)
C
C  Change Oct 19
C
      real*8 evhigh
C
C  End of Oct 19 change
C
*
*    read arrays for current case
*
*     print *,'ntarg = ',ntarg
      read (nfth) (nltarg(i),i=1,ntarg)
*     print *,'reading the H-file'
*     print *,'nltarg'
*     print *,(nltarg(i),i=1,ntarg)
*     print *,'nchan = ',nchan
      read (nfth) (lchl(i),i=1,nchan)
*     print *,'lchl'
*     print *,(lchl(i),i=1,nchan)
*     print *,'lamax = ',lamax
      read (nfth) ((cf(i,j),i=1,n2),j=1,lamax)
*     print *,'cf'
*     print *,((cf(i,j),i=1,n2),j=1,lamax)
*     print *,'nstat = ',nstat
      read (nfth) (eig(i),i=1,nstat)
*     print *,'eig'
*     print *,(eig(i),i=1,nstat)
      read (nfth) ((wmat(i,j),j=1,nchan),i=1,nstat)
*     print *,'wmat'
*     print *,((wmat(i,j),j=1,nchan),i=1,nstat)
C
C  Change Oct 19
C
C  The goal of the following loop is to throw out
C  the spuriously high eigenvalues in the Hamiltonian
C  that originate in calculations with a large number
C  of pseudo states. Typically, these are recognizable
C  in that they form groups of eigenvalues closely
C  together, with big jumps between the groups.
C
      nstout = 0
      do 140 i=1,nstat/2
        if (two*(eig(i)-eig(nstat)).gt.evhigh) nstout = i
140   continue
      if (nstout.ne.0) then
         write (6,2531) nstout,eig(nstout)
2531     format(/,' eigenvalue # ',i5,'  with eigenvalue ',1pd14.6,
     +          /,' is the highest eigenvalue to be kept for the',
     +           ' calculation of the R-matrix',/)
         do 130 i=1,nstout
            do 135 j=1,nchan
               wmat(i,j) = zero
135         continue
130      continue
      endif
C
C  End of Oct 19 change
C
C Change Oct 17
C
      if (nphys.ne.0) then
         nphych = 0
         do 1 n = 1,nphys
            nphych = nphych + nltarg(n)
 1       continue
         do 2 i = nphych+1,nchan
            ii = nchan*(i-1)
            do 3 j = nphych+1,nchan
               do 4 lam = 1,lamax
                  cf(ii+j,lam) = 0.0
 4             continue
 3          continue
 2       continue
      endif
C
C  end of change
C
*
*    divide wmat by sqrt 2*rmatr for consistency with molecular code
*    and  scale for positive ions
*
      if (ion .gt. 0) then
         flion = dfloat (ion)
         fac = one / sqrt (two * rmatr * flion)
      else
         fac = one / sqrt (two * rmatr)
      end if
      do 5 j = 1,nchan
         do 5 i = 1,nstat
            wmat(i,j) = fac * wmat(i,j)
 5    continue
*
      i = 0
      do 10 it  = 1,ntarg
         nt = nltarg(it)
         do 20 n = 1,nt
            i = i+1
            ichl(i) = it
 20      continue
 10   continue
      if (iprint .gt. 0) then
         write (iout, 1000)
         write (iout, 1010) (i,ichl(i),lchl(i),i=1,nchan)
      end if
c
      if (ion .gt. 0) then
         fac = one / (flion * flion)
         do 30  i = 1,nstat
            eig(i) = fac * eig(i)
 30      continue
         fac = one
         do 40  lam = 2,lamax
            fac = fac * flion
            do 50  i = 1,n2
                  cf(i,lam) = fac * cf(i,lam)
 50         continue
 40      continue
      end if
 1000 format (12x,'channel',2x,'target',4x,'small'/
     x        12x,'index ',3x,'index  ',5x,'l'/)
 1010 format (7x,i8,i9,i10)
      return
      end
      subroutine rdhed(nfth,filh,nelc,nz,lrang2,ismax,ntarg,rmatr,
     x                 bbloch)
*
* $Id: rdhed.f,v 2.1 94/09/29 18:44:56 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 rmatr,bbloch
      character*48 filh
c
c     open a serial H file and read the opening few constants
c
      open (nfth, file = filh, access = 'sequential',
     x      status ='old', form = 'unformatted')
      read (nfth) nelc,nz,lrang2,ismax,ntarg,rmatr,bbloch
c
      return
      end
      subroutine rdhslp(nfth,ncases,icase,lrgl,lrgs,lrgpi,lrgl2,
     x                  nspn2,npty2,nchan,nstat,iflg,ix,mcore)
*
* $Id: rdhslp.f,v 2.1 94/09/29 18:44:57 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      dimension lrgl(ncases),lrgs(ncases),lrgpi(ncases),ix(*)
c
c     reads the first few constants for a particular lrgl2,nspn2
c     npty2 when ncases =0 and when ncases > 0
c
 50   read (nfth, end=90) lrgl2,nspn2,npty2,nchan,nstat,more2
      if (ncases .ne. 0) then
         call compc (lrgl2,nspn2,npty2,lrgl,lrgs,lrgpi,
     x               ncases,iflg)
         if (iflg .lt. 0) then
            read (nfth)
            read (nfth)
            read (nfth)
            read (nfth)
            read (nfth)
            go to 50
         end if
      end if
      return
c
 90   iflg=-1
      return
c
      end
      subroutine rdtar(nfth,etarg,ltarg,starg,cfbut,ntarg,lrang2,ion,
     x                 nbutx,maxcas,itlast,ix,mcor1)
*
* $Id: rdtar.f,v 2.1 94/09/29 18:45:04 vmb Exp Locker: vmb $
*
      implicit integer (a-z)
      real*8 etarg, cfbut, atry, etargr, d0, one, flion2, scale,dfloat
      parameter ( one=1.0d0, atry=2.0d0 )
      common /io/ iin,iout
      dimension etarg(ntarg),ltarg(ntarg),starg(ntarg),cfbut(3,lrang2),
     x          ix(mcor1)
*
      maxcas=0
*
*     maxcas is the total number of cases on a direct access H file
*     the value 0 is used to distinguish serial from direct access
*
      read (nfth) (etarg(i),i=1,ntarg)
      read (nfth) (ltarg(i),i=1,ntarg)
      read (nfth) (starg(i),i=1,ntarg)
      read (nfth) ((cfbut(i,j),i=1,3),j=1,lrang2)
C
C  change Oct 18
C
C      do 3 i= 1,3
C         read(37,*) cfbut(i,1)
C  3    continue
C
C  end of change
C
      write (iout,1000)
      if (ion .gt. 0) then
         flion2 = dfloat (ion * ion)
         scale = one / flion2
         do 5  i=1,ntarg
            etarg(i) = scale * etarg(i)
 5       continue
      end if
      d0 = etarg(1)
      do 10  i=1,ntarg
         etargr = atry * (etarg(i) - d0)
         write (iout,1010) i,ltarg(i),starg(i),etargr
 10   continue
      nbutx = 0
      do 20  i=1,lrang2
         nbutx = max( nbutx, -int(cfbut(3,i))/10000)
 20   continue
      if (ion .gt. 0 .and. nbutx .eq. 0) then
         scale = one
         do 25  m = 2,3
            scale = scale * flion2
            do 30  l = 1,lrang2
               cfbut(m,l) = scale * cfbut(m,l)
 30         continue
 25      continue
      end if
 1000 format(20x, 'target states'/20x,'*************',//
     x       10x,'index',5x,'total l',3x,'(2*s+1)',8x,'energy'/
     x       43x,'scaled ryd')
 1010 format(3x,3i10,7x,f12.6)
      return
      end

