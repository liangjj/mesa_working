inprog/bkwrd.f                                                                                      0000644 0000311 0000024 00000002174 05414300776 0013627 0                                                                                                    ustar 00bis                             users                           0000000 0000000                                                                                                                                                                        *deck bkwrd
      subroutine bkwrd(f,wt,intgl,n)
c***begin prologue     bkwrd
c***date written       930608   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose            backward integration of indefinite integral
c***
c***description        performs the set of indefinite integrals which result
c***                   from integrating a function on (a,b) as (a,r1), (a,r2)
c                      (a,r3)....(a,b) where ri are the quadrature points in
c                      the interval. the integral is initialized as intgl(0)
c                      which is either its last value or zero depending on the
c                      interval.   
c               
c***references
c
c***routines called
c
c***end prologue       bkwrd
c
      implicit integer (a-z)
      dimension f(n), wt(n-1,n), intgl(1:n)
      common /io/ inp, iout
      real*8 f, wt, intgl
      do 10 i=n-1,1,-1
         intgl(i)=intgl(i+1)
         do 20 j=1,n
            intgl(i)=intgl(i)+wt(i,j)*f(j)
   20    continue
   10 continue
      return
      end















                   which is either its last value or zero depending on the
c                      interval.   
c               
c***references
c
c***routines called
c
c***end prologue       bkwrd
c
      implicit integer (a-z)
      dimension f(n), wt(n-1,n), intgl(1:n)
      common /io/ inp, iout
      real*8 f, wt, intgl
      do 10 i=n-1,1,-1
         intgl(i)=intgl(i+1)
         doinprog/convr.f                                                                                      0000644 0000311 0000024 00000001404 05400726552 0013637 0                                                                                                    ustar 00bis                             users                           0000000 0000000                                                                                                                                                                        c $Header: convr.f,v 1.2 92/12/31 14:45:30 bis Exp $
*deck convr.f
c***begin prologue     convr
c***date written       921221   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           convr, link m6201, legendre quadrature
c***author             schneider, barry (nsf)
c***source             m6201
c***purpose            convert legendre quadrature from (-1,1) to (a,b)
c***references         none
c
c***routines called
c***end prologue       convr
      subroutine convr (a,b,r,wtr,npt)
      implicit integer (a-z)
      real*8 a, b, r, wtr, afac, bfac
      dimension r(npt), wtr(npt)
      afac=(b-a)*.5d0
      bfac=(b+a)*.5d0      
      do 10 i=1,npt
         r(i)=afac*r(i)+bfac
         wtr(i)=afac*wtr(i)
   10 continue 
      return
      end

   schneider, barry (nsf)
c***source             m6201
c***purpose            convert legendre quadrature from (-1,1) to (a,b)
c***references         none
c
c***routines called
c***end prologue       convr
      subroutine convr (a,b,r,wtr,npt)
      iinprog/drvint.f                                                                                     0000644 0000311 0000024 00000003624 05415103370 0014015 0                                                                                                    ustar 00bis                             users                           0000000 0000000                                                                                                                                                                        *deck drvint.f
c***begin prologue     drvint
c***date written       921223   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           drvint, link 2702, y(l,m), projections
c***author             schneider, barry (nsf)
c***source             m6202
c***purpose            driver for integral equation
c***references         none
c
c***routines called
c***end prologue       drvint
      subroutine drvint (psi,flm,j,y,wt,old,nrad,mval,lndex,nl,nm,ltop,
     1                   ldim,nshell,type)
      implicit integer (a-z)
      real*8 psi, flm, j, y, wt, old
      character*(*) type
      dimension psi(*), flm(*), j(*), y(*), nrad(nshell)
      dimension mval(nm), nl(nm), lndex(*), wt(*), old(*)
      common /io/ inp, iout
      locpsi=1
      locflm=1
      locj=1
      locy=1
      locind=1
      locwt=1
      call rzero(old,ldim)
      do 10 ns=1,nshell
c----------------------------------------------------------------------c
c         for each m, compute the forward radial integral              c
c--------------------------------------------------------------------- c     
         do 20 mu=1,nm
            if (type.eq.'gauss') then
                call fgauss(flm(locflm),psi(locpsi),old,j(locj),
     1                      y(locy),wt(locwt),lndex(locind),nl(mu),
     2                      ltop,nrad(ns))
                locflm=locflm+nl(mu)*nrad(ns)
                locpsi=locpsi+nl(mu)*nrad(ns)
                locind=locind+nl(mu)
            elseif (type.eq.'newton-cotes') then
                call fnewt(flm(locflm),psi(locpsi),j(locj),y(locy),
     1                     wt(locwt),lndex(locind),nl(mu),
     2                     ltop,nrad(ns))
            else
                    call lnkerr('quadrature type error')
            endif
   20    continue              
         locwt=locwt+nrad(ns)
         locj=locj+nrad(ns)*(ltop+1)
         locy=locy+nrad(ns)*(ltop+1)
      return
      end
d(ns)
                locind=locind+nl(mu)
            elseif (type.eq.'newton-cotes') then
                inprog/fgauss.f                                                                                     0000644 0000311 0000024 00000004570 05415101330 0013772 0                                                                                                    ustar 00bis                             users                           0000000 0000000                                                                                                                                                                        *deck fgauss
c***begin prologue     fgauss
c***date written       930608   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose            forward integration of integral equation
c***
c***description        performs the forward integration of the radial
c***                   integral equation involving a greens function and
c                      an inhomogeneity.
c               
c***references
c
c***routines called
c
c***end prologue       fgauss
c
      subroutine fgauss(flm,psilm,oldval,j,y,wt,lndex,nl,ltop,n)
      implicit integer (a-z)
      real*8 flm, psilm, wt, j, y, oldval
      dimension flm(n,nl), wt(n), psilm(n,nl), j(n,0:ltop), y(n,0:ltop)
      dimension lndex(nl), oldval(nl)
      common /io/ inp, iout
c----------------------------------------------------------------------c
c             initialize the integral at its first point using         c
c             the old value.                                           c
c----------------------------------------------------------------------c
      do 10 l=1,nl
         psilm(1,l)=oldval(l)+wt(1)*j(1,lndex(l))*flm(1,l)
   10 continue
c----------------------------------------------------------------------c
c           step up on the integral and get it at all the remaining    c
c           points.                                                    c
c----------------------------------------------------------------------c
      do 20 l=1,nl
         do 30 i=2,n
            psilm(i,l)=psilm(i-1,l)+wt(i)*j(i,lndex(l))*flm(i,l)
   30    continue
   20 continue
c----------------------------------------------------------------------c
c          store the last value in oldval so it can be used to         c
c          initialize the integral at the next call                    c
c----------------------------------------------------------------------c
      do 40 l=1,nl
         oldval(l)=psilm(n,l)
   40 continue
c----------------------------------------------------------------------c
c          now make the actual function by multiplying the intgral     c
c          by the prefactor in front of it.                            c
c----------------------------------------------------------------------c
      do 50 l=1,nl
         call vmul(psilm(1,l),psilm(1,l),y(1,lndex(l)),n)
   50 continue
      return
      end















------------------------c
      do 40 l=1,nl
         oldval(l)=psilm(n,l)
   40 continue
c---------------------------------------------inprog/forwrd.f                                                                                     0000644 0000311 0000024 00000002224 05414562725 0014021 0                                                                                                    ustar 00bis                             users                           0000000 0000000                                                                                                                                                                        *deck forwrd
      subroutine forwrd(flm,psilm,j,y,wt,lndex,nl,ltop,n)
c***begin prologue     forwrd
c***date written       930608   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose            forward integration of indefinite integral
c***
c***description        performs the set of indefinite integrals which result
c***                   from integrating a function on (a,b) as (a,r1), (a,r2)
c                      (a,r3)....(a,b) where ri are the quadrature points in
c                      the interval. the integral is initialized as intgl(0)
c                      which is either its last value or zero depending on the
c                      interval.   
c               
c***references
c
c***routines called
c
c***end prologue       forwrd
c
      implicit integer (a-z)
      dimension flm(n), wt(n-1,n), intgl(0:n-1)
      common /io/ inp, iout
      real*8 f, wt, intgl
      do 10 i=1,n-1
         intgl(i)=intgl(i-1)
         do 20 j=1,n
            intgl(i)=intgl(i)+wt(i,j)*f(j)
   20    continue
   10 continue
      return
      end















                 which is either its last value or zero depending on the
c                      interval.   
c               
c***references
c
c***routines called
c
c***end prologue       forwrd
c
      implicit integer (a-z)
      dimension flm(n), wt(n-1,n), intgl(0:n-1)
      common /io/ inp, iout
      real*8 f, wt, intgl
      do 10 i=1,n-1
         intgl(iinprog/legend.f                                                                                     0000644 0000311 0000024 00000004157 05344542546 0013764 0                                                                                                    ustar 00bis                             users                           0000000 0000000                                                                                                                                                                        c $Header: legend.f,v 1.3 93/01/20 20:59:18 bis Exp $
*deck legend.f
c***begin prologue     legend
c***date written       921223   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           legend, link 6201, legendre functions
c***author             schneider, barry (lanl)
c***source             m6201
c***purpose            legendre functions
c***description        calculation of p(l,m) functions
c***references         none
c
c***routines called
c***end prologue       legend
      subroutine legend (plm,x,dfct,ddfct,npt,lmax,m,maxfac)
      implicit integer (a-z)
      real *8 plm, x, dfct, ddfct, fm, facx, f1, norm
      dimension plm(npt,m:lmax), x(npt)
      dimension  dfct(0:maxfac), ddfct(0:maxfac)
c----------------------------------------------------------------------c
c           start recursion with plm(m,m) and plm(m+1,m)               c
c                      and recur upward                                c
c----------------------------------------------------------------------c
      do 10 i=m,lmax
         do 10 j=1,npt
            plm(j,i)=0.d+00
   10 continue
      fm=.5d+00*m
      do 20 i=1,npt
         facx=1.d+00
         if (fm.ne.0.d+00) then
             facx=(1.d+00-x(i)*x(i))**fm
         endif
         plm(i,m)=ddfct(m)*facx
   20 continue
      if (lmax.ne.m) then
          mm=m+m+1
          mpls1=m+1
          do 30 i=1,npt
             plm(i,mpls1)=mm*x(i)*plm(i,m)
   30     continue
          if (lmax.ne.mpls1) then
              lind=m+2
              n1=2
              n2=m+m+3
              n3=n2-2
              do 50 i=lind,lmax
                 ii=i-1
                 jj=i-2
                 do 40 j=1,npt
                    f1=n2*x(j)*plm(j,ii)-n3*plm(j,jj)
                    f1=f1/n1
                    plm(j,i)=f1
   40            continue
                 n1=n1+1
                 n2=n2+2
                 n3=n3+1
   50       continue
          endif
      endif
c     normalize
      do 60 l=m,lmax
         norm=sqrt((l+l+1)*dfct(l-m)/(2.d0*dfct(l+m)))
         do 70 j=1,npt
            plm(j,l)=norm*plm(j,l)
   70    continue
   60 continue
      return
c
      end


j=i-2
                 do 40 j=1,npt
                    f1=n2*x(j)*plm(j,ii)-n3*plm(j,jj)
                    f1=f1/n1
                    plm(j,i)=f1
   40            continue
                 n1=n1+1
                 n2=n2+2
                 n3=n3+1
   50       continue
          endif
      endif
c     normalize
      do 60 l=m,lmax
         norm=sqrt((l+l+1)*dfct(l-m)/(2.d0*dfct(l+m)))
       inprog/legint.f                                                                                     0000644 0000311 0000024 00000002436 05415033455 0013777 0                                                                                                    ustar 00bis                             users                           0000000 0000000                                                                                                                                                                        c $Header $
*deck legint.f
c***begin prologue     legint
c***date written       921221   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           legint, link m6201, legendre integrals
c***author             schneider, barry (nsf)
c***source             m6201
c***purpose            calculate overlap integrals of legendre functions
c***references         none
c
c***routines called
c***end prologue       legint
      subroutine legint (plm,phifn,scr,nth,nph,lmax,m)
      implicit integer (a-z)
      real*8 plm, phifn, scr, phint
      character*80 title
      character*2 itoc
      dimension plm(nth,m:lmax), phifn(nph)
      dimension scr(m:lmax,m:lmax)
      common /io/ inp, iout
      dim=lmax-m+1
      call rzero(scr,dim*dim)
      do 10 l1=m,lmax
         do 20 l2=l1,lmax
            do 30 n=1,nth
               scr(l1,l2)=scr(l1,l2) + plm(n,l1)*plm(n,l2)
   30       continue
            scr(l2,l1)=scr(l1,l2)          
   20    continue
   10 continue 
      title='legendre overlap integrals for m = '//itoc(m)
      call prntfm(title,scr,dim,dim,dim,dim,iout)
      phint=0.d0
      do 40 n=1,nph
         phint = phint + phifn(n)*phifn(n)
   40 continue
      write (iout,1) m, phint
    1 format(/,5x,'phi overlap integral for m = ',i4,1x,e15.8)    
      return
      end

ax
            do 30 n=1,nth
               scr(l1,l2)=scr(l1,l2) + plm(n,l1)*plm(n,l2)
   30       continue
            scr(l2,l1)=scr(l1,l2)          
   20    continue
   10 continue 
      title='legendre overlap integralsinprog/lmdcmp.f                                                                                     0000644 0000311 0000024 00000003725 05415034322 0013765 0                                                                                                    ustar 00bis                             users                           0000000 0000000                                                                                                                                                                        *deck lmdcmp.f
c***begin prologue     lmdcmp
c***date written       921223   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           lmdcmp, link 2702, y(l,m), projections
c***author             schneider, barry (nsf)
c***source             m6201
c***purpose            decompose a function into its legendre projections
c***description        a three dimensional function is decomposed into
c***                   its projection onto p(l,m) phi(m) angular
c***                   functions. the angular functions in the phi
c***                   variable are the normalized sin(m*phi) and
c***                   cos(m*phi). it is assumed that the function is
c***                   stored with no blank spaces in any column direction.
c***references         none
c
c***routines called
c***end prologue       lmdcmp
      subroutine lmdcmp (f,plm,phifn,r,flm,scr,mval,nl,nr,nth,nph,
     1                   nm,prnt)
      implicit integer (a-z)
      real *8 f, plm, phifn, flm, scr, r
      logical prnt
      character*2 itoc, tmp
      character*80 title
      dimension f(nr,nth,nph), plm(nth,*), phifn(nph,*)
      dimension scr(nr,nth), flm(nr,*), r(nr)
      dimension mval(nm), nl(nm)
      common /io/ inp, iout
      plmcnt=1
      flmcnt=1
      do 10 mu=1,nm
         call ebc(scr,f,phifn(1,mu),nr*nth,nph,1)
         if (prnt) then
             write(iout,1) mval(mu)
             tmp=itoc(mval(mu))
             len=length(tmp)
         endif
         call ebc(flm(1,flmcnt),scr,plm(1,plmcnt),nr,nth,nl(mu))
         call vmmul(r,flm(1,flmcnt),flm(1,flmcnt),nr,
     1              nl(mu))
         if (prnt) then
             title='p(l,'//tmp(1:len)//') decomposition'
             call prntfm(title,flm(1,flmcnt),nr,nl(mu),nr,
     1                   nl(mu),iout)
         endif                    
         flmcnt=flmcnt+nl(mu)     
         plmcnt=plmcnt+nl(mu)
   10 continue
    1 format(/,5x,'the legendre decomposition matrix for m = ',i3)
      return
      end
flm(1,flmcnt),scr,plm(1,plmcnt),nr,nth,nl(minprog/miscfn.f                                                                                     0000644 0000311 0000024 00000000535 05344542546 0014001 0                                                                                                    ustar 00bis                             users                           0000000 0000000                                                                                                                                                                        c $Header: miscfn.f,v 1.1 92/12/23 17:37:45 bis Exp $
*deck miscfn.f
      subroutine miscfn(phi,cphi,sphi,nphi)
      implicit integer (a-z)
      real *8 phi, cphi, sphi
      dimension phi(nphi), cphi(nphi), sphi(nphi)
      do 10 pnt=1,nphi
         cphi(pnt)=cos(phi(pnt))
         sphi(pnt)=sin(phi(pnt))
   10 continue
      return
      end
=plmcnt+nl(mu)
   10 continue
    1 format(/,5x,'the legendre decomposition matrix for m = ',i3)
      return
      end
flm(1,flmcnt),scr,plm(1,plmcnt),nr,nth,nl(minprog/plmtof.f                                                                                     0000644 0000311 0000024 00000003076 05411641564 0014021 0                                                                                                    ustar 00bis                             users                           0000000 0000000                                                                                                                                                                        c $Header: plmtof.f,v 1.2 92/12/31 14:44:00 bis Exp $
*deck plmtof.f
c***begin prologue     plmtof
c***date written       921223   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           plm, link 2702, y(l,m), projections
c***author             schneider, barry (nsf)
c***source             m6201
c***purpose            reconstruct a three dimensional function from
c***                   its radial y(l,m) projection and p(l,m)*phifn(m)
c***description        the input projections and the angular functions
c***                   are summed to recreate the original three 
c***                   dimensional function. this can be used to make
c***                   a pointwise least squares comparison with the
c***                   exact original function.
c***references         none
c
c***routines called
c***end prologue       plmtof
      subroutine plmtof (f,plm,phifn,flm,scr,lmax,m,nr,nth,nph,
     1                   num,prnt)
      implicit integer (a-z)
      real*8 f, plm, phifn, flm, scr
      logical prnt
      character*2 itoc
      character*80 title
      dimension f(nr,nth,nph,2), plm(nth,m:lmax), phifn(nph,num)
      dimension flm(nr,m:lmax,num), scr(nr,nth,num)
      common /io/ inp, iout
      dim=lmax-m+1
      do 10 i=1,num
         call ebct(scr(1,1,i),flm(1,m,i),plm(1,m),nr,dim,nth)
   10 continue
      call apbct(f(1,1,1,2),scr,phifn,nr*nth,num,nph)
      if (prnt) then
          title='exact and fitted function for m = '//itoc(m)
          prd=nr*nth*nph
          call prntfm(title,f,prd,2,prd,2,iout)
      endif
      return
      end
title
      dimension f(nr,nth,nph,2), plm(nth,m:lmax), phifn(nph,num)
      dimension flm(nr,m:lmax,num), scr(nr,nth,num)
      common /io/ inp, iout
      dim=lmax-m+1
      do 10 i=1,num
         call ebct(scr(1,1,i),flm(1,m,i),plm(1,m),nr,dim,nth)
   10 continue
      call apbct(f(1,1,1,2),scr,phifn,nr*nth,num,nph)
      if (prnt) then
          title='exact and fitted function for m = '//itoc(m)
          prd=nr*nth*nph
          call prntfminprog/scalfn.f                                                                                     0000644 0000311 0000024 00000001330 05377745105 0013764 0                                                                                                    ustar 00bis                             users                           0000000 0000000                                                                                                                                                                        c $Header: scalfn.f,v 1.2 92/12/31 14:45:30 bis Exp $
*deck scalfn.f
c***begin prologue     scalfn
c***date written       921221   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           scalfn, link m6201, legendre integrals
c***author             schneider, barry (nsf)
c***source             m6201
c***purpose            scale a matrix of function values
c***references         none
c
c***routines called
c***end prologue       scalfn
      subroutine scalfn (f,wt,n,l,m)
      implicit integer (a-z)
      real*8 f, wt
      dimension f(n,m:l), wt(n)
      common/io/inp,iout
      do 10 i=m,l
         do 20 j=1,n
            f(j,i) = wt(j)*f(j,i)
   20    continue
   10 continue 
      return
      end

 m6201, legendre integrals
c***author             schneider, barry (nsf)
c***source             m6201
c***purpose            scale a matrix of function values
c***references         none
c
c***routines called
c***end prologue       scalfn
      subroutine scalfn (f,wt,n,l,m)
      implicit integinprog/scmphi.f                                                                                     0000666 0000311 0000024 00000001662 05414654004 0014002 0                                                                                                    ustar 00bis                             users                           0000000 0000000                                                                                                                                                                        c $Header: scmphi.f,v 1.2 92/12/24 15:30:50 bis Exp $
*deck scmphi.f
      subroutine scmphi(cphi,sphi,phifn,n,mu)
      implicit integer (a-z)
      common /io/ inp, iout
      dimension cphi(n), sphi(n), phifn(n)
      real *8 cphi, sphi, phifn, const
      real *8 pi, twopi
      complex*16 cmp
      data pi / 3.14159265358979323846d+00 /
      twopi=2.d0*pi
      const=sqrt(1.d0/twopi)
      if(mu.eq.0) then
         do 10 i=1,n
            phifn(i)=const
   10    continue
      else
         const=const*sqrt(2.d0)
         absmu=abs(mu)
         if (mu.gt.0) then
             do 20 i=1,n
                cmp=(dcmplx(cphi(i),sphi(i)))**absmu
                phifn(i)=real(cmp)*const
   20        continue
         else
             do 30 i=1,n
                cmp=(dcmplx(cphi(i),sphi(i)))**absmu                         
                phifn(i)=imag(cmp)*const
   30        continue
         endif
      endif
      return
      end
,n
            phifn(i)=const
   10    continue
      else
         const=consinprog/sqrtwt.f                                                                                     0000664 0000311 0000024 00000001411 05377745503 0014066 0                                                                                                    ustar 00bis                             users                           0000000 0000000                                                                                                                                                                        c $Header: sqrtwt.f,v 1.1 92/12/24 16:21:36 bis Exp $
*deck sqrtwt.f
c***begin prologue     sqrtwt
c***date written       921221   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           sqrtwt, link m6201, legendre integrals
c***author             schneider, barry (nsf)
c***source             m6201
c***purpose            calculate square root of integration weights.
c***references         none
c
c***routines called
c***end prologue       sqrtwt
      subroutine sqrtwt (wtthe,wtphi,nth,nph)
      implicit integer (a-z)
      real*8 wtthe, wtphi
      dimension wtthe(nth), wtphi(nph)
      do 10 n=1,nth
         wtthe(n)=sqrt(wtthe(n))          
   10 continue
      do 20 n=1,nph
         wtphi(n)=sqrt(wtphi(n))
   20 continue
      return
      end

 schneider, barry (nsf)
c***source             m6201
c***purpose            calculate square root of integration weights.
c***references         none
c
c***routines called
c***end prologue       sqrtwt
      subroutine sqrtwt (wtthe,wtphi,nth,nph)inprog/srtleg.f                                                                                     0000666 0000311 0000024 00000001362 05406676344 0014030 0                                                                                                    ustar 00bis                             users                           0000000 0000000                                                                                                                                                                        *deck srtleg.f
c***begin prologue     srtleg
c***date written       930612   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           legend, link 6201, legendre functions
c***author             schneider, barry (nsf)
c***source             m6201
c***purpose            legendre functions
c***description        sort of p(l,m) functions
c***references         none
c
c***routines called
c***end prologue       legend
      subroutine srtleg(plmin,plmout,lval,npt,lmax,m,nl)
      implicit integer (a-z)
      real *8 plmin, plmout
      dimension plmin(npt,m:lmax), plmout(npt,nl), lval(nl)
c
      do 10 i=1,nl
         do 20 j=1,npt
            plmout(j,i)=plmin(j,lval(i))
   20    continue
   10 continue
      return
c
      end




ce             m6201
c***purpose            legendre functions
c***description        sort of p(l,m) functions
c***references         none
c
c***routines called
c***end prologue       legend
      subroutine srtleg(plmin,plmout,lval,npt,lmax,m,nl)
      implicit integerinprog/ylm.f                                                                                        0000644 0000311 0000024 00000043727 05415124526 0013326 0                                                                                                    ustar 00bis                             users                           0000000 0000000                                                                                                                                                                        *deck ylm.f
c***begin prologue     m6201
c***date written       921219   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6201, link 6201, spherical harmonics
c***author             schneider, b. i.(lanl/nsf)
c***source             m6201
c***purpose            solve inhomogeneous wave equation by expansion
c***                   in spherical harmonics.
c***description        solution of the equation ( Del**2 + k**2 ) V = Rho
c***                   is calculated by expansion in spherical harmonics.
c***                   we allow k to be zero to handle the poisson equation.
c*** 
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       m6201
      program ylm 
      implicit integer (a-z)
      parameter ( dimcen=10 , dimshl=100, dimtmp=200, dimene=50, acc=30)
      character*4 chvr1
      character*4 itoc
      character*8 cpass
      character*20 radqud, chrkey
      character*128 fillam
      character*1600 card
      character*4096 ops
      character*80 title
      logical logkey, prnt, prntlm, posinp, tstnrm, tstint
      real*8 z, dummy, center, alpha, rr, rshel, rbox, energy
      real*8 tmp
      dimension z(1), dummy(2)
      dimension nrad(dimshl), rshel(dimshl), alpha(dimcen)
      dimension center(3,dimcen), mval(dimtmp), ibuf(dimtmp), nl(dimtmp)
      dimension energy(dimene), k, tmp
      common a(1)
      common /io/ inp, iout
      common /memory / ioff
      equivalence (z,a)
      call drum
      call iosys ('read character options from rwf',-1,0,0,ops)
      prnt=logkey(ops,'print=m6201=ylm',.false.,' ')
      prntlm=logkey(ops,'print=m6201=lm-decomposition',.false.,' ')
      tstnrm=logkey(ops,'m6201=test-normalization',.false.,' ')
      tstint=logkey(ops,'m6201=test-yukawa-integral',.false.,' ')
      write (iout,1)
      call iosys ('read character "linear algebraic filename" from rwf',
     1                 -1,0,0,fillam)
      call iosys ('open lamdat as old',0,0,0,fillam)
c----------------------------------------------------------------------c
c          if m is positive it denotes cosine like phi behavior.       c
c          negative m is sinelike behavior.                            c                         
c----------------------------------------------------------------------c      
      if ( posinp ('$ylm',cpass) ) then
           call cardin(card)
           nm=intkey(card,'number-of-m-values',1,' ')   
           call iosys ('write integer "number m values" to lamdat',
     1                  1,nm,0,' ')
           call intarr(card,'m-values',mval,nm,' ')
           call iosys ('write integer "m values" to lamdat',
     1                  nm,mval,0,' ')
           call intarr(card,'number-l-values-for-each-m',nl,nm,' ')
           call iosys ('write integer "number of l values for each m"'//
     1                 ' to lamdat',nm,nl,' ')
           lmax=0
           mumax=0
           totl=0
           do 10 mmax=1,nm
              mumax=max(mumax,abs(mval(mmax)))
              totl=totl+nl(mmax)
              call intarr(card,'l-values-for-m-'//itoc(mval(mmax)),
     1                    ibuf,nl(mmax),' ')
              call order(ibuf,ibuf,nl(mmax),'integer')
              call iosys ('write integer "l values for m-'//
     1                    itoc(mval(mmax))//'" to lamdat',nl(mmax),
     2                    ibuf,' ')
              do 20 l=1,nl(mmax)
                 lmax=max(lmax,ibuf(l))
   20         continue
   10      continue
      endif
c
c            we may drop lm pairs in a given region by using zero
c            radial quadrature points for that pair in that region
c
      if ( posinp ('$radial',cpass) ) then
           call cardin(card)
           nshell=intkey(card,'number-of-radial-shells',1,' ')
           radqud=chrkey(card,'type-radial-quadrature',
     1                   'newton-cotes',' ')
           call iosys ('write integer "number of radial shells" '//
     1                 'to lamdat',1,nshell,0,' ')
           write (iout,3) nshell
           call intarr(card,'number-of-radial-points-per-shell',nrad,
     1                 nshell,' ')
           call iosys ('write integer "number radial points per '//
     1                 'shell" to lamdat',nshell,nrad,' ') 
           write (iout,4) ( nrad(ns),ns=1,nshell)
           call fparr(card,'shell-radii',rshel,nshell,' ')
           write (iout,5) (rshel(ns),ns=1,nshell)
           maxr=0
           totr=0
           do 30 ns=1,nshell
              maxr=max(maxr,nrad(ns))
              totr=totr+nrad(ns)
   30      continue
           rbox=rshel(nshell)
           lenwt=totr
           if (radqud.eq.'newton-cotes') then
               lenwt=0
               do 40 ns=1,nshell
                  lenwt=lenwt+nrad(ns)*(nrad(ns)-1)
   40          continue
           endif                 
      endif
      if ( posinp ('$energy',cpass) ) then
           call cardin(card)
           nen=intkey(card,'number-of-energies',1,' ')
           call fparr(card,'energies',energy,nen,' ')
      endif                                 
      write (iout,6) lmax, mumax
      call iosys ('write integer "maximum l in ylm" to lamdat',1,
     1             lmax,0,' ')
      call iosys ('write integer "maximum m in ylm" to lamdat',1,
     1             mumax,0,' ')
      write(iout,8) nen, (energy(ii),ii=1,nen)
c       
c             based on the values of lmax and mumax we can compute
c             the maximum number of quadrature points in theta and
c             phi to integrate a polynomial integrand exactly of a
c             given order.
c
      maxfac=lmax+mumax+10
      maxth=lmax+1
      maxph=2*mumax+1
      maxang=max(maxth,maxph)
      maxall=max(maxang,maxr)
      totpt=maxph*maxth*totr
      if (posinp('$yukawa',cpass)) then
          call cardin(card)
          ncen=intkey(card,'number-yukawa-centers',1,' ')
          do 50 cen=1,ncen
             call intarr(card,'position-center-'//itoc(cen),
     1                   center(1,cen),3,' ')
             alpha(cen)=fpkey(card,'yukawa-exponent-center-'
     1                       //itoc(cen),1.d0,' ')
   50     continue           
      endif
c----------------------------------------------------------------------c
c                     memory allocation                                c
c                          RE DO                                       c
c----------------------------------------------------------------------c
      lplus=lmax+1
      lsq=lplus*lplus
      maxscr=max(maxall,totl,lplus*maxth,2*maxr*maxth,
     1           2*maxr*(ltop1+1))
      if (tstnrm) then
          maxscr=max(maxscr,lsq,4)
      endif               
      tmp=lmax*acc
      ltop=lmax+sqrt(tmp)
      ltop=max(strtl,lmax)
      ltop1=ltop+1
      ioff=1
      do 60 i=1,2
          dfct=ioff
          ddfct=dfct+maxfac
          theta=ddfct+maxfac
          wtthe=theta+maxth
          phi=theta+maxth+maxth
          wtphi=phi+maxph
          rad=wtphi+maxph+maxph
          wtrad=rad+totr+totr
          sphi=wtrad+lenwt+lenwt
          cphi=sphi+maxph
          phifn=cphi+maxph
          plm=phifn+nm*maxph
          yuk=plm+maxth*totl
          j=yuk+totpt
          y=j+totr*ltop1
          wron=y+totr*ltop1
          flm=wron+ltop1
          psilm=flm+totr*totl
          lndex=wpadti(psilm+totr*totl)
          jp=iadtwp(lndex+totl)
          yp=jp+maxr*ltop1
          scr=yp+maxr*ltop1
          plmscr=jp
          words=wpadti(scr+maxscr)
          if (i.eq.1) then
              call iosys ('read integer maxsiz from rwf',1,
     1                     canget,0,' ')
              if (words.gt.canget) then
                  call lnkerr('not enough memory. will quit')
              endif
              call iosys ('write integer maxsiz to rwf',1,
     1                     words,0,' ')
              call getscm(words,z,ngot,'m6201',0)      
          endif
   60 continue
c----------------------------------------------------------------------c          
c             calculate the gauss-legendre quadrature points for the   c
c             plm functions and the simpson quadrature points for      c 
c             the phifn integration.                                   c
c----------------------------------------------------------------------c
      call gaussq('legendre',maxth,0.d0,0.d0,0,dummy,z(scr),
     1             z(theta),z(wtthe))
      call gaussq('simpson',maxph,0.d0,0.d0,0,dummy,z(scr),
     1             z(phi),z(wtphi))
c----------------------------------------------------------------------c     
c     the nature of the weights and points is very different for       c
c     newton-cotes and gauss quadrature. in the newton-cotes           c
c     formulas the end ranges of the points overlap and the weights    c
c     for each region are a matrix of dimension n*(n-1).               c 
c----------------------------------------------------------------------c 
      rr=0.d0
      locr=rad
      locwtr=wtrad
      do 70 ns=1,nshell
         if (radqud.eq.'legendre') then
             call gaussq('legendre',nrad(ns),0.d0,0.d0,0,z(scr),
     1                    z(locr),z(locwtr))
             call convr(rr,rshel(ns),z(locr),z(locwtr),nrad(ns))
             locwtr=locwtr+nrad(ns)
         elseif (radqud.eq.'newton-cotes') then
             call necote(rr,rshel(ns),z(locr),z(locwtr),nrad(ns))
             locwtr=locwtr+nrad(ns)*(nrad(ns)-1)
         else
             call lnkerr('error in radial quadrature type')
         endif
         locr=locr+nrad(ns)
         rr=rshel(ns)
   70 continue     
c----------------------------------------------------------------------c 
c             calculate factorials and double factorials needed        c
c             for computation of legendre functions.                   c   
c----------------------------------------------------------------------c 
      call fact(z(dfct),z(ddfct),maxfac)
c----------------------------------------------------------------------c
c             begin section of code to calculate the spherical         c
c             functions. by setting appropriate flags the              c
c             normalization integrals can be examined.                 c
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c             calculate cosine and sine of phi.                        c
c             they are used over and over again in plm routine.        c
c----------------------------------------------------------------------c
      call miscfn(z(phi),z(cphi),z(sphi),maxph)
c----------------------------------------------------------------------c
c             calculate the yukawa potential                           c
c----------------------------------------------------------------------c
      locr=rad
      locyuk=yuk
      do 80 ns=1,nshell
         call yukawa(z(locyuk),alpha,center,z(locr),z(theta),z(cphi),
     1               z(sphi),ncen,nrad(ns),nthet,nphi)
         locr=locr+nrad(ns)
         locyuk=locyuk+nrad(ns)*nthet*nphi
  80  continue         
c----------------------------------------------------------------------c
c              need copies of weights                                  c
c----------------------------------------------------------------------c
      call copy(z(wtthe),z(wtthe+maxth),maxth)
      call copy(z(wtphi),z(wtphi+maxph),maxph)
c----------------------------------------------------------------------c
c             compute the legendre and phifn functions.                c
c----------------------------------------------------------------------c 
      locplm=plm
      locphm=phifn
      loclnd=lndex
      do 90 m=1,nm
         mu=mval(m)
         absmu=abs(mu)
         call legend(z(plmscr),z(theta),z(dfct),z(ddfct),maxth,
     1               lmax,absmu,maxfac)
         call strleg(z(plmscr),z(locplm),a(loclnd),maxth,lmax,
     1               absmu,nl(m))
         loclnd=loclnd+nl(m)
         title='p(l,m) functions for m = '//itoc(mu)
         lentit=length(title)
         wds=maxth*nl(m)
         if (prnt) then
             call prntfm(title,z(locplm),maxth,nl(m),maxth,nl(m),iout)
         endif
         call scmphi(z(cphi),z(sphi),z(locphm),maxph,mu)
         title='phi angular function for m = '//itoc(mu)
         lentit=length(title)
         if (prnt) then
             call prntfm(title,z(phifn),maxph,1,maxph,
     1                   1,iout)
         endif
         if (tstnrm) then
c
c               calculate square root of weights
c
             call sqrtwt(z(wtthe+maxth),z(wtphi+maxph),maxth,maxph)
c
c               scale functions by the new weights in order to
c               make  normalization test integral easy.
c
             call scalfn(z(locplm),z(wtthe+maxth),maxth,nl(m),1)
             call scalfn(z(locphm),z(wtphi+maxph),maxph,1,1)
             call legint(z(locplm),z(locphm),z(scr),maxth,maxph,
      1                  nl(m),1)
c
c               calculate inverse of new weights
c
             call vinv(z(wtthe+maxth),z(wtthe+maxth),maxth)
             call vinv(z(wtphi+maxph),z(wtphi+maxph),maxph)
c
c               scale functions by inverse of weights to return them
c               to original state.
c
             call scalfn(z(locplm),z(wtthe+maxth),maxth,nl(m),1)
             call scalfn(z(locphm),z(wtphi+maxph),maxph,1,1)
         endif
c
c           scale functions by weights to make p(l,m) projections
c           easily vectorized.
c           write out scaled weights for later use.
c 
         call scalfn(z(locplm),z(wtthe),maxth,nl(m),1)
         call scalfn(z(locphm),z(wtphi),maxph,1,1)
         chvr1=itoc(mu) 
         len=length(chvr1)
         title='"weight scaled p(l,'//chvr1(1:len)//')"'
         lentit=length(title)
         call iosys ('write real '//title//' to lamdat',
     1                wds,z(locplm),0,' ')
         write(iout,*) 'writing iosys file ', title(1:lentit)
         write(iout,*) 'words written = ',wds
         title='"weight scaled phi('//chvr1(1:len)//')"'
         lentit=length(title)
         call iosys ('write real '//title//' to lamdat',maxph,
     1                z(locphm),0,' ')  
         write(iout,*) 'writing iosys file ', title(1:lentit)
         write(iout,*) 'words written = ',maxph
         locplm=locplm+wds
         locphm=locphm+maxph
   90 continue
c
c----------------------------------------------------------------------c   
c             begin energy loop                                        c
c----------------------------------------------------------------------c
      do 100 ien=1,nen
         k=sqrt(energy(ien))
c----------------------------------------------------------------------c
c             begin loop over shells, decomposition of                 c
c             inhomogeneity into spherical harmonics shell by shell.   c 
c----------------------------------------------------------------------c
         locyuk=yuk
         locflm=flm
         locr=rad
         locj=j
         locy=y           
         do 110 ns=1,nshell
c----------------------------------------------------------------------c
c             calculate bessel's of kr                                 c
c----------------------------------------------------------------------c 
            call vscale(z(locr+totr),z(locr),k,nrad(ns))
            call rcbes(z(locr+totr),z(locj),z(jp),z(locy),z(yp),
     1                 z(wron),scr,nrad(ns),lmax,ltop,
     2                 'derivatives',.false.)
c----------------------------------------------------------------------c
c             the next bit of monkey business converts derivatives     c
c             wrt k*r to derivatives wrt r, recomputes the wronskian   c
c             which is needed for the green's function and then        c       
c             scales the both functions so that the proper factor      c
c             will be incorporated into the radial integral.           c
c----------------------------------------------------------------------c  
            call sscal(ltop1,k,z(wron),1)
            call vinv(z(wron),z(wron),ltop1)
            call vsqrt(z(wron),z(wron),ltop1)
            call mvmul(z(wron),z(locj),z(locj),nrad(ns),ltop1)
            call mvmul(z(wron),z(locy),z(locy),nrad(ns),ltop1)
c----------------------------------------------------------------------c 
c             decompose the function into its angular components       c
c----------------------------------------------------------------------c      
            call lmdcmp(z(locyuk),z(plm),z(phifn),z(locr),z(locflm),
     1                  z(scr),mval,nl,nrad(ns),maxth,maxph,nm,prntlm)
            locr=locr+nrad(ns)
            locyuk=locyuk+nrad(ns)*maxth*maxph
            locflm=locflm+nrad(ns)*totl
            locj=locj+nrad(ns)*ltop1
            locy=locy+nrad(ns)*ltop1    
  110    continue
c----------------------------------------------------------------------c
c             from the radial components of the inhomogeneity          c
c             calculate the radial component of the solution           c
c             using the green's function to solve the integral         c
c             equation.                                                c
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c               do the forward integral                                c
c----------------------------------------------------------------------c
         call drvfrd(z(psilm),z(flm),z(j),z(y),z(wtrad),z(scr),nrad,
     1               mval,a(loclnd),nl,nm,ltop,totl,nshell,radqud)
  100 continue
      call iosys ('write integer maxsiz to rwf',1,canget,0,' ')
      call chainx(0)
      stop
    1 format (//,20x,'***** m6201: solve wave equation *****')
    2 format (/,5x,'the number of (l,m) pairs = ',i3,' their values are'
     1         ,/,(5x,5(:,'(',i3,',',i3,')')))
    3 format (/,5x,'number of radial shells = ',i4)
    4 format (/,5x,'number radial points per shell',(/,10(i4,1x)))
    5 format (/,5x,'shell radii',(/,5(1x,e15.8)))
    6 format (/,5x,'maximum l value',1x,i3,5x,'maximum m value',1x,i3)
    7 format (/,5x,'rms error for the ',i2,' region = ',e15.8)
    8 format (/,5x,'number of energies = ',i3,1x,'and their values',
     1              (/,5x,5(e15.8,1x)))
      end
 of (l,m) pairs = ',i3,' their values areinprog/yukawa.f                                                                                     0000644 0000311 0000024 00000003115 05411341424 0014003 0                                                                                                    ustar 00bis                             users                           0000000 0000000                                                                                                                                                                        *deck @(#)yukawa.f	1.1 9/9/91
c***begin prologue     yukawa
c***date written       930524   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           yukawa, link m6201
c***author             schneider, barry (nsf)
c***source             m6201
c***purpose            calculate multicenter yukawa potential
c***references         none
c                      storage of fun is packed.
c***routines called
c***end prologue        yukawa
      subroutine yukawa(fun,alpha,cen,r,theta,cphi,sphi,ncen,nr,
     1                  nthet,nphi)
      implicit integer (a-z)
      real*8 fun, alpha, cen, r, theta, cphi, sphi, snthe, temp
      real*8 x, y, z, d
      dimension fun(nr,nthet,nphi), r(nr), theta(nthet), cphi(nphi)
      dimension sphi(nphi), cen(3,ncen), alpha(ncen)
      common/io/ inp,iout
      call rzero(fun,nr*nthet*nphi)
      do 10 thept=1,nthet
         snthe=sqrt(1.d0-theta(thept)*theta(thept))
         do 20 rpt=1,nr
            temp=r(rpt)*snthe    
            z=r(rpt)*theta(thept)
            do 30 phipt=1,nphi
               x=temp*cphi(phipt)
               y=temp*sphi(phipt)      
               do 40 icen=1,ncen
                  d =sqrt( (x-cen(1,icen))*(x-cen(1,icen))+
     1                     (y-cen(2,icen))*(y-cen(2,icen)) + 
     2                     (z-cen(3,icen))*(z-cen(3,icen)) )
                  fun(rpt,thept,phipt) = fun(rpt,thept,phipt) + 
     1                                   exp(-alpha(icen)*d)/d
  40           continue
  30        continue
  20     continue                                           
  10  continue
      return
      end
=temp*sphi(phipt)      
               do 40 icen=1,ncen
                  d =sqrt( (x-cen(1,icen))*(x-cen(1,icen))+
     1                     (y-cen(2,icen))*(y-cen(2,icen)) + 
     2                     (z-cen(3,icen))*(z-cen(3,icen)) )
                  fun(rpt,thept,phipt) = fun(rpt,thept,phipt) + 
     1                                   exp(-alpha(icen)*d)/d
  40           continue
  30        continue
  20     continue                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     