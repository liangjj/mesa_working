*deck @(#)gencou.f	1.1 9/8/91
c***begin prologue     m7000
c***date written       920525   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m7000, link 7000, spline
c***author             schneider, b. (nsf)
c***source             m7004
c***purpose            regular and irregular coulomb functions on grid
c***
c***description        regular and irregular coulomb function are computed
c***                   by a combination of series expansion and numerical
c***                   integration and then spline fit. the spline
c***                   coefficients are used to calculate the functions and
c***                   their derivatives at arbitrary points using the 
c***                   pp represenations of the b-spline.
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       m7000
      subroutine gencou(a,b,ak,bk,a0,b0,c0,d0,e0,psi0,psi1,dpsi0,
     1                  dpsi1,charge,energy,v,x,xinv,c,scr,fs,xs,
     2                  rmin,rmax,rdel,rswtch,lval,nint,nl,toskp,
     3                  nspln,filnm,nospln,noireg,ntrms,level)
      implicit integer (a-z)
      real*8 a, b, ak, bk, psi0, psi1, dpsi0, dpsi1, charge
      real*8 energy, v, x, xinv, c, scr, rmin, rmax
      real*8 rdel, eta, cl, k2, k, clfun, dum, ddum
      real*8 regexp, fact, fs, xs, a0, b0, c0, d0, e0
      real*8 rswtch, zero, two, dlfun, dl, c0sqfn, c0sq, pre
      real*8 qlplfn, qlpl, wron, gamma, dsqrt2
      real*8 yp1r, ypnr, yp1i, ypni
      character *80 title
      character *(*) filnm
      logical nospln, noireg
      dimension a(*), b(*), ak(*), bk(*), dum(2), ddum(2)
      dimension a0(*), b0(*), c0(*), d0(*), e0(*)
      dimension psi0(nint), psi1(nint), dpsi0(nint), dpsi1(nint)
      dimension v(nint), x(nint), xinv(nint), c(nspln)
      dimension fs(nspln), xs(nspln), lval(nl), scr(nint,2)
      dimension fact(0:100)
      common /io/ inp, iout
      data zero, two, dsqrt2 / 0.d0, 2.d0, .70710678118654752440d0 /
      call rzero(scr,2*nint)
      call factl(fact,100)
c**********************************************************************c
c                    make integration and spline grid                  c
c**********************************************************************c
      call mkgrd(x,xinv,xs,rmin,rdel,nint,toskp)
      if (level.ge.3) then
          title='integration grid'
          write(iout,*) title
          write(iout,5) (x(i),i=1,nint)
      endif 
c**********************************************************************c
c                calculate k or kappa and the value of the             c
c                             coulomb eta                              c
c**********************************************************************c
      if(energy.lt.0.d0) then
         k2=abs(energy)*two
         k=sqrt(k2)
      else
         k=sqrt(energy*two)
      endif
      eta=charge/k
      write(iout,1)
      write(iout,2) energy, eta
c**********************************************************************c
c              find largest x for series expansion                     c
c**********************************************************************c
      do 500 pt=1,nint
         if (x(pt).gt.rswtch) go to 501
  500 continue
  501 ptbeg1=pt     
      ptbeg=ptbeg1-1
      npt=nint-ptbeg+1
      write(iout,7) ptbeg
c**********************************************************************c
c                  loop over l values needed                           c
c**********************************************************************c
      do 200 l=1,nl
         if (energy.ge.zero) then
c**********************************************************************c
c        calculate the coefficients for the series expansion           c
c        at small r and the asymptotic expansion at large r            c
c**********************************************************************c
             call abcoef(a,b,a0,b0,c0,d0,eta,lval(l),ntrms,
     1                   'positive',level)
             call lrcoef(ak,bk,e0,eta,lval(l),ntrms,'positive',level)
c**********************************************************************c
c                use series expansion for small x                      c
c**********************************************************************c
             write(iout,*)
             write(iout,*) '     calculating regular function'
             cl=clfun(eta,lval(l),fact)
             call regexp(psi0,dpsi0,x,a,lval(l),cl,ptbeg1,ntrms,level)
             yp1r=dpsi0(1)
             if (level.ge.2) then
                 write(iout,*)
                 title='series expansion of regular solution'
                 write(iout,*) title
                 write(iout,5) (psi0(i),i=1,ptbeg1,toskp)
             endif
c**********************************************************************c
c              calculate effective potential for this l                c
c**********************************************************************c
             call potntl(v(ptbeg),eta,lval(l),xinv(ptbeg),npt,
     1                   'oscillatory')
c**********************************************************************c
c        integrate the homogeneous equation for the regular            c
c                              and                                     c
c                      irregular coulomb function                      c
c**********************************************************************c
c**********************************************************************c
c                  integrate forward for regular                       c
c                  solution to nint using psi0(ptbeg) and              c
c                  psi0(ptbeg1) as the starting values                 c
c**********************************************************************c
             call coulnt(psi0(ptbeg),v(ptbeg),scr(ptbeg,1),rdel,
     1                   npt,'forward')
             if (level.ge.1) then
                 write(iout,*)
                 title='regular solution'
                 write(iout,*) title
                 write(iout,5) (psi0(i),i=1,nint,toskp)
             endif
             call extphs(psi0(nint),x(nint),lval(l),eta,ak,bk,ntrms,
     1                   level)
             if (.not.noireg) then
                  write(iout,*)
                  write(iout,*) '     calculating irregular function'
c**********************************************************************c
c            get two values of irregular function using asymptotic     c
c                   expansion to start off integration                 c
c**********************************************************************c
                  call asymp(dum,psi1(nint-1),ddum,dpsi1(nint-1),ak,bk,
     1                       x(nint-1),xinv(nint-1),eta,lval(l),2,
     2                       ntrms,wron,level)
c**********************************************************************c
c                 integrate backward to the point ptbeg                c
c**********************************************************************c
                  call coulnt(psi1(ptbeg),v(ptbeg),scr(ptbeg,1),rdel,
     1                        npt,'backward')
c**********************************************************************c
c                finish off using the series from the first point      c
c                             to ptbeg                                 c
c**********************************************************************c
                  dl=dlfun(eta,lval(l),fact)
                  c0sq=c0sqfn(eta)
                  pre=two*eta/c0sq
                  qlpl=qlplfn(eta,lval(l),fact)
                  if (level.ge.2) then
                      write(iout,8) dl, c0sq, pre, qlpl                      
                  endif
                  call iregxp(psi1,dpsi1,psi0,dpsi0,x,b,lval(l),dl,
     1                        pre,qlpl,ptbeg,ntrms,level) 
                  yp1i=dpsi1(1)
                  if (level.ge.2) then
                      write(iout,*)
                      title='series expansion of irregular solution'
                      write(iout,*) title
                      write(iout,5) (psi1(i),i=1,ptbeg,toskp)
                  endif
                  if (level.ge.1) then
                      write(iout,*)
                      title='irregular solution'
                      write(iout,*) title 
                      write(iout,5) (psi1(i),i=1,nint,toskp)
                  endif
             endif
         else
c**********************************************************************c
c        calculate the coefficients for the series expansion           c
c        at small r and the asymptotic expansion at large r            c
c**********************************************************************c
             call abcoef(a,b,a0,b0,c0,d0,eta,lval(l),ntrms,
     1                   'negative',level)
             call lrcoef(ak,bk,e0,eta,lval(l),ntrms,'negative',level)
c**********************************************************************c
c              calculate effective potential for this l                c
c**********************************************************************c
             call potntl(v(ptbeg),eta,lval(l),xinv(ptbeg),npt,
     1                   'exponential')
c**********************************************************************c
c                use series expansion for small x                      c
c**********************************************************************c
             write(iout,*)
             write(iout,*) '     calculating regular function'
             call regexn(psi0,dpsi0,x,xinv,a0,lval(l),ptbeg1,ntrms,
     1                   level)
             yp1r=dpsi0(1)
             if (level.ge.2) then
                 write(iout,*)
                 title='series expansion of regular solution'
                 write(iout,*) title
                 write(iout,5) (psi0(i),i=1,ptbeg1,toskp)
             endif
c**********************************************************************c
c                  integrate forward for regular                       c
c                  solution to nint using psi0(ptbeg) and              c
c                  psi0(ptbeg1) as the starting values                 c
c**********************************************************************c
             call coulnt(psi0(ptbeg),v(ptbeg),scr(ptbeg,1),rdel,
     1                   npt,'forward')
             if (level.ge.1) then
                 write(iout,*)
                 title='regular solution'
                 write(iout,*) title
                 write(iout,5) (psi0(i),i=1,nint,toskp)
             endif
             if (.not.noireg) then
                  write(iout,*)
                  write(iout,*) '     calculating irregular function'
c**********************************************************************c
c            get two values of irregular function using asymptotic     c
c                   expansion to start off integration                 c
c**********************************************************************c
                  call asymn(psi1(nint-1),dpsi1(nint-1),e0,x(nint-1),
     1                       xinv(nint-1),eta,lval(l),2,ntrms,level)
c**********************************************************************c
c                 integrate backward to the point ptbeg                c
c**********************************************************************c
                  call coulnt(psi1(ptbeg),v(ptbeg),scr(ptbeg,1),rdel,
     1                        npt,'backward')
c**********************************************************************c
c                finish off using the series from the first point      c
c                             to ptbeg                                 c
c**********************************************************************c
                  pre=dsqrt2/(gamma(eta+lval(l)+1)*gamma(eta-lval(l)))
                  call iregxn(psi1,dpsi1,a0,b0,c0,pre,x,xinv,lval(l),
     1                        ptbeg,ntrms,level) 
                  yp1i=dpsi1(1)
                  if (level.ge.2) then
                      write(iout,*)
                      title='series expansion of irregular solution'
                      write(iout,*) title
                      write(iout,5) (psi1(i),i=1,ptbeg,toskp)
                  endif
                  if (level.ge.1) then
                      write(iout,*)
                      title='irregular solution'
                      write(iout,*) title 
                      write(iout,5) (psi1(i),i=1,nint,toskp)
                  endif
             endif
         endif 
c**********************************************************************c
c              make the spline coefficients                            c
c**********************************************************************c
         if (.not.nospln) then
             count=0
             do 300 pt=1,nint,toskp
                count=count+1
                fs(count)=psi0(pt)
  300        continue
             call iosys ('write real '//filnm//' to atomci without '//
     1                   'rewinding',nspln,fs,0,' ')
             ypnr=(fs(count)-fs(count-1))/(xs(count)-xs(count-1))
             call splinr(xs,fs,yp1r,ypnr,scr(1,2),c,nspln)
             if (level.ge.2) then
                 write(iout,9) yp1r,ypnr
                 write(iout,10) (c(pt),pt=1,nspln)
             endif
             call iosys ('write real '//filnm//' to atomci without '//
     1                   'rewinding',nspln,c,0,' ')
             if (.not.noireg) then
                 count=0
                 do 400 pt=1,nint,toskp
                    count=count+1
                    fs(count)=psi1(pt)
  400            continue
                 call iosys ('write real '//filnm//' to atomci '//
     1                       'without rewinding',nspln,fs,0,' ')
                 ypni=(fs(count)-fs(count-1))/(xs(count)-xs(count-1))
                 call splinr(xs,fs,yp1i,ypni,scr(1,2),c,nspln)
                 if (level.ge.2) then
                     write(iout,9) yp1i,ypni
                     write(iout,10) (c(pt),pt=1,nspln)
                 endif
                 call iosys ('write real '//filnm//' to atomci '//
     1                       'without rewinding',nspln,c,0,' ')
             endif
         endif
  200 continue
      return
    1 format(//,15x,'integration of coulomb equation')   
    2 format(/,5x,'energy',1x,e15.8,1x,'eta',1x,e15.8)   
    3 format(/,5x,'l value',1x,i2,1x,'cl coefficient',1x,e15.8,/,5x,
     1            'starting values of solution',1x,e15.8,1x,e15.8)
    5 format(5x,5e15.8)   
    7 format(/,5x,'integrate differential equation from point',1x,i5,
     1            1x,'onward')
    8 format(/,5x,'dl = ',e15.8,5x,'c0**2 = ',e15.8,/,5x,
     1            'two*eta/c0**2 = ',e15.8,1x,'ql/pl = ',e15.8)    
    9 format(/,5x,'interpolated first derivatives for spline fit',/,
     1         10x,'df(1) = ',e15.8,1x,'df(n) = ',e15.8)
   10 format(/,5x,'spline coefficients',/,(5x,5e15.8x))
      end
