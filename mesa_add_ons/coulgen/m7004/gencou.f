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
c***                   pp representations of the b-spline.
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       m7000
      subroutine gencou(a,b,ak,bk,a0,b0,c0,d0,e0,psi0,psi1,dpsi0,
     1                  dpsi1,charge,energy,v,x,xinv,break,c,sc,
     2                  scr,fs,xs,rhs,sdrv,rmin,rmax,rdel,rswtch,
     3                  ptbeg,last,lval,nint,nl,order,nbreak,toskp,
     4                  nspln,filnm,nospln,noireg,ntrms,drven,
     5                  prser,prfun)
      implicit integer (a-z)
      real*8 a, b, ak, bk, psi0, psi1, dpsi0, dpsi1, charge
      real*8 energy, v, x, xinv, break, c, sc, scr, rmin, rmax
      real*8 rdel, eta, cl, k2, k, clfun, dum, ddum
      real*8 regexp, fact, fs, xs, a0, b0, c0, d0, e0
      real*8 rswtch, zero, two, dlfun, dl, c0sqfn, c0sq, pre
      real*8 qlplfn, qlpl, wron, gamma, dsqrt2, rhs, sdrv
      character *80 title
      character *16 fptoc, enchr
      character *2 itoc
      character *18 label
      character *(*) filnm
      logical nospln, noireg, prser, prfun, drven
      dimension a(*), b(*), ak(*), bk(*), dum(4), ddum(2)
      dimension a0(*), b0(*), c0(*), d0(*), e0(*)
      dimension psi0(nint), psi1(nint), dpsi0(nint), dpsi1(nint)
      dimension v(nint), x(nint), xinv(nint), c(order,nbreak)
      dimension sc((4*nspln+1)*order,2), break(nbreak+1)
      dimension fs(nspln), xs(nspln), lval(nl), scr(nint), fact(0:100)
      dimension rhs(nint), sdrv(nint)
      common /io/ inp, iout
      data zero, two, dsqrt2 / 0.d0, 2.d0, .70710678118654752440d0 /
      npt=nint-ptbeg+1
      ptbeg1=ptbeg+1
      enchr=fptoc(energy)
      call factl(fact,100)
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
      if (.not.nospln) then
          call copy(sc(1,1),sc(1,2),(4*nspln+1)*order)
      endif
c**********************************************************************c
c                  loop over l values needed                           c
c**********************************************************************c
      do 200 l=1,nl
         call rzero(scr,nint)
         label=enchr//'-'//itoc(lval(l))
         if (energy.ge.zero) then
c**********************************************************************c
c                     positive energy branch                           c
c**********************************************************************c
c**********************************************************************c
c        calculate the coefficients for the series expansion           c
c        at small r and the asymptotic expansion at large r            c
c**********************************************************************c
             call abcoef(a,b,a0,b0,c0,d0,eta,lval(l),ntrms,
     1                   'positive',prser)
             if (.not.nospln) then
                 filnm='"'//label//'-a"' 
                 call iosys ('write real '//filnm//' to atomci',ntrms,
     1                        a,0,' ')             
                 filnm='"'//label//'-b"' 
                 call iosys ('write real '//filnm//' to atomci',ntrms,
     1                        b,0,' ')             
             endif
             call lrcoef(ak,bk,e0,eta,lval(l),ntrms,'positive',prser)
c**********************************************************************c
c                use series expansion for small x                      c
c**********************************************************************c
             write(iout,*)
             write(iout,*) '     calculating regular function'
             cl=clfun(eta,lval(l),fact)
             filnm='"'//label//'-cl"'
             call iosys ('write real '//filnm//' to atomci',1,cl,0,' ')
             call regexp(psi0,dpsi0,x,a,lval(l),cl,ptbeg1,ntrms,prfun)
             if (prfun) then
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
             call coulnt(psi0(ptbeg),v(ptbeg),scr(ptbeg),rdel,
     1                   npt,'forward')
             if (prfun) then
                 write(iout,*)
                 title='regular solution'
                 write(iout,*) title
                 write(iout,5) (psi0(i),i=1,nint,toskp)
             endif
             call extphs(psi0(nint),x(nint),lval(l),eta,ak,bk,ntrms,
     1                   prfun)
             if (.not.noireg) then
                  write(iout,*)
                  write(iout,*) '     calculating irregular function'
c**********************************************************************c
c            get two values of irregular function using asymptotic     c
c                   expansion to start off integration                 c
c**********************************************************************c
                  call asymp(dum,psi1(nint-1),ddum,dpsi1(nint-1),ak,bk,
     1                       x(nint-1),xinv(nint-1),eta,lval(l),2,
     2                       ntrms,wron,prfun)
c**********************************************************************c
c                 integrate backward to the point ptbeg                c
c**********************************************************************c
                  call coulnt(psi1(ptbeg),v(ptbeg),scr(ptbeg),rdel,
     1                        npt,'backward')
c**********************************************************************c
c                finish off using the series from the first point      c
c                             to ptbeg                                 c
c**********************************************************************c
                  dl=dlfun(eta,lval(l),fact)
                  c0sq=c0sqfn(eta)
                  pre=two*eta/c0sq
                  qlpl=qlplfn(eta,lval(l),fact)
                  dum(1)=dl
                  dum(2)=c0sq
                  dum(3)=pre
                  dum(4)=qlpl
                  if (.not.nospln) then
                      filnm='"'//label//'-factors"'
                      call iosys ('write real '//filnm//' to atomci',
     1                             4,dum,0,' ')             
                  endif
                  write(iout,7) dl, c0sq, pre, qlpl                      
                  call iregxp(psi1,dpsi1,psi0,dpsi0,x,b,lval(l),dl,
     1                        pre,qlpl,ptbeg,ntrms,prfun) 
                  if (prfun) then
                      write(iout,*)
                      title='series expansion of irregular solution'
                      write(iout,*) title
                      write(iout,5) (psi1(i),i=1,ptbeg,toskp)
                  endif
                  if (prfun) then
                      write(iout,*)
                      title='irregular solution'
                      write(iout,*) title 
                      write(iout,5) (psi1(i),i=1,nint,toskp)
                  endif
             endif
         else
c**********************************************************************c
c                       negative energy branch                         c
c**********************************************************************c
c**********************************************************************c
c        calculate the coefficients for the series expansion           c
c        at small r and the asymptotic expansion at large r            c
c**********************************************************************c
             call abcoef(a,b,a0,b0,c0,d0,eta,lval(l),ntrms,
     1                   'negative',prser)
             if (.not.nospln) then
                 filnm='"'//label//'-a0"' 
                 call iosys ('write real '//filnm//' to atomci',ntrms,
     1                        a0,0,' ')             
                 filnm='"'//label//'-b0"' 
                 call iosys ('write real '//filnm//' to atomci',ntrms,
     1                        b0,0,' ')             
                 filnm='"'//label//'-c0"' 
                 call iosys ('write real '//filnm//' to atomci',ntrms,
     1                        c0,0,' ')
             endif
             call lrcoef(ak,bk,e0,eta,lval(l),ntrms,'negative',prser)
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
     1                   prfun)
             if (prfun) then
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
             call coulnt(psi0(ptbeg),v(ptbeg),scr(ptbeg),rdel,
     1                   npt,'forward')
             if (prfun) then
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
     1                       xinv(nint-1),eta,lval(l),2,ntrms,prfun)
c**********************************************************************c
c                 integrate backward to the point ptbeg                c
c**********************************************************************c
                  call coulnt(psi1(ptbeg),v(ptbeg),scr(ptbeg),rdel,
     1                        npt,'backward')
c**********************************************************************c
c                finish off using the series from the first point      c
c                             to ptbeg                                 c
c**********************************************************************c
                  pre=dsqrt2/(gamma(eta+lval(l)+1)*gamma(eta-lval(l)))
                  if (.not.nospln) then
                      filnm='"'//label//'-factors"'
                      call iosys ('write real '//filnm//' to atomci',
     1                             1,pre,0,' ')
                  endif
                  call iregxn(psi1,dpsi1,a0,b0,c0,pre,x,xinv,lval(l),
     1                        ptbeg,ntrms,prfun) 
                  if (prfun) then
                      write(iout,*)
                      title='series expansion of irregular solution'
                      write(iout,*) title
                      write(iout,5) (psi1(i),i=1,ptbeg,toskp)
                  endif
                  if (prfun) then
                      write(iout,*)
                      title='irregular solution'
                      write(iout,*) title 
                      write(iout,5) (psi1(i),i=1,nint,toskp)
                  endif
             endif
         endif 
c**********************************************************************c
c                solve inhomogeneous equation if desired               c
c**********************************************************************c
         if (drven) then
             call drver(scr,rhs,psi0,nint)
             call inhomo(scr,psi0,psi1,sdrv,dpsi0,dpsi1,rdel,
     1                   dum(1),dum(2),last,nint)
             filnm='"'//label//'-drints"' 
             call iosys ('write real '//filnm//' to atomci',2,
     1                    dum,0,' ')
             if (prfun) then
                 write(iout,*)
                 title='driven solution'
                 write(iout,*) title 
                 write(iout,5) (sdrv(i),i=1,nint,toskp)
             endif    
         endif
c**********************************************************************c
c              make the spline coefficients                            c
c**********************************************************************c
         if (.not.nospln) then
             count=0
             do 300 pt=ptbeg,nint,toskp
                count=count+1
                fs(count)=psi0(pt)
  300        continue        
             call copy(sc(1,2),sc(1,1),(4*nspln+1)*order) 
             call splcof(nspln,fs,order,nbreak,break,c,sc)
             filnm='"'//label//'-c reg"' 
             call iosys ('write real '//filnm//' to atomci',
     1                    nbreak*order,c,0,' ')
             if (.not.noireg) then
                 count=0
                 do 400 pt=ptbeg,nint,toskp
                    count=count+1
                    fs(count)=psi1(pt)
  400            continue        
                 call copy(sc(1,2),sc(1,1),(4*nspln+1)*order) 
                 call splcof(nspln,fs,order,nbreak,break,c,sc)
                 filnm='"'//label//'-c ireg"' 
                 call iosys ('write real '//filnm//' to atomci',
     1                        nbreak*order,c,0,' ')
             endif
             if (drven) then
                 count=0
                 do 500 pt=ptbeg,nint,toskp
                    count=count+1
                    fs(count)=sdrv(pt)
  500            continue        
                 call copy(sc(1,2),sc(1,1),(4*nspln+1)*order) 
                 call splcof(nspln,fs,order,nbreak,break,c,sc)
                 filnm='"'//label//'-c driven"' 
                 call iosys ('write real '//filnm//' to atomci',
     1                        nbreak*order,c,0,' ')
             endif
         endif
  200 continue
      return
    1 format(//,15x,'integration of coulomb equation')   
    2 format(/,5x,'energy',1x,e15.8,1x,'eta',1x,e15.8)   
    3 format(/,5x,'l value',1x,i2,1x,'cl coefficient',1x,e15.8,/,5x,
     1            'starting values of solution',1x,e15.8,1x,e15.8)
    5 format(5x,5e15.8)   
    7 format(/,5x,'dl = ',e15.8,5x,'c0**2 = ',e15.8,/,5x,
     1            'two*eta/c0**2 = ',e15.8,1x,'ql/pl = ',e15.8)    
      end













