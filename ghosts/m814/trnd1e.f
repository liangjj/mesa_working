*deck @(#)trnd1e.f	1.1  11/30/90
      subroutine trnd1e(c,ao,fxj,fxb,mo,t1,nbf,nnp,ncore,nactiv,
     #                  nder,nnpact,nbfnc,nbfna,lag,dab,nocc)
c
      implicit integer (a-z)
c
      real*8 c(nbf,nbf),ao(nnp,nder),fxj(nbfnc,nder),fxb(nbfna,nder)
      real*8 mo(nnpact,nder),t1(nbf,nbf),lag(nbf,nocc,nder)
      real*8 dab(nactiv,nactiv)
c
      common /io/ inp,iout
c
c     ----- half-transform the derivative core fock over core orbitals
c
c
      if (ncore.gt.0) then
         do 200 der=1,nder
            call trtosq(t1,ao(1,der),nbf,nnp)
            call ebc(fxj(1,der),t1,c,nbf,nbf,ncore)
  200    continue
c
         call iosys('write real der_f(x,j) to ints',
     $               nbfnc*nder,fxj,0,' ')
c
c        ----- add to the lagrangians -----
c
         do 250 der=1,nder
            do 240 i=1,nbfnc
               fxj(i,der)=fxj(i,der)*2.0d+00
  240       continue
            call vadd(lag(1,1,der),lag(1,1,der),fxj(1,der),nbfnc)
ctemp
cps            if (der.lt.9) go to 250
cps            write (6,23) der
cps   23       format (//,' ao lagrang after fxj',i3,/)
cps            call matout(lag(1,1,der),nbf,nocc,nbf,nocc,6)
cend
  250    continue
c
      end if
c
c     ----- half-transform the derivative core fock over active orbitals
c
      do 300 der=1,nder
         call trtosq(t1,ao(1,der),nbf,nnp)
         call ebc(fxb(1,der),t1,c(1,ncore+1),nbf,nbf,nactiv)
  300 continue
c
      call iosys('write real der_f(x,b) to ints',nbfna*nder,fxb,0,' ')
c
c     ----- put this contribution in the lagrangian -----
c
      call iosys('read real mcscf_mo_1pdm from rwf',nnpact,t1,0,' ')
c
c     ----- the density is triangular, with off-diagonals doubled
c           since we need 2*f(x,b)*d(b,a), double diagonals and
c           square up, thus incorporating the factor of 2
c
      call trtosq(dab,t1,nactiv,nnpact)
c
      do 350 i=1,nactiv
         dab(i,i)=dab(i,i)*2.0d+00
  350 continue
cps      do 351 i=1,nactiv
cps         do 352 j=1,nactiv
cps            dab(i,j)=dab(i,j)*2.0d+00
cps  352    continue
cps  351 continue
c
      do 360 der=1,nder
         call apbc(lag(1,ncore+1,der),fxb(1,der),dab,nbf,nactiv,
     #                                    nactiv)
ctemp
cps         if (der.lt.9) go to 360
cps         write (6,24) der
cps   24    format (//,' ao lagrange after fxb',i3,/)
cps         call matout(lag(1,1,der),nbf,nocc,nbf,nocc,6)
cend
  360 continue
c
c     ----- complete the transformation to the active space -----
c
      do 400 der=1,nder
         call ebtc(t1,c(1,ncore+1),fxb(1,der),nactiv,nbf,nactiv)
         call sqtotr(mo(1,der),t1,nactiv,nnpact)
c..bhl
c         if(der.eq.1) then
c          write(iout,*)' ndf=1 der_h mo '
c          call matout(t1,nactiv,nactiv,nactiv,nactiv,iout)
c         endif
c..bhl
  400 continue
c
      call iosys('write real der_h to ints',nnpact*nder,mo,0,' ')
c
c
      return
      end
