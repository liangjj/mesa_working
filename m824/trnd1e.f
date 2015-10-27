*deck @(#)trnd1e.f	5.1  11/6/94
      subroutine trnd1e(c,ao,fxj,fxb,mo,t1,nbf,nnp,ncore,nactiv,
     #                  nder,nnpact,nbfnc,nbfna,lag,dab,nocc)
c
      implicit integer (a-z)
c
      real*8 c(nbf,nbf),ao(nnp),fxj(nbfnc),fxb(nbfna)
      real*8 mo(nnpact),t1(nbf,nbf),lag(nbf,nocc)
      real*8 dab(nactiv,nactiv)
c
      common /io/ inp,iout
c
c     ----- half-transform the derivative core fock over core orbitals
c
c
      if (ncore.gt.0) then
            call trtosq(t1,ao,nbf,nnp)
            call ebc(fxj,t1,c,nbf,nbf,ncore)
c
         call iosys('write real der_f(x,j) to dints without rewinding',
     $     nbfnc,fxj,0,' ')
c
c        ----- add to the lagrangians -----
c
            do 240 i=1,nbfnc
               fxj(i)=fxj(i)*2.0d+00
  240       continue
            call vadd(lag,lag,fxj,nbfnc)
ctemp
cps            if (der.lt.9) go to 250
cps            write (6,23) der
cps   23       format (//,' ao lagrang after fxj',i3,/)
cps            call matout(lag(1,1,der),nbf,nocc,nbf,nocc,6)
cend
c
      end if
c
c     ----- half-transform the derivative core fock over active orbitals
c
         call trtosq(t1,ao,nbf,nnp)
         call ebc(fxb,t1,c(1,ncore+1),nbf,nbf,nactiv)
c
      call iosys('write real der_f(x,b) to dints  without rewinding',
     $   nbfna,fxb,0,' ')
c
c     ----- put this contribution in the lagrangian -----
c
      call iosys('read real mcscf_mo_1pdm from rwf',nnpact,t1,
     # 0,' ')
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
c
         call apbc(lag(1,ncore+1),fxb,dab,nbf,nactiv,nactiv)
ctemp
cps         if (der.lt.9) go to 360
cps         write (6,24) der
cps   24    format (//,' ao lagrange after fxb',i3,/)
cps         call matout(lag(1,1,der),nbf,nocc,nbf,nocc,6)
cend
c
c     ----- complete the transformation to the active space -----
c
         call ebtc(t1,c(1,ncore+1),fxb,nactiv,nbf,nactiv)
         call sqtotr(mo,t1,nactiv,nnpact)
c..bhl
c         if(der.eq.1) then
c          write(iout,*)' ndf=1 der_h mo '
c          call matout(t1,nactiv,nactiv,nactiv,nactiv,iout)
c         endif
c..bhl
c
      call iosys('write real der_h to dints without rewinding',
     $ nnpact,mo,0,' ')
c
c
      return
      end
