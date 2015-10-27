*deck @(#)fab.f	1.1  11/30/90
      subroutine fab(c,fxk,t1,ao,nbf,ncore,nbfnc,nnp,lag,nder,
     #               nocc)
c
      implicit integer (a-z)
c
      real*8 c(nbf,nbf),fxk(nbfnc,nder),t1(nbf,nbf),ao(nnp,nder)
      real*8 lag(nbf,nocc,nder)
c
      do 100 der=1,nder
         call trtosq(t1,ao(1,der),nbf,nnp)
         call ebc(fxk(1,der),t1,c,nbf,nbf,ncore)
  100 continue
c
      call iosys('write real der_fab(x,k) to ints',
     $            nbfnc*nder,fxk,0,' ')
c
c     ----- add to the lagrangian -----
c
      do 200 der=1,nder
         do 190 i=1,nbfnc
            fxk(i,der)=fxk(i,der)*2.0d+00
  190    continue
         call vadd(lag(1,1,der),lag(1,1,der),fxk(1,der),nbfnc)
ctemp
cps         if (der.lt.9) go to 200
cps         write (6,23) der
cps   23    format (//,' ao lagrange after fab(xk)',i3,/)
cps         call matout(lag(1,1,der),nbf,nocc,nbf,nocc,6)
cend
  200 continue
c
c
      return
      end
