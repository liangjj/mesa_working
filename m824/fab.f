*deck @(#)fab.f	5.1  11/6/94
      subroutine fab(c,fxk,t1,ao,nbf,ncore,nbfnc,nnp,lag,nder,
     #               nocc)
c
      implicit integer (a-z)
c
      real*8 c(nbf,nbf),fxk(nbfnc),t1(nbf,nbf),ao(nnp)
      real*8 lag(nbf,nocc)
c
      call trtosq(t1,ao,nbf,nnp)
      call ebc(fxk,t1,c,nbf,nbf,ncore)
c
      call iosys('write real der_fab(x,k) to dints without rewinding',
     $            nbfnc,fxk,0,' ')
c
c     ----- add to the lagrangian -----
c
      do 190 i=1,nbfnc
         fxk(i)=fxk(i)*2.0d+00
  190 continue
c
      call vadd(lag,lag,fxk,nbfnc)
ctemp
cps         if (der.lt.9) go to 200
cps         write (6,23) der
cps   23    format (//,' ao lagrange after fab(xk)',i3,/)
cps         call matout(lag(1,1,der),nbf,nocc,nbf,nocc,6)
cend
c
c
      return
      end
