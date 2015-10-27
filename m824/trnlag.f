*deck @(#)trnlag.f	5.1  11/6/94
      subroutine trnlag(c,lag,molag,nbf,nocc,nder)
c
      implicit integer (a-z)
c
      real*8 c(nbf,nbf),lag(*),molag(nbf,nocc)
c
      common /io/ inp,iout
c
c      real c(nbf,nbf),lag(nbf,nocc,nder),molag(nbf,nocc,nder)
c
c     ----- factor of two to get to the true lagrangian -----
c
      do 1 i=1,nbf*nocc
         lag(i)=lag(i)*2.0d+00
    1 continue
c
c     ----- transform the first index to mo basis -----
c
      call ebtc(molag,c,lag,nbf,nbf,nocc)
c
      call iosys('write real half_der_lagrangian to rwf'//
     #' without rewinding',nbf*nocc,lag,0,' ')
      call iosys('write real mo_der_lagrangian to rwf'//
     #' without rewinding',nbf*nocc,molag,0,' ')
c
cps      do 10 der=1,nder
c..bhl
c         der=1
c         write (iout,23) der
c  23     format (//,' mo_der lagrangian for ',i3,/)
c         call matout(molag(1,1,der),nbf,nocc,nbf,nocc,iout)
c..bhl
cps   10 continue
c
c
      return
      end
