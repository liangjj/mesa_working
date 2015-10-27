*deck @(#)trnds.f	1.1  11/30/90
      subroutine trnds(c,ao,mo,t1,t2,nbf,nnp,nder,nbfsq)
c
      implicit integer (a-z)
c
      real*8 c(nbf,nbf),ao(nnp,nder),mo(nbfsq,nder)
      real*8 t1(nbf,nbf),t2(nbf,nbf)
c
c     ----- read in and transform the derivative overlap integrals -----
c
      call iosys('read real "ao derivative overlap integrals"'//
     #' from ints',nnp*nder,ao,0,' ')
c
      do 100 der=1,nder
         call trtosq(t1,ao(1,der),nbf,nnp)
         call ebc(t2,t1,c,nbf,nbf,nbf)
         call ebtc(mo(1,der),c,t2,nbf,nbf,nbf)
c.nonsense
c.bhl         call ebc(t2,t1,c(1,ncore+1),nbf,nbf,nbf)
c.bhl         call ebtc(mo(1,der),c(1,ncore+1),t2,nbf,nbf,nbf)
c.nonsense
  100 continue
c
      call iosys('write real mo_der_overlap to ints',nbf**2*nder,
     #            mo,0,' ')
c
c
      return
      end
