*deck @(#)trnds.f	5.1  11/6/94
      subroutine trnds(c,ao,mo,t1,t2,nbf,nnp,nder,nbfsq)
c
      implicit integer (a-z)
c
      real*8 c(nbf,nbf),ao(nnp),mo(nbfsq)
      real*8 t1(nbf,nbf),t2(nbf,nbf)
c
c     ----- read in and transform the derivative overlap integrals -----
c
      call iosys('read real "ao derivative overlap integrals"'//
     #           ' from rdints',nnp,ao,(nder-1)*nnp,' ')
c
      call iosys('write real "ao derivative overlap integrals"'//
     #           ' to dints without rewinding',nnp,ao,0,' ')
c
         call trtosq(t1,ao,nbf,nnp)
         call ebc(t2,t1,c,nbf,nbf,nbf)
         call ebtc(mo,c,t2,nbf,nbf,nbf)
c
      call iosys('write real mo_der_overlap to dints'//
     #           ' without rewinding',nbf**2,mo,0,' ')
c
c
      return
      end
