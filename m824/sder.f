*deck @(#)sder.f	5.1  11/6/94
      subroutine sder(s,lag,num,nnp,t1,t2,grad,der)
c
c***begin prologue     sder
c***date written       861209  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c
c
c***keywords
c***author             saxe, paul (lanl)
c***source             @(#)sder.f	5.1   11/6/94
c***purpose            formation of the overlap contribution to the
c                      scf gradients from lagrangian and derivative
c                      overlap integrals.
c***description
c
c***references
c***routines called    (none)
c***end prologue       sder
c
      implicit integer (a-z)
c
      real*8 s(nnp),lag(nnp),t1(num,num),t2(num,num),grad(*)
      real*8 sdot
c
c     ----- read in the ao scf lagrangian and square it up -----
c
      call iosys('read real "scf ao lagrangian" from rwf',nnp,lag,0,' ')
      call trtosq(t1,lag,num,nnp)
c
c     ----- loop through derivative ao overlap integrals, processing
c
         call iosys('read real "ao derivative overlap integrals" '//
     $        'from rdints',nnp,s,(der-1)*nnp,' ')
         call iosys('write real "ao derivative overlap integrals" '//
     $              'to dints without rewinding',nnp,s,0,' ')
         call trtosq(t2,s,num,nnp)
c
         grad(der)=grad(der)-2.0d+00*sdot(num**2,t1,1,t2,1)
c
c
      return
      end
