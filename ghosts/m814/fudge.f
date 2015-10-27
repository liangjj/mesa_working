*deck @(#)fudge.f	1.1  11/30/90
      subroutine fudge(ao,nnp,nder)
c
      real*8 ao(*)
c
      call iosys('read real "kinetic integrals" from rwf',nnp,ao,0,' ')
      call iosys('read real "potential integrals" from rwf',nnp,
     #            ao(nnp+1),0,' ')
      do 1 i=1,nnp
         ao(i)=ao(i)+ao(i+nnp)
    1 continue
      do 2 i=2,nder
         ioff=(i-1)*nnp
         call vmove(ao(ioff+1),ao,nnp)
    2 continue
      call iosys('write real "ao derivative one-electron integrals"'//
     #           ' to ints',nnp*nder,ao,0,' ')
c
      call iosys('read real "sorted ao integrals" from ints',nnp**2,
     #            ao,0,' ')
      do 3 i=2,nder
         ioff=(i-1)*nnp**2
         call vmove(ao(ioff+1),ao,nnp**2)
    3 continue
c
      call iosys('write real "sorted ao derivative integrals"'//
     #           ' to ints',nnp**2*nder,ao,0,' ')
c
      return
      end
