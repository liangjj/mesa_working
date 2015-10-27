*deck @(#)dmat.f	5.1  11/6/94
      subroutine dmat (d,c,num,nnp,ndmat)
      implicit integer(a-z)
      dimension d(nnp,ndmat), c(num,num)
      real*8 fn,fd,an,ad,bn,bd
      common /shell/ noshel(10), orbshel(50,10), nvshel(10),
     1 vshel(200,10), f(10), aij(10,10), bij(10,10),
     2 fn(10),fd(10),an(10),ad(10),bn(10),bd(10)
      real*8 d, c, f, aij, bij
c
c*** driver subroutine to form density matrices for ivo program
c
      mink=1
      do 10 i=1,ndmat
      maxk=mink+noshel(i)-1
      call gdmat (d(1,i),c,num,num,nnp,mink,maxk)
c      mink=mink+noshel(i)+nvshel(i)
      mink=mink+noshel(i)
   10 continue
c
c     ----- write density matrix to rwf -----
c
c      call iosys ('write real density_matrix to rwf',nnp*ndmat,d,0,' ')
c
      return
      end
