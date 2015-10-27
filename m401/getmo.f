*deck @(#)getmo.f	5.1  11/6/94
      subroutine getmo(smhalf,h,c,eigval,triang,u,t1,num,nnp)
c***begin prologue     getmo
c***date written       850601  yymmdd
c***revision date      yymmdd  yymmdd
c***keywords           diagonalization, eigenvalues, eigenvectors
c***author             martin, richard (lanl)
c***source             @(#)getmo.f	5.1   11/6/94
c***purpose            returns eigenvectors and eigenvalues of a real
c                      symmetric matrix.
c***description
c     call getmo(smhalf,h,c,eigval,triang,u,t1,num,nnp)
c
c***references         (none)
c***routines called    trtosq(math), ebc(math), ebtc(math), sqtotr(math),
c                      rsp(clams)
c***end prologue       getmo
      implicit integer(a-z)
      real*8 smhalf(nnp),h(nnp)
      real*8 c(num,num),eigval(num),triang(nnp),u(num,num),t1(num,num)
c
c     return the eigenvalues(eigval) and eigenvectors(c) of the matrix h.
c     h is symmetric(num,num) stored as nnp upper triangular elements,
c     and is left untouched by this routine.
c     smhalf is the symmetric orthogonalization matrix.
c     triang,u,t1 are scratch arrays.
c
c
c     transform to orthonormal basis.  h'=s(-1/2)*h*s(-1/2)
      call trtosq(t1,h,num,nnp)
      call trtosq(u,smhalf,num,nnp)
      call ebc(c,t1,u,num,num,num)
      call ebtc(t1,u,c,num,num,num)
c
c     diagonalize h.  h'u=eu
      call sqtotr(triang,t1,num,nnp)
      call degrsp(num,nnp,triang,eigval,1,u,t1,c)
c
c     back transform the eigenvectors.  c=s(-1/2)*u
      call trtosq(t1,smhalf,num,nnp)
      call ebtc(c,t1,u,num,num,num)
c
c
      return
      end
