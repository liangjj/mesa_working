*deck @(#)corges.f	5.1  11/6/94
      subroutine corges(smhalf,t,v,h,c,eigval,triang,u,t1,num,nnp,
     $                  altges)
c***begin prologue     corges
c***date written       850601  yymmdd
c***revision date      yymmdd  yymmdd
c***keywords           initial guess, core hamiltonian
c***author             martin, richard (lanl)
c***source             @(#)corges.f	5.1   11/6/94
c***purpose            forms and diagonalizes the core hamiltonian.
c***description
c     call corges(smhalf,t,v,h,c,eigval,triang,u,t1,
c                 num,nnp,altges)
c       smhalf  the s**(-1/2) matrix, i.e. the reciprocal square root
c               of the overlap matrix(nnp).
c       t       the kinetic energy matrix (nnp).
c       v       the potential energy matrix(nnp).
c       h       scratch matrix(nnp).
c       c       eigenvector matrix(num,num).
c       eigval  eigenvalue array(num)
c       triang  scratch array(nnp).
c       u       scratch matrix(num,num).
c       t1      scratch matrix(num,num).
c       num     basis set dimension.
c       nnp     num*(num+1)/2 .
c       altges  .true. if guess is to be altered.
c***references
c***routines called    iosys(io), vadd(math), getmo(m401), alter(m401)
c***end prologue       corges
      implicit integer(a-z)
      real*8 smhalf(nnp),t(nnp),v(nnp),h(nnp)
      real*8 c(num,num),eigval(num),triang(nnp),u(num,num),t1(num,num)
      logical altges
c
c     form and diagonalize the core hamiltonian.
c
      call iosys('read real "kinetic integrals" from rwf',nnp,t,0,' ')
      call iosys('read real "potential integrals" from rwf',nnp,v,0,' ')
c     form h=t+v.
      call vadd(h,t,v,nnp)
c
c     get the orbitals.
      call getmo(smhalf,h,c,eigval,triang,u,t1,num,nnp)
c
c     possibly alter the guess.
      if(altges) call alter(c,eigval,triang,num)
c
c
      return
      end
