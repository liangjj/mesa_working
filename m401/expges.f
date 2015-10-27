*deck @(#)expges.f	5.1  11/6/94
      subroutine expges(t,v,h,c,eigval,triang,u,t1,t2,
     $ num,nnp,nsmall)
c***begin prologue     expges
c***date written       890307  yymmdd
c***revision date      yymmdd  yymmdd
c***keywords           initial guess for e-scattering
c***author             lengsfield, byron (llnl)
c***source             @(#)expges.f	5.1   11/6/94
c***purpose            rotate orbitals under core-hamiltonian
c***description
c     call expges(t,v,h,c,eigval,triang,u,t1,
c                 num,nnp,nsmall)
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
c       nsmall  defines a partition of the orbitals space
c               the last num-nsmall orbitals are rotated
c
c***references
c***routines called    iosys(io), vadd(math), getmo(m401), alter(m401)
c***end prologue       expges
c
      implicit integer(a-z)
      real*8 t(nnp),v(nnp),h(nnp),t2(num,num)
      real*8 c(num,num),eigval(num),triang(nnp),u(num,num),t1(num,num)
c
c     form and diagonalize the core hamiltonian.
c
      call iosys('read real "kinetic integrals" from rwf',nnp,t,0,' ')
      call iosys('read real "potential integrals" from rwf',nnp,v,0,' ')
c     form h=t+v.
      call vadd(h,t,v,nnp)
c
      call trtosq(u,h,num,nnp)
c
      ix=num-nsmall+1
      call ebc(t1,u,c(1,ix),num,num,nsmall)
      call ebtc(u,c(1,ix),t1,nsmall,num,nsmall)
c
      nss=nsmall*(nsmall+1)/2
      call sqtotr(triang,u,nsmall,nss)
c
      call degrsp(nsmall,nss,triang,eigval(ix),1,u,t1,t2)
c
c     back transform the eigenvectors.  c=t1*u
c
      ntot=nsmall*num
      call scopy(ntot,c(1,ix),1,t1,1)
      call ebc(c(1,ix),t1,u,num,nsmall,nsmall)
c
c
      return
      end
