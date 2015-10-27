*deck @(#)cmadd.f	1.1  11/30/90
      subroutine cmadd(alpha,a,ia,beta,b,ib,c,ic,n,m)
c***begin prologue     cmadd
c***date written       960903  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           matrix, add
c***author             schneider, barry (nsf)
c***source             @(#)cmadd.f	1.1   9/3/96
c***purpose            vectorized matrix addition:  c=alpha*a+beta*b
c***description
c                        c        output matrix n*m.
c                        a        input matrix n*m.
c                        b        input matrix n*m.
c                        alpha    input scalar.
c                        beta     input scalar.
c
c***references
c***routines called    (none)
c***end prologue       cmadd
      complex*16 a(ia,m), b(ib,m), c(ic,m), alpha, beta
      do 1 i=1,n
         do 2 j=1,m
            c(i,j) = alpha*a(i,j)+beta*b(i,j)
    2    continue
    1 continue
      return
      end
