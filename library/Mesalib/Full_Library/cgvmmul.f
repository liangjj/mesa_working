*deck @(#)cgvmmul.f	1.1  11/30/90
      subroutine cgvmmul(v,a,b,c,n,m)
c***begin prologue     cgvmmul
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           vector, multiply
c***author             schneider, barry (nsf)
c***source             @(#)gvmmul.f	1.1   6/22/93
c***purpose            vectorized vector matrix multiply and add:
c***                                                 c = v*a + b
c***description
c                      call cgvmmul(v,a,b,c,n,m)
c                        v    input vector of length n.
c                        a    input matrix of size (n*m).
c                        b    input matrix of size (n*m).
c                        c    output matrix of size (n*m)
c                        n    length of vector ( matrix row size ).
c                        m    matrix column size. 
c
c***references
c***routines called    (none)
c***end prologue       gvmmul
      complex*16 v(n), a(n,m), b(n,m), c(n,m)
      do 1 i=1,m
         do 2 j=1,n
            c(j,i) = v(j)*a(j,i) + b(j,i)
    2    continue
    1 continue            
      return
      end









