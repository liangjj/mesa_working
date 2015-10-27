*deck @(#)mvmul.f	1.2  10/27/94
      subroutine mvmul(v,matin,matout,n,m)
c***begin prologue     mvmul
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           vector, multiply
c***author             schneider, barry (nsf)
c***source             @(#)mvmul.f	1.2   10/27/94
c***purpose            vectorized matrix vector multiply:  matout=matin*v.
c***description
c                      call mvmul(v,matin,matout,n,m)
c                        v        input vector of length m.
c                        matin    input matrix of size (n*m).
c                        matout   output matrix of size (n*m).
c                        m        length of vector ( matrix column size ).
c                        n        matrix row size. 
c
c***references
c***routines called    (none)
c***end prologue       mvmul
      real*8 v(m), matin(n,m), matout(n,m)
      do 1 i=1,n
         do 2 j=1,m
            matout(i,j)=matin(i,j)*v(j)
    2    continue
    1 continue            
      return
      end
