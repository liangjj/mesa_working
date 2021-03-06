*deck %W%  %G%
      subroutine vmmul(v,matin,matout,n,m)
c***begin prologue     vmmul
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           vector, multiply
c***author             schneider, barry (nsf)
c***source             %W%   %G%
c***purpose            vectorized vector matrix multiply:  matout=v*matin.
c***description
c                      call vmmul(v,matin,matout,n,m)
c                        v        input vector of length n.
c                        matin    input matrix of size (n*m).
c                        matout   output matrix of size (n*m).
c                        n        length of vector ( matrix row size ).
c                        m        matrix column size. 
c
c***references
c***routines called    (none)
c***end prologue       vmmul
      real*8 v(n), matin(n,m), matout(n,m)
      do 1 i=1,m
         do 2 j=1,n
            matout(j,i)=v(j)*matin(j,i)
    2    continue
    1 continue            
      return
      end
