*deck @(#)vimmul.f	1.1  11/30/90
      subroutine vimmul(v,matin,matout,n,m)
c***begin prologue     vimmul
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           vector, multiply
c***author             schneider, barry (nsf)
c***source             @(#)vimmul.f	1.1   6/22/93              -1
c***purpose            vectorized vector matrix multiply:  matout=v *matin.
c***description
c                      call vimmul(v,matin,matout,n,m)
c                        v        input vector of length n.
c                        matin    input matrix of size (n*m).
c                        matout   output matrix of size (n*m).
c                        n        length of vector ( matrix row size ).
c                        m        matrix column size. 
c
c***references
c***routines called    (none)
c***end prologue       vimmul
      real*8 v(n), matin(n,m), matout(n,m)
      do 1 i=1,m
         do 2 j=1,n
            matout(j,i)=matin(j,i)/v(j)
    2    continue
    1 continue            
      return
      end
