      subroutine prrmat(r1,r2,r3,r4,s1,s2,n)
      implicit integer (a-z)
      common /io/ inp, iout
      real *8 s1, s2, r1, r2, r3, r4
      dimension s1(n), s2(n)
      r1=s2(n)
      r2=s1(n)
      r3=s1(1)
      r4=s2(1)
      write (iout,10) r1, r2, r3, r4
      return
   10 format (/,'the four outer solution r matrices',/,10x,
     1        4e15.8)
      end
