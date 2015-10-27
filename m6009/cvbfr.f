*deck @(#)cvbfr.f	1.1 9/8/91
c***begin prologue     cvbfr
c***date written       880423   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           cvbfr, link 1106, kohn variational
c***author             schneider, barry (nsf)
c***source             m6009
c***purpose            copy of a real to real or complex array
c***references         
c***routines called    iosys, util and mdutil
c***end prologue       cvbfr
      subroutine cvbfr(a,b,c,n,m,type)
      implicit integer(a-z)
      real *8 a, b
      complex *16 c
      character *(*) type
      dimension a(n,m), b(n,m), c(n,m)
      common /io/ inp, iout
      if (type.eq.'real') then
          do 10 i=1,n
             do 20 j=1,m
                b(i,j)=a(i,j)
   20        continue 
   10     continue
      elseif (type.eq.'complex') then
          do 30 i=1,n
             do 40 j=1,m
                c(i,j)=a(i,j)
   40        continue 
   30     continue
      endif
      return
      end
