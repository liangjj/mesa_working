*deck @(#)flip.f	5.1  11/6/94
      subroutine flip (c,num,homo,occmo)
      implicit integer(a-z)
      real*8 c, val
      dimension c(num,num)
c
      do 1 i=1,num
      val=c(i,homo)
      c(i,homo)=c(i,occmo)
      c(i,occmo)=val
 1    continue
c
      return
      end
