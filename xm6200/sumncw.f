*deck %W%  %G%
c***begin prologue     sumncw
c***date written       930802   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           sumncw, link m6200
c***author             schneider, barry (nsf)
c***source             m6200
c***purpose            calculate newton cotes weight for entire
c***                   subdomain as sum over the n-1 partial weights.
c
c***routines called
c***end prologue       sumncw
      subroutine sumncw(wt,n)
      implicit integer (a-z)
      real*8 wt, sum
      dimension wt(n,n-1)
      do 10 i=1,n
         sum=0.d0
         do 20 j=1,n-1
            sum=sum+wt(i,j)
   20    continue
         wt(i,1)=sum
   10 continue
      return
      end

