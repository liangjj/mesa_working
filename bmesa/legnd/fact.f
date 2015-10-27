*deck fact.f	1.1 9/8/91
c***begin prologue     fact
c***date written       880721   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           fact, link 2702, factorials
c***author             schneider, barry (lanl)
c***source             m2702
c***purpose            factorials
c***description        calculation of factorials
c***references         none
c
c***routines called
c***end prologue       fact
      subroutine fact(dfct,ddfct,maxfac)
      implicit integer (a-z)
      real *8 dfct, ddfct
      dimension dfct(0:maxfac), ddfct(0:maxfac)
c----------------------------------------------------------------------c
c               calculate factorials                                   c
c----------------------------------------------------------------------c
      dfct(0)=1.d+00
      dfct(1)=1.d+00
      if (maxfac.gt.1) then
          do 10 i=2,maxfac
             dfct(i)=i*dfct(i-1)
   10     continue
      endif
c----------------------------------------------------------------------c
c           calculate (2*m-1) double factorial                         c
c----------------------------------------------------------------------c
      ddfct(0)=1.d+00
      ddfct(1)=1.d+00
      ddfct(2)=3.d+00
      if (maxfac.gt.2) then
         do 20 i=3,maxfac
            ddfct(i)=(i+i-1)*ddfct(i-1)
   20    continue  
      endif
      return
      end
