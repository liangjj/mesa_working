*deck mkyunt.f
c***begin prologue     mkyunt
c***date written       930802   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           mkyunt, link m6200
c***author             schneider, barry (nsf)
c***source             m6200
c***purpose            calculate yukawa integral
c***references         
c
c***routines called
c***end prologue       mkyunt
      subroutine mkyunt (r,wt,yukawa,aprvol,yukint,scr,nr,nthet,
     1                   nphi,typ,nang,nonsep)
      implicit integer (a-z)
      real*8 r, wt, yukawa, yukint, aprvol, rr, scr, sumel, sdot
      dimension  wt(*), yukawa(*), r(nr), scr(*)
      character*(*) typ
      logical nonsep
      common /io/ inp, iout
      if (nonsep) then
          nprd=nang
      else
          nprd=nthet*nphi
      endif              
      if (typ.eq.'newton-cotes') then
          ntot=nr*nprd
          count=0
          do 10 i=1,nr-1
             do 20 j=1,nr
                rr=r(j)*r(j)
                call vscale(scr(count+1),wt(count+1),rr,nprd)
                count=count+nprd
   20        continue
   10     continue
          aprvol=aprvol+sumel(scr,count)
          count=0
          do 30 i=1,nr-1
             ptcnt=0
             yukint=yukint+sdot(ntot,yukawa,1,scr(count+1),1)
             count=count+ntot
   30     continue
      else
          count=0
          do 50 i=1,nr
             rr=r(i)*r(i)
             call vscale(scr(count+1),wt(count+1),rr,nprd)
             count=count+nprd
   50     continue
          yukint=yukint+sdot(count,scr,1,yukawa,1)             
          aprvol=aprvol+sumel(scr,count) 
      endif                            
      return
      end

