*deck %W%  %G%
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
      subroutine mkyunt (r,wt,yukawa,aprvol,yukint,nr,nthet,
     1                   nphi,typ,nang,nonsep)
      implicit integer (a-z)
      real*8 r, wt, yukawa, yukint, aprvol, rr
      dimension  wt(*), yukawa(*), r(nr)
      character*(*) typ
      logical nonsep
      common /io/ inp, iout
      if (nonsep) then
          nprd=nang
      else
          nprd=nthet*nphi
      endif              
      if (typ.eq.'newton-cotes') then
          count=0
          do 10 i=1,nr-1
             ptcnt=0
             do 20 j=1,nr
                rr=r(j)*r(j)
                do 30 k=1,nprd
                   ptcnt=ptcnt+1
                   count=count+1
                   aprvol=aprvol+wt(count)*rr 
                   yukint=yukint+yukawa(ptcnt)*wt(count)*rr
   30           continue
   20        continue
   10     continue
      else
          count=0
          do 50 i=1,nr
             rr=r(i)*r(i)
             do 60 j=1,nprd
                count=count+1
                yukint=yukint+wt(count)*yukawa(count)*rr
                aprvol=aprvol+wt(count)*rr 
   60        continue
   50     continue
      endif                            
      return
      end

