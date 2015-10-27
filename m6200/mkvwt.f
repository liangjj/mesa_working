*deck mkvwt.f
c***begin prologue     mkvwt
c***date written       930802   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           mkvwt, link m6200
c***author             schneider, barry (nsf)
c***source             m6200
c***purpose            calculate voronoi weights
c***references         
c
c***routines called
c***end prologue       mkvwt
      subroutine mkvwt (p,wt,nr,nthet,nphi,pcnt,count,typ,nang,nonsep)
      implicit integer (a-z)
      real*8 p, wt
      dimension p(*), wt(*)
      character*(*) typ
      logical nonsep
      common /io/ inp, iout
      if (nonsep) then
          nprd=nang
      else
          nprd=nthet*nphi
      endif
      pcnt=nr*nprd              
      if (typ.eq.'newton-cotes') then
          count=0
          do 10 i=1,nr-1
             call vmul(wt(count+1),wt(count+1),p,nr*nprd)
             count=count+nr*nprd
   10     continue
      else
          call vmul(wt,wt,p,nr*nprd) 
      endif                            
      return
      end

