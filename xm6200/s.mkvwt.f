h27399
s 00048/00000/00000
d D 1.1 94/02/16 20:35:04 mesa 1 0
c date and time created 94/02/16 20:35:04 by mesa
e
u
U
f e 0
t
T
I 1
*deck %W%  %G%
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
      count=0
      if (typ.eq.'newton-cotes') then
          do 10 i=1,nr-1
             pcnt=0
             do 20 j=1,nr
                do 30 k=1,nprd
                   pcnt=pcnt+1
                   count=count+1
                   wt(count)=wt(count)*p(pcnt)
   30           continue
   20        continue
   10     continue
      else
          do 40 i=1,nr
             do 50 j=1,nprd
                count=count+1
                wt(count)=wt(count)*p(count)
   50        continue
   40     continue
          pcnt=count 
      endif                            
      return
      end

E 1
