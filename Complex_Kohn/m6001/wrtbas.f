*deck @(#)wrtbas.f	1.1 9/7/91
c***begin prologue     wrtbas
c***date written       890404   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           wrtbas, link 6001, basis
c***author             schneider, barry (lanl)
c***source             m6001
c***purpose            output basis information
c*** 
c
c***references       
c
c***routines called  
c***end prologue       wrtbas
      subroutine wrtbas (ntot)
      implicit integer (a-z)
      parameter (dimpr=300)
      real *8 alf, cont
      character *3 aosym
      character *128 itdum
      common /io/ inp, iout
      common /aosi/ npr, ncon, nxyzc(dimpr,4), nprc(dimpr)
      common /aosr/ alf(dimpr), cont(dimpr)
      common/chrpss/ itdum(2), aosym(dimpr)
      write (iout,30)
      ntot=0
      do 10 i=1,ncon
         do 10 j=1,nprc(i)
            ntot=ntot+1
      write (iout,20) i,ntot,aosym(i),nxyzc(ntot,4),alf(ntot),cont(ntot)
   10 continue
      return
c
   30 format(/,5x,'con fn',3x,'ao',2x,'sym',4x,'cen',9x,'exp',11x,'coef'
     1      )
   20 format(/,7x,i3,3x,i3,3x,a3,4x,i2,3x,f12.5,2x,f12.5)
      end
