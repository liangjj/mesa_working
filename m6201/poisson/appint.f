*deck appint
      subroutine appint(y,int0,psilm,temp,n)
c***begin prologue     appint
c***date written       930608   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose            asymptotic evaluation of forward integral
c***
c***references
c
c***routines called
c
c***end prologue       appint
c
      implicit integer (a-z)
      dimension y(n), psilm(n-1), temp(n-1)
      real*8 y, psilm, temp, int0
      common /io/ inp, iout
      call vfill(psilm,int0,n-1)
      call copy(psilm,temp,n-1)
      call vmul(psilm,psilm,y(2),n-1)
      return
      end















