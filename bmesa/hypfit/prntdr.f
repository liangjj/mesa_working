*deck prntdr.f
c***begin prologue     prntdr
c***date written       001230   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            print derivatives in hyperspherical code
c***references         
c
c***routines called    
c***end prologue       prntdr
      subroutine prntdr(exact,approx,var,title,line,n)
      implicit integer (a-z)
      real*8 exact, approx, var 
      character*(*) title, line
      dimension exact(n), approx(n), var(n)
      common/io/inp, iout
      write(iout,*) title
      write(iout,*) line
      do 10 ji=1,n
         write(iout,1) var(ji), exact(ji), approx(ji)
 10   continue   
      return
 1    format(1x,d15.8,1x,d15.8,1x,d15.8)
      end       















