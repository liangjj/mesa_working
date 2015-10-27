*deck prnfn.f
c***begin prologue     prnfn
c***date written       001230   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            print routine for hyperspherical code
c***references         
c
c***routines called    
c***end prologue       prnfn
      subroutine prnfn(f,var,title,label,n1,n2)
      implicit integer (a-z)
      real*8 f, var 
      character*(*) title, label
      dimension f(n2,n1), var(n2,n1)
      common/io/inp, iout
      write(iout,*) title
      write(iout,1) label
      do 10 i=1,n1
         do 20 j=1,n2
            write(iout,2) var(j,i), f(j,i)
 20      continue   
 10   continue   
      return
 1    format(/,'    ',a6,'     ',1x,'    function   ')
 2    format(1x,d15.8,1x,d15.8)
      end       















