*deck gfac.f
c***begin prologue     gfac
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           copy
c***author             schneider, barry (nsf)
c***source             
c***purpose            gaussian time factors
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       gfac
      subroutine gfac(pre,q,width,shift,n)
      implicit integer (a-z)
      real*8 pre, q, width, shift
      dimension pre(n), q(n)
      common/io/inp, iout
      do 10 i=1,n
         pre(i) = pre(i) * exp( - width *
     1                       ( q(i) - shift ) * 
     2                       ( q(i) - shift ) )
 10   continue  
      return
      end       


