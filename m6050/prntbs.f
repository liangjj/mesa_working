*deck @(#)prntbs.f	
c***begin prologue     prntbs
c***date written       920502   (yymmdd)
c***revision date      
c***keywords           prntbs link 6050
c***authors            Schneider,B (NSF)
c***                   
c***source             m6050
c***purpose            print basis information for expansion of
c***                   complex vlamdas in optical potential.
c***references       
c
c***routines called    
c***end prologue       vlm
      subroutine prntbs (lval,mval,gam,del,n)
      implicit real *8 (a-h,o-z)
      complex *16  gam, del
      dimension lval(n), mval(n), gam(n), del(n)
      common /io/ inp, iout
c
      write (iout,1)
      write(iout,2)
      do 10 i=1,n
         write(iout,3) lval(i),mval(i),gam(i),del(i)
   10 continue
    1 format (/,5x,'basis set information for optical potential')
    2 format(//,2x,'l value',2x,'m value',15x,'gamma',27x,'delta',//)
    3 format(2x,i4,5x,i4,4x,f12.6,5x,f12.6,3x,f12.6,5x,f12.6)
      return
      end
