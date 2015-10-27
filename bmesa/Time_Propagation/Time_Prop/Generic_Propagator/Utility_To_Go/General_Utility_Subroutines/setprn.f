*deck setprn.f
c***begin prologue     setprn
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            set information.
c***                   
c***references         
c
c***routines called    
c***end prologue       setprn
      subroutine setprn(prn,n)
      implicit integer (a-z)
      logical prn
      dimension prn(n) 
      common/io/inp, iout
      do 10 i=1,n
         prn(i)=.true.
 10   continue   
      return      
      end       






