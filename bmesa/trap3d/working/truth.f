*deck truth.f
c***begin prologue     truth
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***                   
c***references         
c
c***routines called    
c***end prologue       truth
      subroutine truth(array,tf,nreg)
      implicit integer (a-z)
      logical array, tf
      dimension array(2,nreg)
      common/io/inp, iout
      do 10 i=1,nreg
         array(1,i)=tf
         array(2,i)=tf
 10   continue
      return
      end       
