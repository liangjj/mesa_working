*deck filvt.f
c***begin prologue     filvt
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            vector fill for tdvr.
c***                   
c***references         
c
c***routines called    
c***end prologue       filvt
      subroutine filvt(vout,vin,n,nt,m)
      implicit integer (a-z)
      real*8 vout, vin
      dimension vout(n,nt,m), vin(n,m)
      common/io/inp, iout
      do 10 t=1,nt
         do 20 j=1,m
            call vadd(vout(1,t,j),vout(1,t,j),vin(1,j),n) 
 20      continue   
 10   continue
      return
      end       


