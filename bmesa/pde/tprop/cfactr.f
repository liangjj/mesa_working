*deck cfactr.f
c***begin prologue     cfactr
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            fill gaussian wavepacket 
c***                   
c***references         
c
c***routines called    
c***end prologue       cfactr
      subroutine cfactr(p,eigc,wt,v,t,n,phase)
      implicit integer (a-z)
      real*8 p, eigc, wt
      complex*16 t, v
      complex*16 eye
      logical phase
      dimension p(n,*), t(*), eigc(n), wt(n), v(n)
      data eye/(0.d0,1.d0)/
      common/io/inp, iout
      if(.not.phase) then
         do 10 i=1,n
            t(i) = exp( - eigc(i)*eigc(i) )
            t(i) = t(i) * p(i,i) * wt(i)                
 10      continue
      else
         do 20 i=1,n
            t(i) = exp( - ( eye*cos(eigc(i)) + eigc(i)*eigc(i) ) )
            t(i) = t(i) * p(i,i) * wt(i)                
 20      continue
      endif
      return
      end       



