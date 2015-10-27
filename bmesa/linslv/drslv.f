*deck drslv.f
c***begin prologue     drslv
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***                   
c***references         
c
c***routines called    
c***end prologue       drslv
      subroutine drslv(ham,v,rhs,energy,ipvt,n,m)
      implicit integer (a-z)
      real*8 ham, v, rhs, energy
      character*80 title
      dimension ham(n,n), v(n,n), rhs(n,m), ipvt(n) 
      common/io/inp, iout
      do 10 i=1,n
         do 20 j=1,n
            ham(i,j) = ham(i,j) + v(i,j)
            ham(i,j) = -ham(i,j)
   20    continue 
         ham(i,i) = energy + ham(i,i)        
   10 continue
c      title='matrix'
c      call prntrm(title,ham,n,n,n,n,iout)      
c      title='rhs'
c      call prntrm(title,rhs,n,m,n,m,iout)      
      call sgefa(ham,n,n,ipvt,info)
      do 30 i=1,m
         call sgesl(ham,n,n,ipvt,rhs(1,i),0)
   30 continue
      title='solution matrix'
      call prntrm(title,rhs,n,m,n,m,iout)                
      return
      end       
