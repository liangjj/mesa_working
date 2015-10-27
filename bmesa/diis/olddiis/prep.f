*deck prep.f
c***begin prologue     prep
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       prep
      subroutine prep(ham,temp,energy,n,prnt)
      implicit integer (a-z)
      real*8 ham, temp, energy
      character*80 title
      logical prnt
      dimension ham(n,n), temp(n,n)
      common/io/inp, iout
      call copy(ham,temp,n*n)
      do 10 i=1,n
         do 20 j=1,n
            temp(i,j) = - temp(i,j)
 20      continue
         temp(i,i) = energy + temp(i,i)   
 10   continue   
      if(prnt) then
         title='matrix'
         call prntrm(title,temp,n,n,n,n,iout)
      endif
      return
      end       
