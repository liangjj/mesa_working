*deck gssums.f
c***begin prologue     gssums
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            R-matrix sums for symmetric wavefunction.
c***                   
c***references         
c
c***routines called    
c***end prologue       gssums
      subroutine gssums(vec,srf,sum,n,ni,prn)
      implicit integer (a-z)
      real*8 vec, srf, sum, fac
      character*80 title
      logical prn
      dimension vec(n,n), sum(ni,n)
      common/io/inp, iout
c
c
      write(iout,1) srf
 1    format(/,1x,'srf = ',e15.8)
      title='vectors'
      call prntrm(title,vec,n,n,n,n,iout)
      fac=1.d0/sqrt(2.d0)
      call rzero(sum,ni*n)
      do 10 i=1,ni-1
         index=ni*(ni-1)/2 +i
         do 20 q=1,n 
            sum(i,q) = vec(index,q)*srf*fac            
 20      continue   
         index=ni*(ni+1)/2
         do 30 q=1,n
            sum(i,q) = vec(index,q)*srf
 30      continue   
 10   continue   
      return      
      end       


