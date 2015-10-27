*deck mkwts.f
c***begin prologue     mkwts
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***purpose            integration weights.
c***                   
c***                                                          
c***references         
c
c***routines called    
c***end prologue       mkwts
      subroutine mkwts(ham,wts,c,n,npts)
      implicit integer (a-z)
      real*8 ham, wts, c, sum
      character*80 title
      dimension ham(npts,0:n-1), wts(npts), c(n)
      common/io/inp, iout 
      sum=0.d0 
      do 10 i=1,n
         c(i)=0.d0
         do 20 j=1,npts
            c(i)=c(i)+ham(j,i-1)*wts(j)
 20      continue
         sum=sum+c(i)*c(i)   
 10   continue
      write(iout,1) sum
      title='integral of polynomial'
      call prntrm(title,c,n,1,n,1,iout)
      return
    1 format(/,'sum of square of integral = ',e15.8)      
      end       
