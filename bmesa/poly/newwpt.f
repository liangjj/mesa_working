*deck newwpt.f
c***begin prologue     newwpt
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***purpose            gaussian points and weights based on polynomial 
c***                   recursion coefficients.
c***                   
c***                                                          
c***references         
c
c***routines called    
c***end prologue       newwpt
      subroutine newwpt(a,b,x,wts,norm0,dum,n,prnwpt)
      implicit integer (a-z)
      real*8 a, b, x, wts, norm0, dum, sum
      logical prnwpt
      character*80 title
      dimension a(n), b(n), x(n), wts(n), dum(n)
      common/io/inp, iout 
      wts(1) = 1.0d0
      do 10 i = 2, n
         wts(i) = 0.0d0
 10   continue   
c
      call copy(a,x,n)
      call copy(b,dum,n)
      call gbtql2 (n, x, dum, wts, ierr)
      if (ierr.ne.0) then
          call lnkerr('quit')
      endif
      sum=0.d0
      do 20 i = 1, n
         wts(i) = norm0 * wts(i) * wts(i)
         sum=sum+wts(i)
 20   continue   
      write(iout,1) sum
      if (prnwpt) then
          title='new gausssian points'
          call prntrm(title,x,n,1,n,1,iout)
          title='new gaussian weights'
          call prntrm(title,wts,n,1,n,1,iout)
      endif
      return
 1    format(/,1x,'sum of the weights = ',e15.8)
      end       
