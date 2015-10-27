*deck modcof.f
c***begin prologue     modcof
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***purpose            compute modified coefficients for orthogonal 
c***                   polynomials when end points are fixed.
c***                   
c***                                                          
c***references         
c
c***routines called    
c***end prologue       modcof
      subroutine modcof(a,b,endpts,nfixed,n,prncof)
      implicit integer (a-z)
      real*8 a, b, endpts, gbslve, t1, gam
      logical prncof
      character*80 title
      dimension a(n), b(n), endpts(2)
      common/io/inp, iout 
      if (nfixed.eq.1) then
c
c         only a(n) must be changed
c
          a(n) =gbslve(endpts(1), n, a, b)*b(n-1)**2 + endpts(1)
      endif
      if (nfixed.eq.2) then
c
c         a(n) and b(n-1) must be recomputed
c
          gam =gbslve(endpts(1), n, a, b)
          t1 = ((endpts(1) - endpts(2))/
     1           (gbslve(endpts(2), n, a, b) - gam))
          b(n-1) =  sqrt(t1)
          a(n) = endpts(1) + gam*t1
      endif
      if (prncof) then
          title='modified a coefficients for fixed endpoints'
          call prntrm(title,a(1),n,1,n,1,iout)
          title='modified b coefficients for fixed endpoints'
          call prntrm(title,b(1),n-1,1,n-1,1,iout)
      endif
      call iosys('write real "modified polynomial a coefficients" '//
     1           'to lamdat',n,a,0,' ')
      call iosys('write real "modified polynomial b coefficients" '//
     1           'to lamdat',n,b,0,' ')
      return
      end       
