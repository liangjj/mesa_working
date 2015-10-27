*deck genrwt.f
c***begin prologue     genrwt
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           reference weight
c***author             schneider, barry (nsf)
c***source
c***purpose            calculate ratio of  actual to reference 
c***                   quadrature weight.
c***                   
c***description             
c***                                                                       
c***                                                          
c***references         
c
c***routines called    
c***end prologue       genrwt
      subroutine genrwt(rwt,pt,wt,alpha,beta,n)
      implicit integer (a-z)
      real*8 rwt, pt, alpha, beta
      character*(*) wt
      dimension rwt(n), pt(n)
      common/io/inp, iout 
      if(wt.eq.'legendre') then
         call vfill(rwt,1.d0,n)
      elseif(wt.eq.'hermite') then
         do 10 i=1,n
            rwt(i)=exp(-pt(i)*pt(i))
 10      continue   
      elseif(wt.eq.'chebyshev-1') then
         do 20 i=1,n
            rwt(i)= sqrt ( 1.d0/( 1.d0 -pt(i)*pt(i)) )
 20      continue   
      elseif(wt.eq.'chebyshev-2') then
         do 30 i=1,n
            rwt(i)= sqrt (  1.d0 -pt(i)*pt(i) ) 
 30      continue   
      elseif(wt.eq.'laguerre') then
         do 40 i=1,n
            rwt(i)= pt(i)**alpha*exp(-pt(i)) 
 40      continue   
      elseif(wt.eq.'jacobi') then
         do 50 i=1,n
            rwt(i)= (1.d0-pt(i))**alpha * (1.d0+pt(i))**beta 
 50      continue   
      elseif(wt.eq.'rys') then
         do 60 i=1,n
            rwt(i)=exp(-alpha*pt(i)*pt(i))
 60      continue   
      endif
      return
      end       
