*deck wtfun.f
c***begin prologue     wtfun
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***purpose            weighting function 
c***                   
c***                   
c***                                                          
c***references         
c
c***routines called    
c***end prologue       wtfun
      subroutine wtfun(x,wtsi,wtso,type,n)
      implicit integer (a-z)
      real*8 x, wtsi, wtso
      character*(*) type
      dimension  x(n), wtsi(n), wtso(n)
      common/io/inp, iout 
      if(type.eq.'laguerre') then    
          do 10 i = 1, n
             wtso(i) = wtsi(i)*exp(-x(i))
 10       continue 
      elseif(type.eq.'hermite') then
          do 20 i = 1, n
             wtso(i) = wtsi(i)*exp(-x(i)*x(i))
 20       continue 
      elseif(type.eq.'one') then
           call copy(wtsi,wtso,n)
      else           
          call lnkerr('error in weight function')   
      endif
      return
      end       
