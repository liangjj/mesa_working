*deck dvderr.f
c***begin prologue     dvderr
c***date written       991013   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           error print and exit for davidson code.
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***                   
c***description        
c***references         
c
c***routines called    
c***end prologue       dvderr
      subroutine dvderr(number)
      implicit integer (a-z)
      common/io/inp, iout
      if(number.eq.1) then
         write(iout,1)
      elseif(number.eq.2) then
         write(iout,2)
      elseif(number.eq.3) then
         write(iout,3)
      endif
      call lnkerr('quit davidson')
      return
 1    format(/,5x,'cannot even begin davidson calculation:',/,5x,
     1            'orthonormalization of initial vectors yields null'
     2            '  set')
 2    format(/1x,'no more orthonormal vectors can be added')
 3    format(/,5x,'cannot continue since maxvec is exceeded')
      end       

