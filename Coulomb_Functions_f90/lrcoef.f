*deck lrcoef
      subroutine lrcoef(ak,bk,e0,eta,l,n,energy,level)
c***begin prologue     lrcoef
c***date written       920324   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose            expansion coefficients needed to construct
c***                   asymptotic values for regular or irregular
c***                   coulomb functions at large rho.
c***
c***references         nbs handbook
c
c***routines called
c
c***end prologue       lrcoef
c
      implicit integer (a-z)
      real*8 ak, bk, e0, eta, etasq, numfac
      character *(*) energy
      dimension ak(0:n), bk(0:n), e0(0:n)
      common/io/inp,iout
      if (energy.eq.'positive') then
c**********************************************************************c
c                coefficients for asymptotic expansion                 c
c                of regular and irregular positive energy              c
c                        coulomb functions                             c
c**********************************************************************c
          etasq=eta*eta
          numfac=l*(l+1)+etasq
          do 10 k=0,n
             k1=k+k+1
             k2=k1+1
             kfac=k*(k+1)
             ak(k)=k1*eta/k2
             bk(k)=(numfac-kfac)/k2
   10     continue      
          if (level.ge.2) then
              write(iout,*)
              write(iout,*) '          the a and b long range '//
     1                             'expansion coefficients'
              write(iout,*)
              write(iout,*) '      i          a                 b'
              do 50 i=0,n
                 write(iout,1) i,ak(i), bk(i)
   50         continue
          endif     
      elseif (energy.eq.'negative') then
c**********************************************************************c
c             coefficients for asymtotic expansion of irregular        c
c             negative energy coulomb functions.                       c
c**********************************************************************c
          e0(0)=1.d0
          do 60 i=1,n
             e0(i)=(i+eta+l)*(i+eta-l-1)*e0(i-1)/i
   60     continue
          if (level.ge.2) then
              write(iout,*)
              write(iout,*) '          the e long range expansion '//
     1                      'coefficients'
              write(iout,*)
              write(iout,*) '      i          e'
              do 70 i=0,n
                 write(iout,2) i,e0(i)
   70         continue
          endif
      endif     
      return
    1 format(5x,i3,3x,e15.8,3x,e15.8)   
    2 format(5x,i3,3x,e15.8)   
      end













