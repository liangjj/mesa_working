*deck coefs.f
c***begin prologue     coefs
c***date written       940213   (yymmdd)
c***revision date               (yymmdd)
c***keywords           expansion coefficients
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***description          
c***references       
c
c***routines called
c***end prologue       coefs
      subroutine coefs(c,f,wt1,f1,wt2,f2,n1,n2)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 c, f, wt1, f1, wt2, f2, fac1, fac2
      dimension c(n2,n1), f(n2,n1), wt1(n1), f1(n1,n1)
      dimension wt2(n2), f2(n2,n2)
      do 10 i=1,n1
         fac1=f1(i,i)*wt1(i)
         do 20 j=1,n2
            fac2=f2(j,j)*wt2(j)
            c(j,i)=fac1*fac2*f(j,i)
 20      continue
 10   continue
      return
      end
