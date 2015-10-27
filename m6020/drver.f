*deck @(#)drver.f	1.1 9/8/91
      subroutine drver(driver,expfn,psi,n)
      implicit integer (a-z)
      real *8 driver, expfn, psi
      dimension driver(n), expfn(n), psi(n)
      do 10 i=1,n
         driver(i)=expfn(i)*psi(i)
   10 continue
      return
      end        
