*deck mkylm.f
c***begin prologue     mkylm
c***date written       921223   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           make spherical harmonics from p(l,m) and phifn(m)
c***author             schneider, barry (nsf)
c***source             
c***purpose            y(l,m) functions
c***references         none
c
c***routines called
c***end prologue       mkylm
      subroutine mkylm (plm,phifn,ylm,nang,lmax,m)
      implicit integer (a-z)
      real *8 plm, phifn, ylm
      dimension plm(nang,m:lmax), phifn(nang,2), ylm(nang,m:lmax,2)
      common /io/ inp, iout
      dim=lmax-m+1
      upper=2
      if (m.eq.0) then
          upper=1
      endif     
      do 10 mu=1,upper
          call vmmul(phifn(1,mu),plm(1,m),ylm(1,m,mu),nang,dim)
   10 continue               
      return
c
      end


