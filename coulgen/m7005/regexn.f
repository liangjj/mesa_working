*deck regexn
      subroutine regexn(fl,dfl,r,rinv,a0,l,npt,n,level)
c***begin prologue     regexn
c***date written       920324   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose            series expansion for negative energy
c***                   regular coulomb function at small rho.
c***description
c***references         NBS handbook
c
c***routines called
c
c***end prologue       regexn
c
      implicit integer (a-z)
      real*8 fl, dfl, r, rinv, a0, pre, add, one, two, dsqrt2, tol
      real*8 tworho
      parameter ( tol=1.d-14)
      dimension a0(0:n), fl(npt), dfl(npt), r(npt), rinv(npt)
      common/io/inp,iout
      data one, two, dsqrt2 / 1.d0, 2.d0, .70710678118654752440d0/
      call rzero(fl,npt)
      call rzero(dfl,npt)
      do 10 pt=1,npt
         pre=one
         tworho=two*r(pt)
         do 20 i=0,n
            add=a0(i)*pre
c**********************************************************************c
c            this funny business is done because accidental zeros      c
c            can occur in the a coefficients and at least a            c
c            reasonable result is likely to come out with 25 terms     c
c**********************************************************************c
            if (abs(add).lt.tol.and.i.ge.25) go to 30
            fl(pt)=fl(pt)+add
            dfl(pt)=dfl(pt)+i*add
            pre=pre*tworho
   20    continue
   30    if (level.ge.3) then
             write (iout,*) 'power series for regular solution '//
     1                      'converged in ',i,' terms'
         endif   
   10 continue  
      do 40 pt=1,npt
         tworho=two*r(pt)
         pre=dsqrt2*(tworho)**(l+1)
         add=exp(-r(pt))
         fl(pt)=fl(pt)*pre*add
         dfl(pt)=fl(pt)*( -one +(l+1)*rinv(pt) ) + 
     1                     pre*add*dfl(pt)*rinv(pt)
   40 continue     
      return
      end



