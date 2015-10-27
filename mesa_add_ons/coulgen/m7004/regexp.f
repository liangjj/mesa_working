*deck regexp
      subroutine regexp(fl,dfl,r,a,l,cl,npt,n,prnt)
c***begin prologue     regexp
c***date written       920324   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose            series expansion for positive energy
c***                   regular coulomb function at small rho.
c***description        
c***references         NBS handbook
c
c***routines called
c
c***end prologue       regexp
c
      implicit integer (a-z)
      real*8 fl, dfl, r, a, cl, pre, add, one, tol
      logical prnt
      parameter ( tol=1.d-30)
      dimension a(l+1:l+1+n), fl(npt), dfl(npt), r(npt)
      common/io/inp,iout
      data one/1.d0/ 
      call rzero(fl,npt)
      call rzero(dfl,npt)
      do 10 pt=1,npt
         pre=one
         count=0 
         do 20 i=l+1,l+1+n
            count=count+1 
            add=a(i)*pre
c**********************************************************************c
c            this funny business is done because accidental zeros      c
c            can occur in the a coefficients and at least a            c
c            reasonable result is likely to come out with 25 terms     c
c**********************************************************************c
            if (abs(add).lt.tol.and.count.ge.25) go to 30
            fl(pt)=fl(pt)+add
            dfl(pt)=dfl(pt)+i*add
            pre=pre*r(pt)
   20    continue
   30    if (prnt) then
             write (iout,*) 'power series converged in ',count,' terms'
         endif
   10 continue  
      do 40 pt=1,npt
         pre=r(pt)**l
         fl(pt)=cl*fl(pt)*pre*r(pt)
         dfl(pt)=cl*dfl(pt)*pre 
   40 continue     
      return
      end



