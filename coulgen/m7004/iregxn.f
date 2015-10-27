*deck iregxn
      subroutine iregxn(gl,dgl,a0,b0,c0,prefac,r,rinv,l,npt,n,prnt)
c***begin prologue     iregxn
c***date written       920324   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose            series expansion for negative energy
c***                   irregular coulomb function at small rho.
c***description
c***references
c
c***routines called
c
c***end prologue       iregxn
c
      implicit integer (a-z)
      real*8 gl, dgl, r, rinv, a0, b0, c0, tworho, mtworo
      real*8 pre, add, zero, half, one, two, tol, prefac
      real*8 expfac, dxpfac, lgfac, dlgfac, rhofac, trofc, dtrofc
      real*8 csum, dcsum, asum, dasum, bsum, dbsum, term1, dterm1
      logical prnt
      dimension gl(npt), dgl(npt), a0(0:n), b0(0:n), c0(0:n), r(npt)
      dimension rinv(npt)
      common/io/inp,iout
      parameter (tol=1.d-30)
      data zero, half, one, two / 0.d0, .5d0, 1.d0, 2.d0 / 
c**********************************************************************c
c      gl = -1./sqrt(2.) * series * exp(-rho) * (2.*rho)**(-l) /       c
c                    (gamma(eta+l+1)*gamma(eta-l))                     c
c                                                                      c
c      series = sum(n=0 to n=2*l) * c0(n) * (-2.*rho)**n               c
c               - (-2.*rho)**(2*l+1) * { ln(2.*rho) *                  c        
c                                        sum(n=0 to n=infinity)        c
c                                            a0(n) * (2.*rho)**n -     c
c                                        sum(n=0 to n=infinity)        c
c                                            b0(n) * (2.*rho)**n }     c
c                                                                      c
c      * Henry has two misprints in his paper. the c0-sum has a        c
c        + 2.*rho and the last sum is added instead of subtracted.     c
c        these errors were discovered by comparing with the NBS        c
c        handbook and using the Henry normalization.                   c
c**********************************************************************c
    
      twoel1=l+l+1
      call rzero(gl,npt)
      call rzero(dgl,npt)
      do 10 pt=1,npt
c**********************************************************************c
c        exponential prefactor and its derivative                      c
c**********************************************************************c
         expfac=prefac*exp(-r(pt))
         dxpfac=-expfac
c**********************************************************************c
c             (2.*rho)**(-l) and derivative                            c
c             then combine with exponential factor                     c 
c**********************************************************************c
         rhofac=one
         if (l.ne.0) then
             rhofac=(half*rinv(pt))**l
             expfac=rhofac*expfac
             dxpfac=rhofac*dxpfac-l*expfac*rinv(pt)
         endif
c**********************************************************************c
c                prefactors in front of sums                           c
c**********************************************************************c
         tworho=two*r(pt)
         trofc = (-tworho)**twoel1
         dtrofc = twoel1*trofc*rinv(pt)
         lgfac=log(tworho)
         dlgfac=rinv(pt)
c**********************************************************************c
c                     do the six sums                                  c
c**********************************************************************c
         csum=zero
         dcsum=zero
         pre=one
         mtworo=-tworho
         do 20 i=0,l+l
            csum=csum+c0(i)*pre
            dcsum=dcsum+i*pre
            pre=pre*mtworo
   20    continue
         dcsum=dcsum*rinv(pt)
         pre=one
         asum=zero
         dasum=zero
         do 30 i=0,n
            add=a0(i)*pre
c**********************************************************************c
c            this funny business is done because accidental zeros      c
c            can occur in the a coefficients and at least a            c
c            reasonable result is likely to come out with 25 terms     c
c**********************************************************************c
            if (abs(add).lt.tol.and.i.ge.25) go to 40
            asum=asum+add
            dasum=dasum+i*add 
            pre=pre*tworho
   30    continue
   40    imax=i
         dasum=dasum*rinv(pt)
         pre=one
         bsum=zero
         dbsum=zero
         do 50 i=0,n
            add=b0(i)*pre
            if (abs(add).lt.tol.and.i.ge.25) go to 60
            bsum=bsum+add
            dbsum=dbsum+i*add 
            pre=pre*tworho
   50    continue
   60    imax=max(imax,i)
         if (prnt) then
             write (iout,*) 'power series for irregular function '//
     1                      'converged in ',imax,' terms'
         endif
         dbsum=dbsum*rinv(pt)
c**********************************************************************c
c        combine all the factors to get the irregular function         c
c                                and                                   c
c                             its   derivative                         c
c**********************************************************************c
         term1= csum - trofc * ( asum*lgfac - bsum )
         dterm1=dcsum - trofc  * ( dasum*lgfac + asum*dlgfac - dbsum )
     1                - dtrofc * ( asum*lgfac - bsum )
         gl(pt) = -expfac * term1
         dgl(pt) = -expfac * dterm1 -dxpfac*term1
   10 continue  
      return
      end


