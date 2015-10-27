*deck iregxp
      subroutine iregxp(gl,dgl,fl,dfl,r,b,l,dl,prefac,
     1                  qlpl,npt,n,prnt)
c***begin prologue     iregxp
c***date written       920324   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose            series expansion for positive energy
c***                   irregular coulomb function at small rho.
c***description
c***references         NBS handbook
c
c***routines called
c
c***end prologue       iregxp
c
      implicit integer (a-z)
      real*8 gl, dgl, fl, dfl, r, b, dl, prefac, pre, qlpl
      real*8 add, one, two, argm, tol, wron
      logical prnt
      dimension gl(npt), dgl(npt), fl(npt), dfl(npt)
      dimension b(-l:-l+n), r(npt)
      common/io/inp,iout
      parameter (tol=1.d-30)
      data one, two/ 1.d0, 2.d0 /
c**********************************************************************c
c        gl = 2. * eta * fl * { ln(2.*rho) + ql/pl } / c0**2           c
c                      -dl * rho**(-l) * sum (n=0 to n=infinity)       c
c                                        b(n) * rho**n                 c
c**********************************************************************c
      call rzero(gl,npt)
      call rzero(dgl,npt)
c**********************************************************************c
c                  do the theta series sum                             c
c**********************************************************************c
      do 10 pt=1,npt
         pre=one
         count=0
         do 20 i=-l,-l+n
            count=count+1
            add=b(i)*pre
c**********************************************************************c
c            this funny business is done because accidental zeros      c
c            can occur in the a coefficients and at least a            c
c            reasonable result is likely to come out with 25 terms     c
c**********************************************************************c
            if (abs(add).lt.tol.and.count.ge.25) go to 30
            gl(pt)=gl(pt)+add
            dgl(pt)=dgl(pt)+i*add
            pre=pre*r(pt)
   20    continue
   30    if (prnt) then
             write (iout,*) 'power series for irregular function'//
     1                      ' converged in ',count,' terms'
         endif
   10 continue  
      do 40 pt=1,npt
         pre=r(pt)**(-l-1)
         gl(pt)=dl*gl(pt)*pre*r(pt)
         dgl(pt)=dl*dgl(pt)*pre
   40 continue
c**********************************************************************c     
c               continue with rest of wavefunction                     c
c**********************************************************************c
      do 50 pt=1,npt
         argm=log(two*r(pt))+qlpl
         gl(pt)=gl(pt)+prefac*fl(pt)*argm
         dgl(pt)=dgl(pt)+prefac*(dfl(pt)*argm+fl(pt)/r(pt))
         wron=fl(pt)*dgl(pt)-dfl(pt)*gl(pt)
   50 continue     
      return
      end


