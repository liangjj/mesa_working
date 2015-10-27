*deck asymp
      subroutine asymp(fl,gl,dfl,dgl,ak,bk,r,rinv,eta,l,npt,
     1                 n,wron,level)
c***begin prologue     asymp
c***date written       920324   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose            asymptotic expansion for regular or irregular
c***                   positive coulomb function at large rho.
c***
c***description        f  = g * cos ( theta ) + f * sin ( theta )
c                       l                  l                   l
c***                   the leading term of f is 1.0
c***                   for the irregular function switch f and g and replace
c***                   the plus by a minus sign.
c***                   g, f and their derivatives are defined in the NBS
c***                   mathematical handbook.
c***references       
c
c***routines called
c
c***end prologue       asymp
c
      implicit integer (a-z)
      real*8 fl, gl, dfl, dgl, pi, carg, eta, sigmal, wron
      real*8  sn, cn, tol, thetal, r, rinv
      real*8 f, g, fd, gd, ak, bk, fold, gold, fnew, gnew
      real*8 foldd, goldd, fdnew, gdnew
      real*8  zero, half, one, two
      complex *16 cgamma, eye
      parameter (tol=1.d-14)
      dimension fl(npt), gl(npt), dfl(npt), dgl(npt)
      dimension ak(0:n), bk(0:n), r(npt), rinv(npt)
      common/io/inp,iout
      data pi /3.14159265358979323846d+00/
      data zero, half, one, two/0.d0,.5d0,1.d0,2.d0/
      data eye/(0.d0,1.d0)/
c**********************************************************************c
c        calculate the coulomb phase shift and the sin and cos         c
c        of the argument of the asymptotic thetal                      c
c**********************************************************************c
      sigmal=carg(cgamma(l+one+eye*eta))
      if(level.ge.1) then
         write(iout,1) sigmal
      endif
      do 10 pt=1,npt
         thetal=r(pt)-eta*log(two*r(pt))-l*pi*half+sigmal
         sn=sin(thetal)
         cn=cos(thetal)
c**********************************************************************c
c                do the series expansion to the level tol              c
c**********************************************************************c
         fold=one
         gold=zero
         foldd=zero
         goldd=one-eta*rinv(pt)
         f=fold
         g=gold
         fd=foldd
         gd=goldd
         do 20 i=0,n
            fnew=(ak(i)*fold-bk(i)*gold)*rinv(pt)
            gnew=(ak(i)*gold+bk(i)*fold)*rinv(pt)
            fdnew=(ak(i)*foldd-bk(i)*goldd-fnew)*rinv(pt)
            gdnew=(ak(i)*goldd+bk(i)*foldd-gnew)*rinv(pt)
            if (abs(fnew).lt.tol.and.abs(gnew).lt.tol) go to 30
            f=f+fnew
            g=g+gnew
            fd=fd+fdnew
            gd=gd+gdnew
            fold=fnew
            gold=gnew
            foldd=fdnew
            goldd=gdnew
   20    continue     
   30    if (level.ge.3) then
             write(iout,*) 'asymptotic expansion converged in ',i,
     1                     ' terms'
         endif
         fl(pt)=f*sn+g*cn
         dfl(pt)=fd*sn+gd*cn
         gl(pt)=f*cn-g*sn
         dgl(pt)=fd*cn-gd*sn
         wron=fl(pt)*dgl(pt)-dfl(pt)*gl(pt)
   10 continue     
      return
    1 format(/,5x,'coulomb phase shift',1x,e15.8,/)   
      end

