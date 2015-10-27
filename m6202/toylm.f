*deck toylm.f
c***begin prologue     toylm
c***date written       921223   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           projection, ylm
c***author             schneider, barry (nsf)
c***source
c***purpose            compute projection of function onto spherical
c***                   harmonics for a given m value where the angular
c***                   quadrature in non-separable in theta and phi.              
c***references         none
c
c***routines called
c***end prologue       toylm
      subroutine toylm (ylm,flm,f,wtang,l,m,no,nang,nr,title,prnt)
      implicit integer (a-z)
      real*8 ylm, flm, wtang, f
      character*80 title
      logical prnt
      dimension ylm(nang,m:l,no), flm(nr,m:l,no), wtang(nang,2)
      dimension f(nr,nang)
      common /io/ inp, iout
      dim=l-m+1
      do 10 mtot=1,no
         call vmmul(wtang,ylm(1,m,mtot),ylm(1,m,mtot),nang,dim)
         call ebc(flm(1,m,mtot),f,ylm(1,m,mtot),nr,nang,dim)
         call setzro(flm(1,m,mtot),nr*dim)
         if (prnt) then
             write(iout,1) mtot, m
             call prntfm(title,flm(1,m,mtot),nr,dim,nr,dim,iout)
         endif
         call vmmul(wtang(1,2),ylm(1,m,mtot),ylm(1,m,mtot),nang,dim)
   10 continue         
      return
    1 format(/,'the ',i2,' ylm decomposition matrix for m = '
     1        ,i3)
c
      end


