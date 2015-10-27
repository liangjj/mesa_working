*deck toplm.f
c***begin prologue     toplm
c***date written       921223   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           projection, ylm, plm
c***author             schneider, barry (nsf)
c***source
c***purpose            compute projection of function onto spherical
c***                   harmonics for a given m value where the angular
c***                   quadrature is separable in theta and phi.             
c***references         none
c
c***routines called
c***end prologue       toplm
      subroutine toplm (plm,phifn,flm,f,wtthe,wtphi,scr,l,m,no,nthet,
     1                  nphi,nr,title,prnt)
      implicit integer (a-z)
      real*8 plm, phifn, flm, f, wtthe, wtphi, scr
      character*80 title
      logical prnt
      dimension plm(nthet,m:l), phifn(nphi,no), flm(nr,m:l,no)
      dimension f(nr,nthet,nphi), wtthe(nthet,2), wtphi(nphi,2)
      dimension scr(nr,nthet,no)
      common /io/ inp, iout
      dim=l-m+1
      call vmmul(wtthe,plm(1,m),plm(1,m),nthet,dim)
      call vmmul(wtphi,phifn,phifn,nphi,no)
      call ebc(scr,f,phifn,nr*nthet,nphi,no)
      do 10 mtot=1,no
         call ebc(flm(1,m,mtot),scr(1,1,mtot),plm(1,m),nr,nthet,dim)
         call setzro(flm(1,m,mtot),nr*dim)
         if (prnt) then
             write(iout,1) mtot, m
             call prntfm(title,flm(1,m,mtot),nr,dim,nr,dim,iout)
         endif
   10 continue         
      call vmmul(wtthe(1,2),plm(1,m),plm(1,m),nthet,dim)
      call vmmul(wtphi(1,2),phifn,phifn,nphi,no)    
      return
    1 format(/,'the ',i2,' ylm decomposition matrix for m = '
     1        ,i3)
c
      end


