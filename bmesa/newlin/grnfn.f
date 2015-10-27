c $Header$
*deck grnfn.f
c***begin prologue     grnfn
c***date written       930201   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           grnfn, link 6201, wavefunction
c***author             schneider, barry (nsf)
c***source             m6203
c***purpose            oscillatory or exponential solutions to the free
c***                   wave schroedinger equation.
c***description        regular and irregular spherical or modified
c***                   spherical bessel functions are computed using 
c***                   recursion relations.  
c***                    
c***                   wronskian of basic solutions is defined as:
c***                   wron = ( regular*deriv. of irregular
c***                                  - deriv. of regular*regular)
c***
c***references         
c***routines called
c***end prologue       grnfn
      subroutine grnfn (greg,gireg,energy,pt,scr,nr,lmax,ldim,prnt)
      implicit integer(a-z)
      parameter (acc=30 )
      real*8 greg, gireg, energy, pt, scr, dum, top, k, wron
      real*8 kappa
      character*80 title
      character*3 itoc 
      logical prnt
      common /io/ inp, iout
      dimension greg(nr,0:ldim), gireg(nr,0:ldim), pt(nr), scr(nr,4)
*
*
*          the radial functions are defined as solutions
*          of the equation;
*            g'' + 2*abs(e) - l*(l+1)/(r+r) = delta(r-r')
*                -
*
*
c     calculate the top l value needed to perform backward recursion
c     to accuracy acc. see numerical recipes for discussion. 
      top=lmax+sqrt(dfloat(lmax*acc))
      ltop=top+1
      if (ldim.lt.ltop) then
          call lnkerr('quit. ldim too small in grnfn')
      endif 
      if (energy.ge.0.d0) then
c                   positive energy code
c         make the variable k*r
          k=sqrt(energy)
          call copy(pt,scr(1,1),nr)
          call sscal(nr,k,scr(1,1),1)
c         calculate bessel functions at all k*r points
          call rcbes(scr(1,1),greg,dum,gireg,dum,wron,scr(1,2),nr,lmax,
     1               ltop,'no derivatives',.false.)
c         the wronskian is one when derivatives are defined wrt k*r
c
c         convert wronskian to factor needed to make greens
c         function have unit jump in derivative. we must associate
c         this with one of the two solutions or carry it around as
c         an additional factor. I will use it to scale the iregular
c         function.
          wron=1.d0/k
          call sscal(nr*(lmax+1),wron,gireg,1)
c                   negative energy code
      elseif(energy.lt.0.d0) then
c         make the variable kappa*r
          kappa=sqrt(abs(energy))
          call copy(pt,scr(1,1),nr)
          call sscal(nr,kappa,scr(1,1),1)
c         calculate modified bessel functions at all kappa*r points
          call modbes(scr(1,1),greg,dum,gireg,dum,wron,scr(1,2),nr,lmax,
     1                ltop,'no derivatives',.false.)
c         the wronskian is minus one when derivatives are defined wrt kappa*r
c
c         convert wronskian to factor needed to make greens
c         function have unit jump in derivative. we must associate
c         this with one of the two solutions or carry it around as
c         an additional factor. I will use it to scale the irregular
c         function.
          wron=-1.d0/kappa
          call sscal(nr*(lmax+1),wron,gireg,1)
      endif
      if (prnt) then
          title='regular solution for l = 0 to l ='//itoc(lmax)
          call prntrm (title,greg,nr,lmax+1,nr,lmax+1,iout)
          title='irregular boundary solution for l = 0 to l ='//
     1           itoc(lmax)
          call prntrm (title,gireg,nr,lmax+1,nr,lmax+1,iout)
      endif
      return
      end


