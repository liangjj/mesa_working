*deck radial.f
c***begin prologue     radial
c***date written       001230   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            radial functions
c***references         
c
c***routines called    
c***end prologue       radial
      subroutine radial(rho,jmk,djmk,nmk,dnmk,jb,djb,yb,dyb,
     1                  fact,top,m,type)
      implicit integer (a-z)
      real*8 rho, jmk, djmk, nmk, dnmk, jb, djb, yb, dyb
      real*8 ainv, sqainv, a2inv, a32inv, a52inv, a72inv, r1mach
      real*8 zero, half, one, two, valm0
      real*8 atan, rhomin, fac, pi, pi2, fact
      character*16 dir
      character*(*) type
      dimension jb(0:*), djb(0:*), yb(0:*), dyb(0:*)
      dimension fact(0:*)
      common/io/inp, iout
      data pi  / 3.141592653589793238462643d+00 /
      data pi2 / 1.5707963267948966192313215d+00 /
      data zero, half, one, two / 0.d0, .5d0, 1.d0, 2.d0 /
      data valm0 / .125d0 / 
      data dir /'derivatives'/
      data rhomin / 1.d-14 /
      ainv=one/rho
      sqainv=sqrt(ainv)
      a2inv=ainv*ainv
      a32inv=ainv*sqainv
      a52inv=a2inv*sqainv
      a72inv=a52inv*ainv
      fac=sqrt(2.d0/pi)
      twom=m+m+2
      if(type.eq.'hyperspherical') then
         if(rho.gt.rhomin) then
            call bessel(rho,jb,djb,yb,dyb,fact,1,top,.false.)
            jmk = jb(twom)*a2inv
            nmk = yb(twom)*a2inv
            djmk= djb(twom)*a2inv - two*jmk*ainv
            dnmk= dyb(twom)*a2inv - two*nmk*ainv
         else
            jmk=r1mach(1)
            djmk=r1mach(1)
            nmk=r1mach(2)
            dnmk=r1mach(2)
            if(m.eq.0) then
               jmk=valm0
            endif
	 endif   
      elseif(type.eq.'bessel') then
         call bessel(rho,jb,djb,yb,dyb,fact,1,top,.false.)
         jmk = jb(m)
         nmk = yb(m)
         djmk= djb(m)
         dnmk= dyb(m)
      elseif(type.eq.'asymptotic-bessel') then
         jmk = cos(rho-.5d0*m*pi - .25d0*pi)
         nmk = sin(rho-.5d0*m*pi - .25d0*pi)
         djmk = - nmk*sqainv - .5d0*jmk*a32inv
         dnmk =   jmk*sqainv - .5d0*nmk*a32inv
         jmk = jmk*sqainv*fac
         nmk = nmk*sqainv*fac
         djmk = djmk*fac
         dnmk = dnmk*fac
      elseif(type.eq.'asymptotic') then
         jmk = cos(rho-(m+1)*pi - .25d0*pi)
         nmk = sin(rho-(m+1)*pi - .25d0*pi)
         djmk = - nmk*a52inv - 2.5d0*jmk*a72inv
         dnmk =   jmk*a52inv - 2.5d0*nmk*a72inv
         jmk = jmk*a52inv*fac
         nmk = nmk*a52inv*fac
         djmk = djmk*fac
         dnmk = dnmk*fac
      else
         call lnkerr('screw up')
      endif
      return
      end
