*deck @(#)genbes.f	1.1 9/8/91
c***begin prologue     m7000
c***date written       920525   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m7000, link 7000, spline
c***author             schneider, b. (nsf)
c***                   derived from earlier routines by t. n. rescigno (llnl)
c***source             m7000
c***purpose            regular and regularized bessel functions on grid
c***
c***description        a universal (as function of rho=k*r) regular
c***                   and regularized irregular function are computed,
c***                   and spline fit. the spline coefficients are used
c***                   to calculate the free functions for each l and energy
c***                   at the integration mesh using the routine mkbes.
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       m7000
      subroutine genbes(x,xinv,expfn,v,psi0,psi1,driver,j,jp,y,yp,norm,
     1                 fun,funout,dfnout,coefr,coefc,scr,rmin,rmax,
     2                 rdel,smldel,alpha,nint,subint,nr,lmax,ltop)
      implicit integer (a-z)
      real *8 x, xinv, expfn, v, psi0, psi1, driver, j, jp, y, yp
      real *8 norm, rmin, rmax, rdel, smldel, dl
      real *8 alpha, jnint, ynint, jnintm, ynintm, scale
      real *8 dfnout, coefr, energy, r0, r1
      complex *16 fun, funout, coefc, scr
      dimension x(nint), xinv(nint), expfn(nint), v(nint)
      dimension j(2,0:ltop), jp(2,0:ltop), y(2,0:ltop), yp(2,0:ltop)
      dimension norm(2), psi0(nint), psi1(nint), fun(nint)
      dimension funout(nr,0:lmax), dfnout(nr,0:lmax), scr(nr)
      dimension coefr(nr,0:lmax), coefc(nr,0:lmax), driver(nint)
      common /io/ inp, iout
      data energy / .5d+00 /
c---------------------------------------------------------------------c
c   note that rmin is the first spline knot point. if the real        c
c    grid goes below rmin, the cubic spline is being extrapolated     c
c                         (should be o.k.)                            c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                      generate spline data                           c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c             lmax = max l for splined bessel fcns                    c
c             rmax(rmin) = max(min) argument for                      c
c                          splined bessel fcns                        c
c             nper is the number of points to use per unit interval   c
c                         in computing the spline coefficients        c
c---------------------------------------------------------------------c  
c----------------------------------------------------------------------c
c           make integration grid and exponential function             c
c----------------------------------------------------------------------c
      call mkxexp(x,xinv,expfn,alpha,rmin,smldel,nint)
c
      do 200 l=0,lmax
c----------------------------------------------------------------------c
c              calculate effective potential for this l                c
c----------------------------------------------------------------------c
         call potntl(v,rmin,l,smldel,nint)
c----------------------------------------------------------------------c
c        integrate the homogeneous equation for the bessel             c
c                        function                                      c
c----------------------------------------------------------------------c
         call rzero(driver,nint)
         r0=rmin**(l+1)
         r1=(rmin+smldel)**(l+1)
         call numerv(psi0,r0,r1,driver,energy,v,smldel,nint)
c----------------------------------------------------------------------c
c              calculate regular and irregular bessel                  c
c              functions to normalize wavefunction                     c
c----------------------------------------------------------------------c
         dl=l
         if ( dl.gt.rmax) then
c----------------------------------------------------------------------c
c             backward recursion                                       c
c----------------------------------------------------------------------c
              call rcbesb(x(nint-1),xinv(nint-1),j,jp,y,yp,norm,2,2,l,
     1                    ltop,.false.)
         else
c----------------------------------------------------------------------c
c              forward recursion                                       c
c----------------------------------------------------------------------c
              call rcbesf(x(nint-1),xinv(nint-1),j,jp,y,yp,2,2,l,ltop,
     1                    .false.)
         endif
         jnint=j(2,l)
         ynint=y(2,l) 
         jnintm=j(1,l)
         ynintm=y(1,l)    
c----------------------------------------------------------------------c
c             calculate new driver                                     c
c----------------------------------------------------------------------c
         call drver(driver,expfn,psi0,nint)
c----------------------------------------------------------------------c
c               integrate inhomogeneous differential equation          c
c----------------------------------------------------------------------c
         call numerv(psi1,r0,r1,driver,energy,v,smldel,nint)
c----------------------------------------------------------------------c
c              compute the constants necessary to make                 c
c              the solution go to a pure neumann function.             c
c              then construct the final complex function.              c
c----------------------------------------------------------------------c
         call fnlfun(fun,scale,psi0,psi1,x,jnint,ynint,jnintm,ynintm,
     1               nint)
c----------------------------------------------------------------------c
c         put function and second derivative to be interpolated        c
c                        on to small grid                              c
c----------------------------------------------------------------------c
         call tolggr(fun,x,driver,funout(1,l),dfnout(1,l),l,scale,
     1               nint,subint,nr)

  200 continue
c----------------------------------------------------------------------c
c               calculate x on small grid                              c
c----------------------------------------------------------------------c
      call mkx(x,rmin,rdel,nr)
c----------------------------------------------------------------------c
c              make the spline coefficients                            c
c----------------------------------------------------------------------c
      do 300 l=0,lmax      
         call splinc(x,funout(1,l),rdel,nr,scr,coefc(1,l))
         call splinr(x,dfnout(1,l),rdel,nr,scr,coefr(1,l))
  300 continue
      return
      end


