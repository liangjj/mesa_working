*deck fcn
      subroutine fcn ( x, bsav, phi )
      implicit real *8 (a-h,o-z)
      parameter (npmax=20, nptmx=100)
      dimension x(1), bsav(1)
      common / pdata / np, npbck, elist(256), plist(256), sreal(256),
     >                 simag(256), phase(256), pabs(256), pres(256),
     >                 pbck(256), edup(256), ilist(256)
      dimension coef(npmax,npmax), xn(npmax,nptmx), ipvt(npmax)
c
c
      data pi2 / 1.5707963268d0 /, ncall /0/
c
      ncall = ncall + 1
      pmax = x(1)
      eres = x(2)
      gres = x(3)
      nlpar = npbck + 1
c
      do 1 i = 1, np
      pres(i) = pmax*(pi2 + atan(2.d0*(elist(i)-eres)/gres))
      pbck(i) = pabs(i) - pres(i)
 1    continue
c
      phi = 0d0
      rms = 0d0
c     call lstsq ( pbck, elist, x(4), np, nlpar, rms, v2fit1 )
      call lstsq ( pbck, elist, x(4), xn, coef, ipvt, np, nlpar, rms)
c
      do 2 i = 1, np
      pfit = pres(i)
      ep = 1d0
c
      do 3 j = 1, nlpar
      pfit = pfit + x(j+3)*ep
      ep = ep*elist(i)
 3    continue
c
      bsav(i) = pabs(i) - pfit
      phi = phi + bsav(i)**2
 2    continue
c
      return
      end

















