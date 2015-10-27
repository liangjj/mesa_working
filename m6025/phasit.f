*deck phasit
      subroutine phasit ( e, pin, n, const, nf, sdf, a, pout )
      implicit real *8 (a-h,o-z)
      dimension e(n), pin(n), pout(n), a(n)
      common / opts / dbgphs
      logical dbgphs
      data pi / 3.14159265358979d0 /
c
      a(1) = 0d0
      pout(1) = pin(1)
      if ( n .eq. 1 ) return
c      const = npi*180.d0
      cinv = 1d0/const
c
      if ( dbgphs ) print 1001
 1001 format ( " phasit -- flattening the data" )
      isml = 1
      esml = e(2) - e(1)
      do 1 i = 2, n
         im = i - 1
         pout(i) = pin(i)
         pdif = pout(i) - pout(im)
         ndif = pdif*cinv + sign(.5,pdif)
         pout(i) = pout(i) - ndif*const
         de = e(i) - e(im)
         a(i) = (pout(i)-pout(im))/de*cinv
         if ( dbgphs ) print 1000, i, pin(i),pout(i),ndif,a(i)
         if ( de .gt. esml ) go to 1
         isml = i
         esml = de
 1    continue
c
         if ( nf .lt. 1 .or. sdf .lt. 0.d0 ) go to 8
      if ( dbgphs ) print 1002
 1002 format ( " phasit -- smoothing the data" )
      nt = 2*nf
      do 2 i = 2, n
         im = i - 1
         i1 = max(2, i-nf)
         i2 = i1 + nt
         if ( i2 .le. n ) go to 3
         i2 = n
         i1 = max(2, n-nt)
 3       ntot = i2 - i1 + 1
         call stest ( a(i),a(i1),ntot, sdf, iflg, avg, sdev )
c
         if ( iflg .eq. 0 ) go to 5
         xns = cinv*(pout(i) - pout(im)) - avg*(e(i) - e(im))
         xn = xns + sign(.5,xns)
         nc = xn
         if ( nc .ne. 0 ) go to 6
         dif1 = a(i) - avg
         nct = sign(1.01d0,xn)
         atest = (pout(i)-pout(im)-nct*const)*cinv/(e(i)-e(im))
         dif2 = atest - avg
         if ( abs(dif2/dif1) .gt. 1.5d0 ) go to 5
         nc = nct
 6       continue
         padd = nc*const
         if ( dbgphs ) print 1003, padd, xn, i
c
         do 4 j = 1, im
            pout(j) = pout(j) + padd
 4       continue
c
         a(i) = (pout(i) - pout(im))/(e(i) - e(im))*cinv
 5       continue
         if ( dbgphs ) print 1005, i, pin(i), pout(i), a(i), avg, sdev
 2    continue
 8    continue
c
      if ( dbgphs ) print 1004, isml
      i = 1
      if ( dbgphs ) print 1005, i, e(1), pout(1)
 1004 format ( " phasit -- final derivative check", /,
     >            "     assume derivative correct at i = ", i5 )
c
c     choose closest derivative working backwards
c
      ismlm = isml - 1
      if ( ismlm .lt. 2 ) go to 17
c
      do 15 ir = 2, ismlm
         i = ismlm - ir + 2
         im = i - 1
         ip = i + 1
         aok = a(ip)
c
         pdif = pout(i) - pout(im)
         einv = cinv/(e(i) - e(im))
c
         atp = einv*(pdif - const)
         atm = einv*(pdif + const)
         if ( dbgphs ) print 1005, i, e(i), pout(i), atm, a(i), atp
         ad1 = abs(aok - atp)
         ad2 = abs(aok - a(i))
         ad3 = abs(aok - atm)
         adm = min(ad1,ad2,ad3)
c
         if ( adm .eq. ad2 ) go to 15
         padd = -const
         if ( adm .eq. ad1 ) padd = const
         call vas ( im, pout, 1, padd )
         a(i) = einv*(pout(i) - pout(im))
         if ( dbgphs ) print 1006, padd, i
 15   continue
 17   continue
c
c     choose closest derivative working forwards
c
      ismlp = isml + 1
      if ( ismlp .gt. n ) go to 21
c
      do 20 i = ismlp, n
         im = i - 1
         aok = a(im)
c
         pdif = pout(i) - pout(im)
         einv = cinv/(e(i) - e(im))
c
         atp = einv*(pdif - const)
         atm = einv*(pdif + const)
         if ( dbgphs ) print 1005, i, e(i), pout(i), atm, a(i), atp
         ad1 = abs(aok - atp)
         ad2 = abs(aok - a(i))
         ad3 = abs(aok - atm)
         adm = min(ad1,ad2,ad3)
c
         if ( adm .eq. ad2 ) go to 20
         padd = -const
         if ( adm .eq. ad1 ) padd = const
         call vas ( im, pout, 1, padd )
         a(i) = einv*(pout(i) - pout(im))
         if ( dbgphs ) print 1006, padd, i
 20   continue
 21   continue
c
c
 1005 format ( 1x,i5, 2f12.2, 3g12.4 )
 1006 format ( "   adding ", g12.4, " at i = ", i5 )
      nc = cinv*pout(n) + sign(.5,pout(n))
      padd = -nc*const
         call vas ( n, pout, 1, padd )
c
      return
 1000 format ( 1x,i5, 2f10.2, i5, g12.4 )
 1003 format ( "   adding ", 2g12.4, " to data at i = ", i5 )
      end
