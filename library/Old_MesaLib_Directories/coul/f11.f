*deck f11
      function f11(x,eta,zl,p,eps,limit,kind,err,nits,fpmax,acc8,acc16)
      common /io/ inp, iout
      complex x,eta,zl,p,aa,bb,z,f11,cdigam,ci
       complex dd,g,f,ai,bi,t
      logical zllin
      real*16 ar,br,gr,gi,dr,di,tr,ti,ur,ui,fi,fi1,den
      data zero,one,two / 0e+0, 1e+0, 2e+0 /, ci / (0e+0, 1e+0) /
      absc(aa) = abs(real(aa)) + abs(aimag(aa))
      nintc(aa) = nint(real(aa))
c
c *** evaluate the hypergeometric function 1f1
c                                        i
c            f (aa;bb; z) = sum  (aa)   z / ( (bb)  i| )
c           1 1              i       i            i
c
c     to accuracy eps with at most limit terms.
c  if kind = 0 : using extended precision but real arithmetic only,
c            1 : using normal precision in complex arithmetic,
c   or       2 : using normal complex arithmetic, but with cdigam factor
c
c  where
         aa = zl+one - eta*p
         bb = two*(zl+one)
c  and
         z  = two*p*x
c
         zllin = real(bb).le.zero .and. abs(bb-nintc(bb)).lt.acc8**0.25
             if(.not.zllin.or.real(bb)+limit.lt.1.5) go to 10
                nits = -1
                return
   10 if(limit.le.0) then
         f11 = zero
         err = zero
         nits= 1
         return
         endif
      ta = one
      rk = one
      if(kind.le.0.and.absc(z)*absc(aa).gt.absc(bb) * 1.0) then
         dr = one
         di = zero
         gr = one
         gi = zero
         ar = real(aa)
         br = real(bb)
         fi = zero
      do 20 i=2,limit
         fi1 = fi + one
         tr = br * fi1
         ti = aimag(bb) * fi1
         den= one / (tr*tr + ti*ti)
         ur = (ar*tr + aimag(aa)*ti) * den
         ui = (aimag(aa)*tr - ar*ti) * den
         tr = ur*gr - ui*gi
         ti = ur*gi + ui*gr
         gr = real(z) * tr - aimag(z)*ti
         gi = real(z) * ti + aimag(z)*tr
         dr = dr + gr
         di = di + gi
            err = abs(gr) + abs(gi)
               if(err.gt.fpmax) go to 60
            rk  = abs(dr) + abs(di)
            ta = max(ta,rk)
         if(err.lt.rk*eps .or. i.ge.4.and.err.lt.acc16) go to 30
         fi = fi1
         ar = ar + one
   20    br = br + one
c
   30    f11 = dr + ci * di
         err = acc16 * ta / rk
c
      else
c* ---------------------------------- alternative code
c*    if real*16 arithmetic is not available, (or already using it|),
c*    then use kind > 0
         g = one
          f = one
          if(kind.ge.2) f = cdigam(aa) - cdigam(bb) - cdigam(g)
         dd = f
         do 40 i=2,limit
            ai = aa + (i-2)
            bi = bb + (i-2)
            r  = i-one
         g = g * z * ai / (bi * r)
         if(kind.ge.2)
c                              multiply by (psi(a+r)-psi(b+r)-psi(1+r))
     x        f = f + one/ai - one/bi - one/r
         t  = g * f
         dd = dd + t
            err = absc(t)
               if(err.gt.fpmax) go to 60
            rk = absc(dd)
         ta = max(ta,rk)
         if(err.lt.rk*eps.or.err.lt.acc8.and.i.ge.4) go to 50
   40    continue
 
   50    err = acc8 * ta / rk
         f11 = dd
c* ------------------------------------------- end of alternative code
      endif
   60    nits = i
      return
      end
