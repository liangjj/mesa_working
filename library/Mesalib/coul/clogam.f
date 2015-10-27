*deck clogam
      function clogam(z)
      common /io/ inp, iout
c
c     this routine computes the logarithm of the gamma function gamma(z)
c     for any complex argument 'z' to any accuracy preset by call logam
c
      complex z,u,v,h,r,clogam,cdigam,ser
      dimension b(15),bn(15),bd(15)
c
      data lerr /6/, nx0 /6/, nb /15/,
     x  zero,one,two,four,half,quart /0e+0,1e+0,2e+0,4e+0,.5e+0,.25e+0/
      data bn(1),bd(1)    / +1e+0,   6e+0 /,
     x     bn(2),bd(2)    / -1e+0,  30e+0 /,
     x     bn(3),bd(3)    / +1e+0,  42e+0 /,
     x     bn(4),bd(4)    / -1e+0,  30e+0 /,
     x     bn(5),bd(5)    / +5e+0,  66e+0 /,
     x     bn(6),bd(6)    /          -691e+0,  2730e+0/,
     x     bn(7),bd(7)    /          +  7e+0,     6e+0/,
     x     bn(8),bd(8)    /         -3617e+0,   510e+0/,
     x     bn(9),bd(9)    /         43867e+0,   798e+0/,
     x     bn(10),bd(10)  /       -174611e+0,   330e+0/,
     x     bn(11),bd(11)  /        854513e+0,   138e+0/,
     x     bn(12),bd(12)  /    -236364091e+0,  2730e+0/,
     x     bn(13),bd(13)  /     + 8553103e+0,     6e+0/,
     x     bn(14),bd(14)  /  -23749461029e+0,   870e+0/,
     x     bn(15),bd(15)  / 8615841276005e+0, 14322e+0/
      data fplmin / -850e+0 /
c
      x=real(z)
      t=aimag(z)
      mx = int(real(accur*100 - x))
      if(abs(abs(x)-mx) + abs(t).lt.accur*100) go to 60
      f=abs(t)
      v=cmplx(x,f)
      if(x .lt. zero) v=one-v
      h=zero
      c=real(v)
      n=nx0-int(c)
      if(n .lt. 0) go to 30
      h=v
      d=aimag(v)
      a=atan2(d,c)
      if(n .eq. 0) go to 20
      do 10 i = 1,n
      c=c+one
      v=cmplx(c,d)
      h=h*v
   10 a=a+atan2(d,c)
   20 h=cmplx(half*log(real(h)**2+aimag(h)**2),a)
      v=v+one
   30 r=one/v**2
      ser = b(nt)
      do 40 j=2,nt
        k = nt+1 - j
   40 ser = b(k) + r*ser
      clogam = hl2p+(v-half)*log(v)-v + ser/v - h
      if(x .ge. zero) go to 50
c
      a= int(x)-one
      c=pi*(x-a)
      d=pi*f
c     e=exp(-two*d)
        e = zero
        f = -two*d
        if(f.gt.fplmin) e = exp(f)
      f=sin(c)
      e= d + half*log(e*f**2+quart*(one-e)**2)
      f=atan2(cos(c)*tanh(d),f)-a*pi
      clogam=alpi-cmplx(e,f)-clogam
c
   50 if(sign(one,t) .lt. -half) clogam=conjg(clogam)
      return
c
   60 write(iout,1000) 'clogam',x
 1000 format(1x,a6,' ... argument is non positive integer = ',f20.2)
      clogam = zero
      return
c
      entry cdigam(z)
c
c     this routine computes the logarithmic derivative of the gamma
c     function  psi(z) = digamma(z) = d (ln gamma(z))/dz  for any
c     complex argument z, to any accuracy preset by call logam(acc)
c
      u=z
      x=real(u)
      a=abs(x)
      if(abs(aimag(u)) + abs(a + int(x)) .lt. accur) go to 110
      if(x .lt. zero) u=-u
      v=u
      h=zero
      n=nx0-int(a)
      if(n .lt. 0) go to 90
      h=one/v
      if(n .eq. 0) go to 80
      do 70 i = 1,n
      v=v+one
   70 h=h+one/v
   80 v=v+one
   90 r=one/v**2
      ser = b(nt) * (2*nt-1)
      do 100 j=2,nt
        k = nt+1 - j
  100 ser = b(k)*(2*k-1) + r*ser
      cdigam = log(v) - half/v - r*ser - h
      if(x .ge. zero) return
      h=pi*u
      cdigam = cdigam + one/u + pi*cos(h)/sin(h)
      return
c
  110 write(iout,1000) 'cdigam',x
      cdigam=zero
      return
c
      entry logam(acc)
c
c      initialisation call for calculations to accuracy 'acc'
c
      nx0 = 6
      x0  = nx0 + one
      pi = four*atan(one)
      alpi = log(pi)
      hl2p = log(two*pi) * half
      accur = acc
      do 120 k=1,nb
       f21 = k*2 - one
       b(k) = bn(k) / (bd(k) * k*two * f21)
       err = abs(b(k)) * k*two / x0**f21
  120 if(err.lt.acc) go to 130
       nx0 = int((err/acc)**(one/f21) * x0)
       k = nb
  130 nt = k
c     print *,' logam requires k = ',k ,' with cutoff at x =',nx0+1
      return
      end
