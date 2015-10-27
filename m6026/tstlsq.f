*deck tstlst
      program tstlst
      implicit integer (a-z)
      character*1600 card
      character*16 ftyp, chrkey
      character*8 cpass
      real*8 z, fpkey, rmin, rmax, start, del
      character*4096 ops
      dimension z(1), power(0:20)
      common a(1)
      common /io/ inp, iout
      common / memory / ioff
      equivalence (z,a)
      call drum
      call iosys ('read character options from rwf',-1,0,0,ops)
      call posinp('$tstlst',cpass)
      call cardin(card)
      npts=intkey(card,'no-points',100,' ')
      rmin=fpkey(card,'rmin',0.d0,' ')
      rmax=fpkey(card,'rmax',1.d0,' ')
      ftyp=chrkey(card,'function-type','exponential',' ')
      npar=intkey(card,'no-linear-parameters',2,' ')
      do 5 i=0,npar-1
         power(i)=i
 5    continue
      nwords=2*npts+2*npar+npar*npar+npar*npts
      call getscm(nwords,z,ngot,'tstlst',0)
      f=ioff
      pt=f+npts
      rhs=pt+npts
      xn=rhs+npar
      coef=xn+npar*npts
      ipvt=wpadti(coef+npar*npar)
      icnt=f
      jcnt=pt
      start=rmin
      del=(rmax-rmin)/(npts-1)
      do 10 i=1,npts
         z(icnt)=exp(start)
         z(jcnt)=start
         icnt=icnt+1
         jcnt=jcnt+1
         start=start+del
   10 continue
      call plyfit(z(f),z(coef),z(rhs),z(xn),z(pt),power,a(ipvt),
     1            npar,npts,.true.,'first')
      call plyfit(z(f),z(coef),z(rhs),z(xn),z(pt),power,a(ipvt),
     1            npar,npts,.true.,'second')
      call chainx(0)
      stop
      end
