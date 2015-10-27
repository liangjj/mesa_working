*deck m6004
c***begin prologue     m6004
c***date written       890404   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m6004, link 6004, spline
c***author             rescigno, t. n.(llnl)
c***source             m6004
c***purpose            spline fit of bessel functions
c***description        bessel functions calculated and then converted
c***                   to appropriate cutoff functions which are then
c***                   spline fit, coefficients put on file. this is a
c***                   stand alone code.  
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       m6004
      program fitbes
      implicit integer (a-z)
      parameter (acc=30)
      character *8 cpass, bessfl, chrkey, filkne
      character *4096 ops
      character *800 card
      logical logkey, prnbes, prnspn, green
      real *8 z, fpkey, alpha, rmin, rmax, gamma, rdel, rd26
      real *8 rl
      dimension z(1)
      common a(1)
      common /io/ inp, iout
      common /memory/ ioff
      equivalence(a,z)
c----------------------------------------------------------------------c
c  lbig = biggest physical l value allowed                             c 
c  nblokk = # pts. in large argument region for splines - don't change c
c  nsblok = # pts. in small argument region for splines - don't change c
c  ltop = top l for backwards recursion (maximum value)                c
c----------------------------------------------------------------------c
      call drum
      call iosys ('read character "kohn data filename" from rwf',0,0,0,
     1             filkne)
      call iosys ('read character options from rwf',-1,0,0,ops)
      call iosys ('open kohndt as old',0,0,0,filkne)
      prnbes=logkey(ops,'print=m6004=bessel',.false.,' ')
      prnspn=logkey(ops,'print=m6004=spline',.false.,' ')
      green=logkey(ops,'m6004=no-greens-function',.false.,' ')
      call posinp('$bess',cpass)
      call cardin(card)
      bessfl=chrkey(card,'bessel-file-name','bessel',' ')
      call locase(bessfl,bessfl)
      call iosys ('open bessel as new on ssd',262144,0,0,bessfl)
      call iosys ('write character "bessel function filename" to rwf',
     1             0,0,0,bessfl)
c---------------------------------------------------------------------c
c   note that rmin is the first spline knot point. if the real        c
c    grid goes below rmin, the cubic spline is being extrapolated     c
c                         (should be o.k.)                            c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                      generate spline data                           c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c             llmax = max l for splined bessel fcns                   c
c             alpha is for inner cut-off (y fcns)                     c
c             gamma and ncut are for outer cut-off                    c
c             rmax(rmin) = max(min) argument for                      c
c                          splined bessel fcns                        c
c             nper is the number of points to use per unit interval   c
c                         in computing the spline coefficients        c
c---------------------------------------------------------------------c  
      write (iout,100)
      lmax=intkey(card,'maximum-l-for-bessel',10,' ')
      alpha=fpkey(card,'inner-cutoff',1.0d0,' ')
      call iosys ('read real rmin from kohndt',1,rmin,0,' ')
      call iosys ('read real rmax from kohndt',1,rmax,0,' ')
      gamma=fpkey(card,'outer-exponential-cutoff',1.0d0,' ')
      ncut=intkey(card,'outer-n-cutoff',6,' ')
      nper=intkey(card,'points-per-interval-for-spline-coefficients',
     1            4,' ')      
      ltop=intkey(card,'top-recursion-l',100,' ')       
c---------------------------------------------------------------------c
c            estimate starting l                                      c
c---------------------------------------------------------------------c
      rl=lmax*acc
      strtl=lmax+sqrt(rl)
      strtl=max(strtl,lmax+1)
      ltop=min(ltop,strtl)
c-----------------------------------------------------------------------c
c compute total number of points,nr, and number of points in backward   c
c                       recursion call, ns.                             c
c-----------------------------------------------------------------------c
      nr=(rmax-rmin)*nper
      rdel=(rmax-rmin)/(nr-1)
      rd26=rdel*rdel/6.
      ns=lmax/rdel+1
      if(ns.gt.nr)ns=nr
      write(iout,10) rmin,rmax,rdel,lmax,alpha,gamma,ncut
      write (iout,20) nr, ns, ltop
      dimf1=nr*(ltop+1)
      dimf2=nr*(lmax+1)
      words=iadtwp(lmax+1)+6*nr+4*dimf1+8*dimf2
      write (iout,30) words
      call iosys ('read integer maxsiz from rwf',1,canget,0,' ')
      if (words.gt.canget) then
          call lnkerr ('cannot get required memory:will quit')
      endif
      call iosys ('write integer maxsiz to rwf',1,words,0,' ')
      call getscm(words,z,ngot,'m6004',0)
      ipow=wpadti(ioff)
      x=iadtwp(ipow+lmax+1)
      xinv=x+nr
      j=xinv+nr
      jp=j+dimf1
      y=jp+dimf1
      yp=y+dimf1
      norm=yp+dimf1
      toim=norm
      toim2=toim+nr
      hs=toim2+nr
      hsder=hs+2*dimf2
      cj=hsder+2*dimf2
      cy=cj+2*dimf2
      d=cy+2*dimf2
      call mkpow(a(ipow),lmax)
c----------------------------------------------------------------------c
c                 make grid                                            c
c----------------------------------------------------------------------c
      call mkx(z(x),z(xinv),rmin,rdel,nr)
c----------------------------------------------------------------------c
c                write out necessary stuff to bessel                   c
c----------------------------------------------------------------------c      
      call iosys ('write integer maximum-l to bessel',1,lmax,0,' ')
      call iosys ('write integer "total pts" to bessel',1,nr,0,' ')
      call iosys ('write real points to bessel',nr,z(x),0,' ')
      call iosys ('write integer "no. small pts" to bessel',1,ns,
     1            0,' ')
      call iosys ('write real rmin to bessel',1,rmin,0,' ')
      call iosys ('write real rmax to bessel',1,rmax,0,' ')
      call iosys ('write real "r spacing" to bessel',1,rdel,0,' ')
      call iosys ('write real "inner cutoff" to bessel',1,alpha,0,' ')
      call iosys ('write real "outer exp cutoff" to bessel',1,
     1            gamma,0,' ')
      call iosys ('write integer "outer n cutoff" to bessel',1,
     1            ncut,0,' ')
c----------------------------------------------------------------------c
c             small r, backward recursion                              c
c----------------------------------------------------------------------c
      call rcbesb(z(x),z(xinv),z(j),z(jp),z(y),z(yp),z(norm),ns,nr,lmax,
     1            ltop,prnbes)
c----------------------------------------------------------------------c
c             large r, forward recursion                               c
c----------------------------------------------------------------------c
      if (nr.gt.ns) then
          nptlft=nr-ns
          call rcbesf(z(x+ns),z(xinv+ns),z(j+ns),z(jp+ns),z(y+ns),
     1                z(yp+ns),nptlft,nr,lmax,ltop,prnbes)
      endif
c----------------------------------------------------------------------c
c              make the spline coefficients                            c
c----------------------------------------------------------------------c
      call mkspln(z(x),z(xinv),z(j),z(jp),z(y),z(yp),z(hs),z(hsder),
     1            z(cj),z(cy),z(d),z(toim),z(toim2),rdel,alpha,gamma,
     2            ncut,a(ipow),nr,lmax,ltop)
      nowds=2*dimf2
      if (prnspn) then
          call prspln(z(x),z(hs),z(hsder),z(cj),z(cy),nr,lmax)
      endif
      call iosys('write real hs to bessel',nowds,z(hs),0,' ')            
      call iosys('write real hsder to bessel',nowds,z(hsder),0,' ')            
      call iosys('write real cj to bessel',nowds,z(cj),0,' ')            
      call iosys('write real cy to bessel',nowds,z(cy),0,' ')            
   10 format(/,5x,'rmin spline',1x,f10.5,1x,'rmax spline',1x,f10.5,1x,
     1            'step size',1x,f10.5,/,30x,'max l',1x,i3,/,5x,
     2            'cutoff parameters:     alpha',1x,f6.3,1x,'gamma',1x,
     3            f6.3,1x,'ncut',1x,i3)
   20 format(/,5x,'total no. points',1x,i5,1x,/,5x,
     1       'no. points small r region',1x,i4,/,5x,
     2       'top recursion l',1x,i4)
   30 format (/,5x,'get',1x,i8,1x,'words of memory')
  100 format(//,20x,'***** m6004:splined ricatti-bessel functions *****'
     1)
      call iosys ('write integer maxsiz to rwf',1,canget,0,' ')
      call chainx(0)
      stop
      end
      subroutine mkpow(ipow,lmax)
      implicit integer (a-z)
      dimension ipow(0:lmax)
      data zero, one/ 0, 1/
      ipow(0)=zero
      do 10 i=1,lmax
         ipow(i)=one
   10 continue
      return
      end
*deck mkspln
c***begin prologue     mkspln
c***date written       xxxxxx   (yymmdd)
c***revision date      890422   (yymmdd)
c***keywords           m6004, link 6004, bessel, spline
c***author             rescigno, t. n.(llnl)
c***source             m6004
c***purpose            spline fit free functions
c***description        regular and irregular functions spline fit
c***references       
c
c***routines called
c***end prologue       mkspln
      subroutine mkspln(x,xinv,j,jp,y,yp,hs,hsder,cj,cy,d,toim,toim2,
     1                  rdel,alpha,gamma,ncut,ipow,nr,lmax,ltop)
      implicit integer (a-z)
      real *8 x, alpha, gamma, toim, toim2, term, term2, cr, crp
      real *8 crpp, gr, grp, grpp, gc, gcp, gcpp, j, jp, y, yp
      real *8 rdel, xinv, rdeli
      complex *16 hs, hsder, yp1, yp2, ai, precon, cj, d, cy
      dimension x(nr), xinv(nr), toim(nr), toim2(nr), j(nr,0:ltop)
      dimension jp(nr,0:ltop), y(nr,0:ltop), yp(nr,0:ltop)
      dimension hs(nr,0:lmax), hsder(nr,0:lmax), cj(nr,0:lmax)
      dimension cy(nr,0:lmax), d(nr), ipow(0:lmax)
      ai=cmplx(0.d0,1.d0)
      precon = -0.5d0*ai
      rdeli=1.d+00/rdel
      do 10 i=1,nr
         toim(i) = exp(-alpha*x(i))
         toim2(i)=exp(-gamma*x(i))
   10 continue
      do 20 l=0,lmax
         do 30 i=1,nr
            term = 1.0d0-toim(i)
            term2=1.0d0-toim2(i)
            cr = term**(2*l+3)
            crp = (2*l+3)*alpha*toim(i)*term**(2*l+2)
            crpp = crp*((2*l+2)*alpha*toim(i)/term - alpha)
            gr = term2**ncut
            grp=ncut*toim2(i)*gamma*term2**(ncut-1)
            grpp=grp*((ncut-1)*gamma*toim2(i)/term2 - gamma)
            gc=cr*gr
            gcp=cr*grp+crp*gr
            gcpp= cr*grpp+2.d0*grp*crp + gr*crpp
            hs(i,l) = gr*(ai*j(i,l) -cr*y(i,l))*xinv(i)
            hsder(i,l)=precon*(j(i,l)*grpp + 2.d0*grp*jp(i,l) +
     1                ai*(y(i,l)*gcpp + 2.d0*gcp*yp(i,l)))*
     2                                  xinv(i)**ipow(l)
   30    continue
         yp1=(hs(2,l)-hs(1,l))*rdeli
         yp2 = (hs(nr,l)-hs(nr-1,l))*rdeli
         call spline(x,hs(1,l),nr,yp1,yp2,d,cj(1,l))
         yp1=(hsder(2,l)-hsder(1,l))*rdeli
         yp2 = (hsder(nr,l)-hsder(nr-1,l))*rdeli
         call spline(x,hsder(1,l),nr,yp1,yp2,d,cy(1,l))
   20 continue
      return
      end
      subroutine mkx(x,xinv,rmin,rdel,n)
      implicit integer (a-z)
      real *8 x, xinv, rmin, rdel, one
      dimension x(n), xinv(n)
      data one / 1.0d+00 /
      do 10 i=1,n
         x(i)=rmin+(i-1)*rdel
         xinv(i)=one/x(i)
   10 continue
      return
      end
      subroutine prspln(x,hs,hsder,cj,cy,nr,lmax)
      implicit integer (a-z)
      common/ io/ inp, iout
      complex *16 hs, hsder, cj, cy
      real *8 x
      dimension x(nr), hs(nr,0:lmax), hsder(nr,0:lmax)
      dimension cj(nr,0:lmax), cy(nr,0:lmax)
      do 10 l=0,lmax
         write (iout,100) l
         write (iout,200) (hs(i,l), i=1,nr)
         write (iout,110) l
         write (iout,200) (hsder(i,l), i=1,nr)
         write (iout,120) l
         write (iout,200) (cj(i,l), i=1,nr)
         write (iout,130) l
         write (iout,200) (cy(i,l), i=1,nr)
   10 continue
  100 format(/,5x,'hs for l',1x,i3)
  110 format(/,5x,'hsder for l',1x,i3)
  120 format(/,5x,'cj for l',1x,i3)
  130 format(/,5x,'cy for l',1x,i3)
  200 format(/,5x,4d15.8)
      return
      end
*deck rcbesb
c***begin prologue     rcbesb
c***date written       xxxxxx   (yymmdd)
c***revision date      890422   (yymmdd)
c***keywords           m6004, link 6004, bessel, spline
c***author             schneider, barry (lanl)
c***source             m6004
c***purpose            generate bessel functions
c***description        bessel functions calculated backward
c***                   recursion. if greater accuracy needed change
c***                   parameter statement.
c***references       
c
c***routines called
c***end prologue       rcbesb
      subroutine rcbesb(x,xinv,j,jp,y,yp,norm,np,ptdim,lmax,ltop,prnt)
      implicit integer (a-z)
      parameter (acc=30)
      common /io/ inp, iout
      real *8 one, two, x, j, jp, y, yp, xinv, norm, rl
      logical prnt
      dimension x(ptdim), xinv(ptdim), j(ptdim,0:ltop), jp(ptdim,0:ltop)
      dimension y(ptdim,0:ltop), yp(ptdim,0:ltop), norm(ptdim)
      data one, two/ 1.d+00, 2.d+00 /
c---------------------------------------------------------------------c
c            estimate starting l                                      c
c---------------------------------------------------------------------c
      rl=lmax*acc
      strtl=lmax+sqrt(rl)
      strtl=max(strtl,lfinal)
      if (strtl.gt.ltop) then
          call lnkerr('starting l bigger than ltop')
      endif
c----------------------------------------------------------------------c
c               make first two bessel functions                        c
c----------------------------------------------------------------------c
      do 10 i=1,np
         j(i,0)=sin(x(i))
         norm(i)=j(i,0)
         y(i,0)=-cos(x(i))
         jp(i,0)=-y(i,0)
         yp(i,0)=j(i,0)
   10 continue
      if (lmax.eq.0) then
          return
      endif
      do 20 i=1,np
         j(i,1)=j(i,0)*xinv(i)+y(i,0)
         y(i,1)=y(i,0)*xinv(i)-j(i,0)
         jp(i,1)=j(i,0)-one*j(i,1)*xinv(i)
         yp(i,1)=y(i,0)-one*y(i,1)*xinv(i)    
   20 continue
      if (lmax.eq.1) then
          return
      endif       
c----------------------------------------------------------------------c
c                   lmax is greater than one                           c
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c              calculate y by upward recursion                         c
c----------------------------------------------------------------------c
       lfinal=min(ltop,lmax+1)     
       do 30 l=1,lfinal-1
          ll1=l+l+1
          lm1=l-1
          lp1=l+1
          do 40 i=1,np
             y(i,lp1)=ll1*y(i,l)*xinv(i)-y(i,lm1)
   40     continue
   30  continue
c----------------------------------------------------------------------c
c             calculate j by downward recursion                        c
c----------------------------------------------------------------------c
       call rzero(j(1,strtl),np)
       onelss=strtl-1
       do 50 i=1,np
          j(i,onelss)=one
   50  continue
       do 60 l=onelss-1,0,-1
          ll3=l+l+3
          lp1=l+1
          lp2=l+2
          do 70 i=1,np
             j(i,l)=ll3*j(i,lp1)*xinv(i)-j(i,lp2)
   70     continue
   60  continue
c----------------------------------------------------------------------c
c                normalize the j                                       c
c----------------------------------------------------------------------c
       do 80 i=1,np
          norm(i)=norm(i)/j(i,0)
  80   continue
       do 90 l=0,lfinal
          do 100 i=1,np
             j(i,l)=j(i,l)*norm(i)
  100     continue
   90  continue
c---------------------------------------------------------------------c
c             finish calculation by getting derivatives               c
c---------------------------------------------------------------------c
       do 200 l=0,lmax
          lp1=l+1
          do 300 i=1,np
             jp(i,l)=lp1*j(i,l)*xinv(i)-j(i,lp1)            
             yp(i,l)=lp1*y(i,l)*xinv(i)-y(i,lp1) 
  300    continue
  200 continue
      if (prnt) then
          do 800 l=0,lmax
             write (iout,900) l
             write (iout,1000) (j(i,l),i=1,np)
             write (iout,1010) l
             write (iout,1000) (jp(i,l),i=1,np)
             write (iout,1020) l
             write (iout,1000) (y(i,l),i=1,np)
             write (iout,1030) l
             write (iout,1000) (yp(i,l),i=1,np)
  800     continue
      endif
  900 format(/,5x,'regular function l = ',1x,i3)
 1000 format( (/,5x,5(d15.8,1x) ) )
 1010 format(/,5x,'derivative regular function l = ',1x,i3)
 1020 format(/,5x,'irregular function l = ',1x,i3)
 1030 format(/,5x,'derivative irregular function l = ',1x,i3)
      return
      end
*deck rcbesf
c***begin prologue     rcbesf
c***date written       xxxxxx   (yymmdd)
c***revision date      890422   (yymmdd)
c***keywords           m6004, link 6004, bessel, spline
c***author             schneider, barry (lanl)
c***source             m6004
c***purpose            generate bessel functions
c***description        bessel functions calculated by forward
c***                   recursion. use only when lmax lt arg
c***references       
c
c***routines called
c***end prologue       rcbesf
      subroutine rcbesf(x,xinv,j,jp,y,yp,nr,ptdim,lmax,ltop,prnt) 
      implicit integer (a-z)
      common /io/ inp, iout
      logical prnt
      real *8 x, j, jp, y, yp, xinv
      dimension j(ptdim,0:lmax), jp(ptdim,0:lmax), x(*), xinv(*)
      dimension y(ptdim,0:lmax), yp(ptdim,0:lmax)
      do 10 i=1,nr
         j(i,0)=sin(x(i))
         y(i,0)=-cos(x(i))
         jp(i,0)=-y(i,0)
         yp(i,0)=j(i,0)
   10 continue
      if (lmax.eq.0) then
          return
      else
          do 20 i=1,nr
             y(i,1)=y(i,0)*xinv(i)-j(i,0)
             j(i,1)=j(i,0)*xinv(i)+y(i,0)
   20     continue
          if (lmax.ne.1) then
              do 30 l=1,lmax-1
                 lp1=l+1
                 ll1=l+l+1
                 lm1=l-1
                 do 40 i=1,nr
                    j(i,lp1)=ll1*j(i,l)*xinv(i)-j(i,lm1)
                    y(i,lp1)=ll1*y(i,l)*xinv(i)-y(i,lm1)
   40            continue
   30         continue
          endif
          do  50 l=1,lmax
              lm1=l-1
              do 60 i=1,nr
                 jp(i,l)=j(i,lm1)-l*j(i,l)*xinv(i)
                 yp(i,l)=y(i,lm1)-l*y(i,l)*xinv(i)
   60         continue
   50     continue
      endif
      if (prnt) then
          do 800 l=0,lmax
             write (iout,900) l
             write (iout,1000) (j(i,l),i=1,nr)
             write (iout,1010) l
             write (iout,1000) (jp(i,l),i=1,nr)
             write (iout,1020) l
             write (iout,1000) (y(i,l),i=1,nr)
             write (iout,1030) l
             write (iout,1000) (yp(i,l),i=1,nr)
  800     continue
      endif
  900 format(/,5x,'regular function l = ',1x,i3)
 1000 format( (/,5x,5(d15.8,1x) ) )
 1010 format(/,5x,'derivative regular function l = ',1x,i3)
 1020 format(/,5x,'irregular function l = ',1x,i3)
 1030 format(/,5x,'derivative irregular function l = ',1x,i3)
      return
      end
*deck spline
c***begin prologue     spline
c***date written       890404   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           spline, link 6004, orbital decomposition
c***author             rescigno, t. n.(llnl)
c***source             spline
c***purpose            spline fit
c***references       
c
c***routines called    none
c***end prologue       spline
      subroutine spline(x,y,n,yp1,ypn,scr,y2)
      implicit integer (a-z)
      real *8 x
      complex*16 y, yp1, ypn, scr, y2, p, qn, un, sig
      dimension x(n), y(n), scr(n), y2(n)
      y2(1)=-.5d0
      scr(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      do 10 i=2,n-1
         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
         p=sig*y2(i-1)+2.d0
         y2(i)=(sig-1.d0)/p
         scr(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     1           /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*scr(i-1))/p
   10 continue
      qn=.5d0
      un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      y2(n)=(un-qn*scr(n-1))/(qn*y2(n-1)+1.)
      do 20 k=n-1,1,-1
         y2(k)=y2(k)*y2(k+1)+scr(k)
   20 continue
      return
      end
