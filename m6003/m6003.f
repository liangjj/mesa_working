*deck m6003
c***begin prologue     m6003
c***date written       881006   (yymmdd)
c***revision date      890418   (yymmdd)
c***keywords           m6003, link 6003, spherical harmonics
c***author             rescigno, t. n.(llnl)
c***source             m6003
c***purpose            spherical harmonic generation
c***description        calculates spherical harmonics on grid
c***                   for integral calculation. fully vectorized
c*** 
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       m6003
      program ylm 
      implicit integer (a-z)
      character *3 ans
      character *8 cpass, itp, chrkey
      character *1600 card
      character *4096 ops
      character *13 grdnam
      logical grdtyp, logkey, prnt
      real *8 z
      dimension itp(2), z(1)
      common a(1)
      common /io/ inp, iout
      common /memory / ioff
      equivalence (z,a)
      call drum
      call iosys ('read character options from rwf',-1,0,0,ops)
      prnt=logkey(ops,'print=m6003=ylm',.false.,' ')
      write (iout,100)
      call posinp ('$ylm',cpass)
      call cardin(card)
      grdtyp=logkey(card,'untransformed-grid',.false.,' ')
      grdnam='"trns grid"'
      if (grdtyp) then
          grdnam='"untrns grid"'
      endif
      lmax=intkey(card,'maximum-l-value',30,' ')
      mumax=intkey(card,'maximum-m-value',5,' ')
      write (iout,300) lmax, mumax
      call iosys ('does "grid filename" exist on rwf',0,0,0,ans)
      if (ans.ne.'no') then
          call iosys ('read character "grid filename" from rwf',-1,0,0,
     1                 itp(1))
      else
          itp(1)=chrkey(card,'grid-file-name','grid',' ')
      endif
      call iosys ('does "spherical harmonic filename" exist on rwf',0,
     1             0,0,ans)
      if (ans.ne.'no') then
          call iosys ('read character "spherical harmonic filename" '//
     1                'from rwf',-1,0,0,itp(2))
      else
          itp(2)=chrkey(card,'ylm-file-name','ylms',' ')
      endif
      call locase(itp(2),itp(2))
      call iosys ('open grid as old',0,0,0,itp(1))
      call iosys ('read integer "no. grid pts" from grid',1,ngrid,0,
     1            ' ')
      call iosys ('read integer "point buffer" from grid',1,
     1            pntbuf,0,' ')
      write (iout,400) ngrid, pntbuf
      call iosys ('open ylms as new on ssd',262144,0,0,itp(2))
      call iosys ('write character "grid type" to ylms',0,0,0,grdnam)
      call iosys ('write character "spherical harmonic filename" to '//
     1            'rwf',0,0,0,itp(2))
      call iosys ('write integer "no. grid pts" to ylms',1,ngrid,0,
     1            ' ')
      call iosys ('write integer "point buffer" to ylms',1,pntbuf,0,
     1            ' ')
      call iosys ('write integer "max l in ylm" to ylms',1,lmax,0,' ')
      call iosys ('write integer "max m in ylm" to ylms',1,mumax,0,' ')
      call iosys ('create real ylm on ylms',-1,0,0,' ')
c----------------------------------------------------------------------c
c              calculate number of passes needed                       c
c----------------------------------------------------------------------c
      npass=ngrid/pntbuf
      nwleft=ngrid-pntbuf*npass
      if (nwleft.ne.0) then
          npass=npass+1
      else
          nwleft=pntbuf
      endif
      lmmax=lmax+mumax
c----------------------------------------------------------------------c
c                     memory allocation                                c
c----------------------------------------------------------------------c
      call iosys ('read integer maxsiz from rwf',1,canget,0,' ')
      maxfac=lmmax+10
      words=2*maxfac+11*pntbuf+3*pntbuf*(lmax+1)
      if (canget.lt.words) then
          call lnkerr('not enough memory:decrease point buffer')
      endif
      call iosys ('write integer maxsiz to rwf',1,words,0,' ')
      call getscm(words,z,ngot,'m6003',0)      
      dfct=ioff
      ddfct=dfct+maxfac
      grid=ddfct+maxfac
      x=grid+4*pntbuf
      cphi=x+pntbuf
      sphi=cphi+pntbuf
      plm=sphi+pntbuf
      cmp=plm+pntbuf*(lmax+1)
      cmphi=cmp+pntbuf+pntbuf
      smphi=cmphi+pntbuf
      ylmp=smphi+pntbuf
      ylmm=ylmp+pntbuf*(lmax+1)
      words=ylmm+pntbuf*(lmax+1)
 
      call fact(z(dfct),z(ddfct),maxfac)
c----------------------------------------------------------------------c
c                    loop over points                                  c
c----------------------------------------------------------------------c
      do 10 ipass=1,npass
         nwread=pntbuf
         if (ipass.eq.npass) then
             nwread=nwleft
         endif
c----------------------------------------------------------------------c
c                    read buffer of points                             c
c----------------------------------------------------------------------c        
         call iosys ('read real '//grdnam//' from grid without '//
     1               'rewinding',4*nwread,z(grid),0,' ')
c----------------------------------------------------------------------c
c                 calculate some necessary functions                   c
c                 used over and over again in plm routine              c
c----------------------------------------------------------------------c 
         call miscfn(z(grid),z(x),z(cphi),z(sphi),nwread)
c----------------------------------------------------------------------c
c             calculate the legendre functions                         c
c----------------------------------------------------------------------c
         do 20 mu=0,mumax
            call legend(z(plm),z(x),z(dfct),z(ddfct),nwread,lmax,mu,
     1                  maxfac)
            call sph(z(plm),z(cphi),z(sphi),z(cmp),z(cmphi),
     1               z(smphi),z(ylmp),z(ylmm),z(dfct),nwread,mu,lmax,
     2               maxfac,ipass)
            if (prnt) then
                call prnty(z(plm),z(x),nwread,mu,lmax)
            endif
   20    continue
   10 continue
c----------------------------------------------------------------------c
c              close random file                                       c
c----------------------------------------------------------------------c
      call iosys ('endfile ylm on ylms',0,0,0,' ')
      call iosys ('rewind all on ylms read-and-write',0,0,0,' ')
      call iosys ('rewind all on grid read-and-write',0,0,0,' ')
      call iosys('close ylms',0,0,0,' ')
      call iosys ('close grid',0,0,0,' ')
      call iosys ('write integer maxsiz to rwf',1,canget,0,' ')
      call chainx(0)
      stop
  100 format (//,20x,'***** m6003: calculate spherical harmonics *****')
  300 format(/,5x,'maximum l value',1x,i3,5x,'maximum m value',1x,i3)
  400 format(/,5x,'no. points',1x,i8,5x,'point buffer',1x,i8)
      end
      subroutine miscfn(grid,x,cphi,sphi,npnt)
      implicit integer (a-z)
      real *8 grid, r, x, cphi, sphi
      dimension grid(4,npnt), x(npnt), cphi(npnt)
      dimension sphi(npnt)
      do 10 pnt=1,npnt
         r=sqrt(grid(1,pnt)*grid(1,pnt)+grid(2,pnt)*grid(2,pnt)+
     1             grid(3,pnt)*grid(3,pnt))
         x(pnt)=grid(3,pnt)/r
         cphi(pnt)=grid(1,pnt)/(r*sqrt(1.d+00-x(pnt)*x(pnt)))
         sphi(pnt)=grid(2,pnt)/(r*sqrt(1.d+00-x(pnt)*x(pnt)))
   10 continue
      return
      end
      subroutine sph(plm,cphi,sphi,cmp,cmphi,smphi,ylmp,ylmm,dfct,n,
     1               mu,lmax,maxfac,pass)
      implicit integer (a-z)
      common /io/ inp, iout
      dimension plm(n,0:lmax), ylmp(n,0:lmax), ylmm(n,0:lmax)
      dimension dfct(0:maxfac)
      dimension cphi(n), sphi(n), cmp(n), cmphi(n), smphi(n)
      real *8 plm, phi, cphi, sphi, cmphi, smphi, const, ylmp, ylmm
      real *8 sq2, dfct, pi, fourpi
      complex*16 cmp
      data pi, fourpi /3.14159265358979323846d+00,
     1                12.56637061435917295384d+00/
      sq2=sqrt(2.d+00)
      if(mu.eq.0) then
         do 10 l=0,lmax
            const=sqrt((2*l+1)/fourpi)
            do 20 i=1,n
               ylmp(i,l)=plm(i,l)*const
               plm(i,l)=ylmp(i,l)
   20       continue
   10    continue
         words=n*(lmax+1)
         call iosys ('write real ylm to ylms without rewinding',words,
     1               ylmp,0,' ')
         write (iout,100) pass,mu,words
         return
      else
         do 30 i=1,n
            cmp(i)=(cmplx(cphi(i),sphi(i)))**mu
            cmphi(i)=real(cmp(i))
            smphi(i)=aimag(cmp(i))
   30    continue
         do 40 l=mu,lmax
         const=sqrt((2*l+1)/fourpi*dfct(l-mu)/dfct(l+mu))
c  sqrt(2) factor added to normalize "real valued" ylms
         const=const*sq2
            do 50 i=1,n
               ylmp(i,l)=plm(i,l)*cmphi(i)*const
               ylmm(i,l)=plm(i,l)*smphi(i)*const
   50       continue
   40    continue
         words=n*(lmax+1)
         call iosys ('write real ylm to ylms without rewinding',words,
     1                  ylmp(1,0),0,' ')
         call iosys ('write real ylm to ylms without rewinding',words,
     1                  ylmm(1,0),0,' ')
         words=words+words
         write (iout,100) pass,mu,words
      endif
  100 format(/,5x,'pass',1x,i3,2x,'m value',1x,i3,2x,'words written',1x,
     1             i8)
      return
      end
      subroutine prnty(plm,x,n,mu,lmax)
      dimension x(n), plm(n,0:lmax)
      common /io/ inp, iout         
      write (iout,100) mu
      do 200 i=1,n
         write (iout,300) x(i)
         write (iout,400) (plm(i,j),j=mu,lmax)
  200 continue
  100 format(/,5x,'plm for mu',1x,i2)
  300 format(/,5x,'arg',1x,e15.8)
  400 format( (/,5x,5e15.8) )
      return
      end     
*deck legend
c***begin prologue     legend
c***date written       880721   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           legend, link 2702, legendre functions
c***author             schneider, barry (lanl)
c***source             m2702
c***purpose            legendre functions
c***description        calculation of p(l,m) functions
c***references         none
c
c***routines called
c***end prologue       legend
      subroutine legend (plm,x,dfct,ddfct,npt,lmax,m,maxfac)
      implicit integer (a-z)
      real *8 plm, x, dfct, ddfct, fm, facx, f1, anorm
      dimension plm(npt,0:lmax), x(npt)
      dimension  dfct(0:maxfac), ddfct(0:maxfac)
c----------------------------------------------------------------------c
c           start recursion with plm(m,m) and plm(m+1,m)               c
c                      and recur upward                                c
c----------------------------------------------------------------------c
      do 10 i=m,lmax
         do 10 j=1,npt
            plm(j,i)=0.d+00
   10 continue
      fm=.5d+00*m
      do 20 i=1,npt
         facx=1.d+00
         if (fm.ne.0.d+00) then
             facx=(1.d+00-x(i)*x(i))**fm
         endif
         plm(i,m)=ddfct(m)*facx
   20 continue
      if (lmax.ne.m) then
          mm=m+m+1
          mpls1=m+1
          do 30 i=1,npt
             plm(i,mpls1)=mm*x(i)*plm(i,m)
   30     continue
          if (lmax.ne.mpls1) then
              lind=m+2
              n1=2
              n2=m+m+3
              n3=n2-2
              do 50 i=lind,lmax
                 ii=i-1
                 jj=i-2
                 do 40 j=1,npt
                    f1=n2*x(j)*plm(j,ii)-n3*plm(j,jj)
                    f1=f1/n1
                    plm(j,i)=f1
   40            continue
                 n1=n1+1
                 n2=n2+2
                 n3=n3+1
   50       continue
          endif
      endif
      return
c
      end
*deck fact
c***begin prologue     fact
c***date written       880721   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           fact, link 2702, factorials
c***author             schneider, barry (lanl)
c***source             m2702
c***purpose            factorials
c***description        calculation of factorials
c***references         none
c
c***routines called
c***end prologue       fact
      subroutine fact(dfct,ddfct,maxfac)
      implicit integer (a-z)
      real *8 dfct, ddfct
      dimension dfct(0:maxfac), ddfct(0:maxfac)
c----------------------------------------------------------------------c
c               calculate factorials                                   c
c----------------------------------------------------------------------c
      dfct(0)=1.d+00
      dfct(1)=1.d+00
      if (maxfac.gt.1) then
          do 10 i=2,maxfac
             dfct(i)=i*dfct(i-1)
   10     continue
      endif
c----------------------------------------------------------------------c
c           calculate (2*m-1) double factorial                         c
c----------------------------------------------------------------------c
      ddfct(0)=1.d+00
      ddfct(1)=1.d+00
      ddfct(2)=3.d+00
      if (maxfac.gt.2) then
         do 20 i=3,maxfac
            ddfct(i)=(i+i-1)*ddfct(i-1)
   20    continue  
      endif
      return
      end
