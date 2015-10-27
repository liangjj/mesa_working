*deck  @(#)drvintl.f	2.1 10/10/91
      subroutine drvint(z,a)
c***begin prologue     drvint
c***date written       920324   (yymmdd)
c***revision date               (yymmdd)
c
c
c***keywords
c***author             schneider, barry  (nsf)
c***source             %W% %G% 
c***purpose
c                      fft for gaussian integrals 
c***description
c***references
c
c***routines called
c
c***end prologue       kintgl
c
      implicit integer (a-z)
c
      parameter ( bfmx=100, cnmx=10, qdmx=1000)
      integer a(1)
      real*8 z(1), alpha, center, rcen, kpt, kwt, pi
      real *8 fpkey, kstr, kfnr, kfs
      dimension alpha(bfmx), center(3,bfmx), rcen(3,cnmx)
      dimension kpt(qdmx), kwt(qdmx)
c     
      character*3 itoc
      character *4096 ops
      character *10 cpass
      character *1600 card
      logical posinp
      common/io/inp,iout
      data pi /3.14159265358979323846d+00/
c
      ncen=1
      call iosys ('read character options from rwf',-1,0,0,ops)      
      if ( posinp('$orbs',cpass) ) then
           call cardin(card)
           ncen=intkey(card,'number-centers',1,' ')
           nrg=intkey(card,'number-integration-regions',1,' ')
           do 10 i=1,ncen
              call fparr(card,'center-'//itoc(i),rcen(1,i),3,' ')
   10      continue     
      endif
      count=0
      do 20 i=1,ncen
         if( posinp('$basis-'//itoc(i),cpass) ) then
             call cardin(card)
             numb=intkey(card,'number-basis-functions',1,' ')
             call fparr(card,'exponents',alpha(count+1),numb,' ')
             do 30 j=1,numb
                count=count+1
                center(1,count)=rcen(1,i)
                center(2,count)=rcen(2,i)
                center(3,count)=rcen(3,i)
   30        continue
         endif
   20 continue     
      nbf=count
      if (nbf.eq.0) then
          call lnkerr('error in basis set number')
      endif
      write(iout,100) nbf
      ntri=nbf*(nbf+1)/2
      npts=0
      do 40 iz=1,nrg
         if ( posinp('$region-'//itoc(iz),cpass) ) then
              call cardin(card)
              mshtyp=chrkey(card,'type-quad','legendre',' ')
              npnk=intkey (card,'number-points',3,' ')
              kstr=fpkey (card,'starting-k',-50.d0,' ')
              kfnr=fpkey (card,'ending-k',50.d0,' ')
              kfs=kfnr-kstr
              do 50 jz=1,npnk
                 npts=npts+1
                 call lgndrx (npnk,jz,kwt(npts),kpt(npts))
                 kwt(npts)=kwt(npts)*kfs
                 kpt(npts)=kpt(npts)*kfs+kstr
   50         continue
         endif
   40 continue
      write(iout,200) npts
      alf12=1
      cen12=alf12+ntri
      pre12=cen12+ntri
      ftran=pre12+ntri
      words=ftran+ntri*npts*npts*npts
c**********************************************************************c
c                 get new centers, prefactors and exponents            c
c**********************************************************************c
      call prod(alpha,center,z(alf12),z(cen12),z(pre12),nbf,ntri)
c**********************************************************************c
c             calculate fourier transform of density                   c
c**********************************************************************c
      call gtrans(z(alf12),z(cen12),z(pre12),kpt,kpt,kpt,z(ftran),ntri,
     1            npts)    
c**********************************************************************c
c                     do the integrals                                 c
c**********************************************************************c
      call chainx(0)
c
c
  100 format(//,10x,'number basis functions',1x,i5)
  200 format(//,10x,'number of k-space integration points in each dimens
     1ion',1x,i4) 
      stop
      end
