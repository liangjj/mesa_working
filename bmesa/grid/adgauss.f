c $Header: adgauss.f,v 1.1 92/12/31 15:23:07 bis Exp $
*deck @(#)adgauss.f	1.1 9/9/91
      program adgauss 
      parameter ( nocen=10 , nptdim=100 , numshl=50 )
      implicit real*8 (a-h,o-z)
      character *4096 ops
      character *8 cpass
      character *128 fillam
      character *2 itoc
      character *1600 card
      character*24 chrkey, intgrd
      integer pntbuf
      logical logkey, movech, posinp, prnt, atgrid
      common z(1)
      dimension ia(1)
      equivalence (z,ia)
      common /memory / ioff
      common /gauss/ x(nptdim,3), w(nptdim,3), scr(nptdim), dummy(2)
      common /centsi/ ncent
      common /centsr/ cc(3), a(3,nocen), eta(nocen),
     1                znuc(nocen)
      common /spheri/ nshell, nr(numshl), nthell(numshl), 
     1                nphell(numshl), ntheta(numshl,numshl),
     2                nphi(numshl,numshl)
      common /spherr/ r(numshl), theta(numshl,numshl),
     1                phi(numshl,numshl)
      common /quadpt/ nsub(numshl), nth(numshl), nph(numshl)
      common/quadtp/ itype, jtype, ktype
      common/io/ inp,iout
      dimension cdum(3,nocen)
      data pi/3.14159265358d+00/
c
      call drum
      call iosys ('read character options from rwf',-1,0,0,ops)
      atgrid=logkey(ops,'generate-atomic-grids',.false.,' ')
      pntbuf=intkey(ops,'point-buffer',10000,' ')
      movech=logkey(ops,'to-center-of-charge',.false.,' ')
      write (iout,100)
      pntbuf=10000
      call iosys ('read character "linear algebraic filename" from rwf',
     1                 -1,0,0,fillam)
      call iosys ('open lamdat as new',0,0,0,fillam)
c read centers and the parameters for the yukawa potential
      zsum=0.d+00
      if ( posinp('$centers',cpass) ) then
           call cardin(card)
           ncent=intkey(card,'no-centers',1,' ')
           intgrd=chrkey(card,'integrand-type',
     1                               'shielded-coulomb', ' ')
           do 10 i=1,ncent
              eta(i)=fpkey(card,'exponent-center-'//itoc(i),0.d+00,' ')
              znuc(i)=fpkey(card,'charge-center-'//itoc(i),1.d+00,' ')
              call fparr(card,'position-center-'//itoc(i),cdum(1,i),
     1                   3,' ')
              cc(i)=0.d+00
              zsum=zsum+znuc(i)
              write(iout,110) eta(i),(cdum(iii,i),iii=1,3),znuc(i)
   10      continue
      endif
      call iosys ('create real "center positions" on lamdat',
     1             3*ncent,0,0,' ')
      call iosys ('write integer "number of centers" to lamdat',1,
     1             ncent,0,' ')
      call iosys ('write real "yukawa exponents" to lamdat',ncent,
     1             eta,0,' ')
      call iosys ('write real "nuclear charges" to lamdat',ncent,
     1             znuc,0,' ')
      do 20 i=1,ncent
         call iosys ('write real "center positions" to lamdat '//
     1               'without rewinding',3,cdum(1,i),0,' ')
         do 30 j=1,3
            a(j,i)=cdum(j,i)
   30    continue
   20 continue
c compute center of charge
      call rzero(cc,3)
      if (movech) then
          write (iout,120)
          do 40 i=1,ncent
             do 50 j=1,3
                cc(j)=cc(j)+znuc(i)*a(j,i)/zsum
   50        continue
   40     continue
c
c compute distance to farthest nucleus from center of charge
c
          dist=0.d+00
          write(iout,130)
          do 60 i=1,ncent
             temp=(cc(1)-a(1,i))**2+(cc(2)-a(2,i))**2+(cc(3)-a(3,i))**2
             temp=sqrt(temp)
             write(iout,140) temp
             if(dist.lt.temp) then
                dist=temp
             endif
   60     continue
c
c move center of coordinates to cc
c
      else
          write (iout,150)
      endif
      call iosys ('write real "center of charge" to lamdat',3,cc,0,' ')
      call iosys ('create real "center positions in center of charge '//
     1            'coordinates" on lamdat',3*ncent,0,0,' ')
      do 70 i=1,ncent
         do 80 j=1,3
            a(j,i)=a(j,i)-cc(j)
   80    continue
         call iosys ('write real "center positions in center of '//
     1               'charge coordinates" to lamdat without rewinding',
     2                3,a(1,i),0,' ')  
   70 continue    
      write(iout,160) (( a(i,j),i=1,3), j=1,ncent )
      write(iout,170)
c     generate main grid and perform integral tests
      write(iout,*) '    generate main integration grid'
      call setsph(ops)
      nrmax=0
      nthmax=0
      nphmax=0
      ntot=0
      do 90 i=1,nshell
         nrmax=max(nrmax,nr(i))
         nthmax=max(nthmax,nth(i))
         nphmax=max(nphmax,nph(i))
         ntot=ntot+nsub(i)
   90 continue
      pntbuf=min(ntot,pntbuf)
      call iosys ('write integer "point buffer" to lamdat',1,
     1             pntbuf,0,' ')
c           get memory
      call getscm(need,z(1),maxcor,'grid',0)
      i1=ioff
      i2=i1+nrmax*nshell
      i3=i2+nthmax*nshell
      i4=i3+nphmax*nshell
      i5=i4+nrmax*nshell
      i6=i5+nthmax*nshell
      i7=i6+nphmax*nshell
      prnt=logkey(ops,'print=m6200',.false.,' ')
      call shells(z(i1),z(i2),z(i3),z(i4),z(i5),z(i6),nshell,nrmax,
     1            nthmax,nphmax,'all',prnt)
c
c     do an integral numerically and analytically
c
      call trnsf(z(i1),z(i2),z(i3),z(i4),z(i5),z(i6),z(i7),volume,
     1           yukawa,intgrd,nshell,nrmax,nthmax,nphmax,'all')
      exact=0.d+00
      if (intgrd.eq.'shielded-coulomb') then
          do 200 i=1,ncent
             exact=exact+4.d+00*pi/eta(i)**2
  200     continue
      elseif (intgrd.eq.'exponential') then
          do 300 i=1,ncent
             exact=exact+8.d+00*pi/eta(i)**3
  300     continue
      endif
      write(iout,180) volume
      write(iout,190) exact
      write(iout,191) yukawa
c     if atomic grids are needed generate them and the model atomic
c                                potentials
      if (atgrid) then
          write(iout,*) '    generate atomic grids'
          do 500 i=1,ncent
             write(iout,*)
             write(iout,*) '     generating grid for atom = ',i
             call satshl(i)
             nrmax=0
             nthmax=0
             nphmax=0
             do 95 j=1,nshell
                nrmax=max(nrmax,nr(j))
                nthmax=max(nthmax,nth(j))
                nphmax=max(nphmax,nph(j))
                ntot=ntot+nsub(j)
   95        continue
             j1=ioff
             j2=j1+nrmax*nshell
             j3=j2+nthmax*nshell
             j4=j3+nphmax*nshell
             j5=j4+nrmax*nshell
             j6=j5+nthmax*nshell
             j7=j6+nphmax*nshell
             call shells(z(j1),z(j2),z(j3),z(j4),z(j5),z(j6),nshell,
     1                   nrmax,nthmax,nphmax,itoc(i),prnt)
             call trnsf(z(j1),z(j2),z(j3),z(j4),z(j5),z(j6),z(j7),
     1                  volume,yukawa,intgrd,nshell,nrmax,nthmax,
     2                  nphmax,itoc(i))
             if (intgrd.eq.'shielded-coulomb') then
                 exact=4.d+00*pi/eta(i)**2
             elseif (intgrd.eq.'exponential') then
                 exact=8.d+00*pi/eta(i)**3
             endif
             write(iout,190) exact
             write(iout,191) yukawa
  500     continue
      endif
      call chainx(0)
  100 format (//,15x,'***** m6000:grid generation program *****',//)
  110 format(' yukawa exponent ',f10.5,'  center  ',3f10.5
     1, ' charge ',f10.5)
  120 format(/,10x,'move to center of charge co-ordinates')
  130 format(/' distances from cc to nuclei :')
  140 format(f10.5)
  150 format(/,10x,'do not move to center of charge co-ordinates')         
  160 format(//,' nuclei in cc coordinates',/,(3(2x,f10.5)))
  170 format(/' spherical polar coordinates of nuclei:'/)
  180 format(/,5x,'volume enclosed by integration region = ',e15.8)
  190 format(/,5x,'exact     value of yukawa integral = ',e15.8)
  191 format(/,5x,'numerical value of yukawa integral = ',e15.8)
      stop
      end
