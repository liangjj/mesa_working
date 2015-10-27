h56460
s 00232/00000/00000
d D 1.1 94/02/16 20:35:06 mesa 1 0
c date and time created 94/02/16 20:35:06 by mesa
e
u
U
f e 0
t
T
I 1
*deck %W%  %G%
      subroutine pm6200(z,ia)
c***begin prologue     %M%
c***date written       930802   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           gnquad, link m6200
c***author             schneider, barry (nsf)
c***source             m6200
c***purpose            calculate atomic grids and quadrature weights
c***references         
c
c***routines called    satshl, shells, voronoi
c***end prologue       %M%
      implicit real*8 (a-h,o-z)
      parameter ( nocen=10 , numshl=50 )
      character*4096 ops
      character*30 cpass, str
      character*128 fillam
      character*3 itoc, chra, yn
      character*1600 card
      integer wpadti,cskipb
      logical logkey, positn, prnt, voron, noscat, nonsep
      real*8 z(1)
      integer ia(1)
      common/io/ inp,iout
      dimension a(3,nocen), eta(nocen), znuc(nocen)
      dimension nr(numshl), r(numshl+1), ngrid(nocen), nwts(nocen)
      data pi/3.14159265358d+00/
c
      call iosys ('read character options from rwf',-1,0,0,ops)
      niter=intkey(ops,'number-of-voronoi-iterations',3,' ')
      alpha=fpkey(ops,'exponential-scattering-grid-cutoff',1.d0,' ')
      prnt=logkey(ops,'print=m6200=grid',.false.,' ')
      voron=logkey(ops,'voronoi=off',.false.,' ')
      cutoff=fpkey(ops,'scattering-grid-cutoff',10.d0,' ')
      noscat=logkey(ops,'no-scattering-center',.false.,' ')
      yn='yes'
      if (noscat) then
          yn='true'
      endif    
      write (iout,100)
      call iosys ('read character "linear algebraic filename" from rwf',
     1                 -1,0,0,fillam)
      call iosys ('open lamdat as new',0,0,0,fillam)
      call iosys('write character "scattering center" to lamdat',0,0,
     1            0,yn)
c read centers and the parameters for the yukawa potential
      if ( positn('$centers',card,inp) ) then
           call cardin(card)
           ncent=intkey(card,'no-atomic-centers',1,' ')
           do 10 i=1,ncent
              eta(i)=fpkey(card,'exponent-center-'//itoc(i),0.d+00,' ')
              znuc(i)=fpkey(card,'charge-center-'//itoc(i),1.d+00,' ')
              call fparr(card,'position-center-'//itoc(i),a(1,i),
     1                   3,' ')
              write(iout,110) eta(i),(a(iii,i),iii=1,3),znuc(i)
   10      continue
           ncplus=ncent
c scattering center coordinates are read in here. the default is (0.,0.,0.)
           if (.not.noscat) then
                ncplus=ncent+1         
                a(1,ncplus)=0.d0
                a(2,ncplus)=0.d0
                a(3,ncplus)=0.d0 
                call fparr(card,'scattering-center-position',
     1                     a(1,ncplus),3,' ')
                write(iout,115) (a(iii,ncplus),iii=1,3)
                call iosys ('write real "scattering center position" '//
     1                      'to lamdat',3,
     2                       a(1,ncplus),0,' ')
           endif     
      endif
      call iosys ('create real "atomic center positions" on lamdat',
     1             3*ncent,0,0,' ')
      call iosys ('write integer "number of atomic centers" to lamdat',
     1            1,ncent,0,' ')
      call iosys ('write real "yukawa exponents" to lamdat',ncent,
     1             eta,0,' ')
      call iosys ('write real "nuclear charges" to lamdat',ncent,
     1             znuc,0,' ')
      do 20 i=1,ncent
         call iosys ('write real "atomic center positions" to lamdat '//
     1               'without rewinding',3,a(1,i),0,' ')
   20 continue
      nrmax=0
      nthmax=0
      nphmax=0
      maxgrd=0
      maxshl=0
      ntotal=0
      nwttot=0
      ltop=0
      mtop=0
      ntrad=0
      maxang=0
      do 500 i=1,ncplus 
         if (i.le.ncent) then
             cpass='atom'
             write(iout,*)
             write(iout,*) '     reading grid data for atom = ',i
         else
             write(iout,*)
             write(iout,*) '     reading grid data for scattering '//
     1                     'center'
             cpass='scattering'
         endif
         call satshl(cpass,i,nr,r,numshl,ngrid(i),nwts(i),
     1               nrmax,nthmax,nphmax,maxgrd,maxshl,ntrad,
     2               ltop,mtop,nleb,nang,nonsep)
         maxang=max(maxang,nang)
         ntotal=ntotal+ngrid(i)
         nwttot=nwttot+nwts(i)
  500 continue
      call iosys ('write integer "max l value" to lamdat',1,ltop,0,' ')
      call iosys ('write integer "max m value" to lamdat',1,mtop,0,' ')
      call iosys('write integer "max radial points in shell" to '//
     1           'lamdat',1,nrmax,0,' ')
      call iosys('write integer "max radial points" to lamdat',1,
     1            ntrad,0,' ')
      call iosys ('write integer "max theta points" to lamdat',1,
     1             nthmax,0,' ')
      call iosys ('write integer "max phi points" to lamdat',1,
     1             nphmax,0,' ')
      call iosys ('write integer "max grid points" to lamdat',1,
     1             maxgrd,0,' ')
      call iosys ('write integer "biggest l" to lamdat',1,
     1             ltop,0,' ')   
      call iosys ('write integer "biggest m" to lamdat',1,
     1             mtop,0,' ')
      if (nonsep) then
          call iosys ('write integer "max number lebedev points" '//
     1                'to lamdat',1,nang,0,' ')          
      endif                
      write(iout,190) nrmax, nthmax, nphmax
      nwords=nrmax*maxshl+3*nthmax+4*nphmax+nrmax*(nrmax-1)*maxshl+
     1          maxgrd+6*ntotal+nwttot+max(nrmax,nthmax,nphmax,ncplus)+
     2          ncplus*ncplus
      if (nonsep) then
          nwords=nwords+3*maxang
      endif    
      nwords=wpadti(nwords)
      write(iout,*) 'get ',nwords,' integer words of core'
      call iosys ('read integer maxsiz from rwf',1,nget,0,' ')
      if (nwords.gt.nget) then
          call lnkerr('not enough memory. will quit')
      endif
      call iosys ('write integer maxsiz to rwf',1,nwords,0,' ')
      call getscm(nwords,z,ngot,'grid',0)
      j1=ioff
      j2=j1+nrmax*maxshl
      j3=j2+max(nthmax,maxgrd)
      j4=j3+nphmax
      j5=j4+nrmax*(nrmax-1)*maxshl
      j6=j5+nthmax
      j7=j6+nphmax
      j8=j7+nthmax
      j9=j8+nphmax
      j10=j9+nphmax
      j11=j10+maxgrd
      j12=j11+max(nrmax,nthmax,nphmax,ncplus)
      j13=j12+ncplus*ncplus
      j14=j13+3*ntotal
      j15=j14+nwttot
      j16=j15+3*ntotal
      jhold=j13
      khold=j14
      lhold=j15
      do 600 i=1,ncplus
         if (i.le.ncent) then
             cpass='atom'
             chra=itoc(i)
             len=cskipb(chra,' ')
             str='atom-'//chra(1:len)
         else
             str='scattering center'
             cpass='scattering'
         endif         
         call iosys ('read integer "number of shells '//
     1                str//'" from lamdat',1,nshell,0,' ')
         call iosys ('read integer "number radial points per '//
     1               'shell '//str//'" from lamdat',
     2                nshell,nr,0,' ')
         if (.not.nonsep) then   
             call iosys('read integer "theta quadrature order '//
     1                   str//'" from lamdat',1,nthet,0,' ')  
             call iosys('read integer "phi quadrature order '//
     1                   str//'" from lamdat',1,nphi,0,' ')
         else
             call iosys ('read integer "number lebedev points '//
     1                    str//'" from lamdat',1,nang,0,' ') 
             nphi=nang
             nthet=nang
         endif         
c         
c        we are actually storing all of the atomic grids since we
c        need them for the voroni transformation 
c        
         call shells(cpass,i,a,r,nr,eta,z(j1),z(j2),z(j3),z(j4),
     1               z(j5),z(j6),z(j16),z(j5),z(j7),z(j8),z(j9),z(j10),
     2               z(j11),z(j13),z(j14),z(j15),nshell,nrmax,nthet,
     3               nphi,ngrid(i),nwts(i),ncent,ncplus,prnt,
     4               nleb,nang,nonsep)
         j13=j13+3*ngrid(i)
         j14=j14+nwts(i)
         j15=j15+3*ngrid(i)
  600 continue
      j13=jhold
      j14=khold
      if (.not.voron) then
          write(iout,*) 'compute voronoi weights'
          call voronoi (eta,a,znuc,z(j13),z(j14),z(j10),z(j1),z(j12),
     1                  z(j11),z(j2),alpha,cutoff,ncent,ncplus,ngrid,
     2                  nwts,niter,nr,nrmax,maxshl,noscat,
     3                  nonsep,prnt)  
      else
          write(iout,*) 'do not compute voronoi weights'
      endif
      call iosys('rewind all on lamdat read-and-write',0,0,0,' ')
      call chainx(0)
  100 format (//,15x,'***** m6200:grid generation program *****',//)
  110 format(' yukawa exponent = ',f10.5,' center = (',f10.5,',',f10.5,
     1       ',',f10.5')'/,' charge = ',f10.5)
  115 format(' scattering center ',3f10.5)
  140 format(f10.5)
  170 format(/' spherical polar coordinates of nuclei:'/)
  180 format(/,'volume enclosed by integration region = ',e15.8)
  190 format(/,'maximum radial points in a shell ',i4,
     1         1x,'maximum theta points ',i4,/,
     2         'maximum phi points ',i4)
  195 format(/,'total number of grid points in this atom ',i6)
      return
      end
E 1
