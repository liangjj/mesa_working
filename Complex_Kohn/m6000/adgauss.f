*deck @(#)adgauss.f	1.1 9/9/91
      program adgauss 
      parameter ( nocen=10 , nptdim=100)
      implicit real*8 (a-h,o-z)
      character *4096 ops
      character *8 cpass, chrkey
      character *128 filgrd, filkne
      character *2 itoc
      character *1600 card
      integer pbuf, pntbuf
      logical itran,iplot,logkey,movech,nokohn
      common z(1)
      dimension ia(1)
      equivalence (z,ia)
      common /memory / ioff
      common /gauss/ x(nptdim,3),w(nptdim,3),scr(nptdim),dummy(2)
      common /centsi/ncent,iplot
      common /centsr/cc(3),a(3,nocen),eta(nocen),
     1              znuc(nocen),amu(nocen), dum(nocen)
      common/io/ inp,iout
      dimension cdum(3,nocen)
      data pi/3.14159265358d+00/
c
      call drum
      call iosys ('read character options from rwf',-1,0,0,ops)
      nokohn=logkey(ops,'m6000=no-kohndt',.false.,' ')
      write (iout,2)
      call posinp('$grid',cpass)
      call cardin(card)
      pbuf=intkey(card,'point-buffer',5000,' ')
      if (.not.nokohn) then
          call iosys ('read character "kohn data filename" from rwf',
     1                 -1,0,0,filkne)
          call iosys ('open kohndt as old',0,0,0,filkne)
          call iosys ('read integer "point buffer" from kohndt',1,
     1                 pbuf,0,' ')
      endif
      iplot=logkey(ops,'plots',.false.,' ')
      itran=logkey(ops,'no-co-ordinate-transform',.false.,' ')
      movech=logkey(ops,'to-center-of-charge',.false.,' ')
      call iosys ('read character "grid filename" from rwf',-1,0,0,
     1             filgrd)
c read centers and the parameters for the yukawa potential
c
      write(iout,432)
  432 format(' type of background grid is spheres ')
      zsum=0.d+00
      call posinp('$centers',cpass)
      call cardin(card)
      ncent=intkey(card,'no-centers',1,' ')
      do 400 i=1,ncent
      eta(i)=fpkey(card,'exponent-center-'//itoc(i),0.d+00,' ')
      znuc(i)=fpkey(card,'charge-center-'//itoc(i),1.d+00,' ')
      call fparr(card,'position-center-'//itoc(i),cdum(1,i),3,' ')
      amu(i)=fpkey(card,'test-exponent-center-'//itoc(i),0.d+00,' ')
      cc(i)=0.d+00
      zsum=zsum+znuc(i)
      write(iout,821) eta(i),(cdum(iii,i),iii=1,3),znuc(i)
  821 format(' yukawa exponent ',f10.5,'  center  ',3f10.5
     1, ' charge ',f10.5)
c
c  read the parameters for the transformation
c   amu = exponent
c
      write(iout,822) amu(i)
  400 continue
      do 55 i=1,ncent
      do 55 j=1,3
   55 a(j,i)=cdum(j,i)
822   format(' test exponent ',f10.5,
     1  /,' transformation centered at yukawa singularity')
c
c
c compute center of charge
c
      call rzero(cc,3)
      if (movech) then
          write (iout,77)
   77 format(/,10x,'move to center of charge co-ordinates')
          do 41 i=1,ncent
          do 41 j=1,3
41           cc(j)=cc(j)+znuc(i)*a(j,i)/zsum
c
c compute distance to farthest nucleus from center of charge
c
          dist=0.d+00
          write(iout,639)
639       format(/' distances from cc to nuclei :')
          do 42 i=1,ncent
             temp=(cc(1)-a(1,i))**2+(cc(2)-a(2,i))**2+(cc(3)-a(3,i))**2
             temp=sqrt(temp)
             write(iout,641)temp
641   format(f10.5)
42        if(dist.lt.temp)dist=temp
c
c move center of coordinates to cc
c
      else
          write (iout,78)
   78 format(/,10x,'do not move to center of charge co-ordinates')          
      endif
          do 43 i=1,ncent
          do 43 j=1,3
43           a(j,i)=a(j,i)-cc(j)
          write(iout,371)((a(i,j),i=1,3),j=1,ncent)
371   format(//,' nuclei in cc coordinates'
     1       ,/,(3(2x,f10.5)))
      write(iout,52)
52    format(/' spherical polar coordinates of nuclei:'/)
c     do 50 i=1,ncent
c     rr=sqrt(a(1,i)**2+a(2,i)**2+a(3,i)**2)
c     theta=acos(a(3,i)/rr)
c     ss=sin(theta)
c     phi=acos(a(1,i)/rr/ss)
c     theta=theta*180.d+00/pi
c     phi=phi*180.d+00/pi
c     write(iout,51)rr,theta,phi
50    continue
51    format(3f10.5)
c
      call setsph(ops,igrid,isave)
      if (isave.eq.0) then
          isave=1
      endif
      pntbuf=min(igrid,pbuf)
      call iosys ('write integer "point buffer" to rwf',1,pntbuf,0,' ')
c           get memory
      nwords=wpadti(4*(igrid+isave+1))
      call iosys ('read integer maxsiz from rwf',1,maxsiz,0,' ')
      if (nwords.gt.maxsiz) then
          call lnkerr ('insufficient memory')
      endif
      call iosys ('write integer maxsiz to rwf',1,nwords,0,' ')
      call getscm(nwords,z(1),ngot,'grid',0)
      igr=ioff
      isve=igr+4*igrid
      call shells(z(igr),z(isve),igrid,isave,ops)
c
c  do an integral and apply condensing transformation
c
      call iosys ('open grid as new on ssd',4*igrid,0,0,filgrd)
      call iosys ('write integer "point buffer" to grid',1,pntbuf,0,' ')
      call trnsf(z(igr),z(isve),sum,suma,sumb,sumc,igrid,isave,
     1           itran)
      exact=0.d+00
      do 57 i=1,ncent
      exact=exact+4.d+00*pi/eta(i)**2
   57 continue
      write(iout,939) exact
  939 format('        exact value of yukawa integral ',e15.8)
      write(iout, 936) suma
  936 format(' integral using untransformed coords   ',e15.8)
      write(iout,937) sumb
 937  format(' integral using transformed coords     ',e15.8)
      write(iout,951) sum
  951  format(' volume =',e15.8)
      write(iout,938)sumc
  938  format(' volume from transformed coords =',e15.8)
c
c plot grid in x y plane
c
c     if(iplot)then
c     call posinp('$plots,cpass)
c     call cardin(card)
c     xmn=fpkey(card,'x-min',0.d+00,' ')
c     xmx=fpkey(card,'x-max',0.d+00,' ')
c     ymn=fpkey(card,'y-min',0.d+00,' ')
c     ymx=fpkey(card,'y-max',0.d+00,' ')
c     read(inp,*)xmn,xmx,ymn,ymx
c     call maps(xmn,xmx,ymn,ymx)
c     if(igrid.le.isave) call points(save(1,2),save(1,4),igrid)
c     if(igrid.gt.isave) call points(save(1,2),save(1,4),isave)
c     call frame
c     call maps(xmn,xmx,ymn,ymx)
c     if(igrid.le.isave) call points(save(1,1),save(1,3),igrid)
c     if(igrid.gt.isave) call points(save(1,1),save(1,3),isave)
c     call frame
c     endif
c     call iosys ('rewind all on kohndt read-and-write',0,0,0,' ')
      call iosys ('rewind all on grid read-and-write',0,0,0,' ')
      call iosys ('close grid',0,0,0,' ')     
c     call iosys ('close kohndt',0,0,0,' ')
      call iosys ('write integer maxsiz to rwf',1,maxsiz,0,' ')
      call chainx(0)
    1 format(a80)
    2 format (//,15x,'***** m6000:grid generation program *****',//)
      stop
      end
