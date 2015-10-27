*deck @(#)xyzout.f	4.1  7/7/93
      subroutine xyzout(phix,phiy,phiz,ngridx,ngridy,ngridz,
     $                 nat,npf,nmo,c,lowdin,coord,ian,
     $                 xmin,ymin,zmin,incrx,incry,incrz,
     $                 nbasis,charge,multip,nel,nalphe,nbetae)
      implicit integer(a-z)
      integer ian(nat)
      real*8 xmin,ymin,zmin,incrx,incry,incrz,coord(3,nat)
      real*8 phix(ngridx,npf),phiy(ngridy,npf),phiz(ngridz,npf)
      real*8 c(npf,nmo),lowdin(npf,nmo)
      real*8 zan
      character*128 plotnm
      character*4 ptgrp
      character*80 title
      data ptgrp/'c1  '/
      data title/'test'/
      data naxis/0/
      save ptgrp,title,naxis
c
c
 1000 format('mesa')
 1010 format(a80)
 1015 format(1x,'group =',a4,'  naxis =',i5,' natoms=',i5)
 1020 format(1x,'charge=',i5,' multip=',i5,' nbasis=',i5,
     $         ' nelect=',i5,' nbetae=',i5,' nalphe=',i5)
 1025 format(1x,i5,1x,f4.1,3(1x,e16.9))
 1028 format(1x,'ngridx=',i5,' ngridy=',i5,' ngridz=',i5)
 1029 format(1x,'xmin  =',e20.10,' ymin  =',e20.10,' zmin  =',e20.10)
 1030 format(1x,'incrx =',e20.10,' incry =',e20.10,' incrz =',e20.10)
 1031 format(1x,'nprim =',i5,' nmo  =',i5)
 1032 format(1x,'x amplitudes')
 1033 format(1x,'y amplitudes')
 1034 format(1x,'z amplitudes')
 1035 format(1x,'mo:',i5)
 1036 format(1x,'lowdin mo:',i5)
 1039 format(6f20.10)
      iplot=3
      plotnm='plot'
      call comand(1,plotnm)
      if (plotnm.eq.'plot') plotnm='mesa.plot'
c
      open (unit=iplot,file=plotnm,access='sequential',
     $      form='formatted',err=2000,status='unknown')
c
      write(iplot,1000)
      write(iplot,1010) title
      write(iplot,1015) ptgrp,naxis,nat
      write(iplot,1020) charge,multip,nbasis,nel,nbetae,nalphe
      do 10 i=1,nat
         zan=float(ian(i))
         write(iplot,1025) i,zan,(coord(j,i),j=1,3)
   10 continue
      write(iplot,1028) ngridx,ngridy,ngridz
      write(iplot,1029) xmin,ymin,zmin
      write(iplot,1030) incrx,incry,incrz
      write(iplot,1031) npf,nmo
      write(iplot,1032)
      do 100 i=1,ngridx
         write(iplot,1039) (phix(i,pf),pf=1,npf)
  100 continue
      write(iplot,1033)
      do 110 i=1,ngridy
         write(iplot,1039) (phiy(i,pf),pf=1,npf)
  110 continue
      write(iplot,1034)
      do 120 i=1,ngridz
         write(iplot,1039) (phiz(i,pf),pf=1,npf)
  120 continue
      do 150 j=1,nmo
         write(iplot,1035) j
         write(iplot,1039) (c(i,j),i=1,npf)
  150 continue
      do 200 j=1,nmo
         write(iplot,1036) j
         write(iplot,1039) (lowdin(i,j),i=1,npf)
  200 continue
c
      close (iplot)
      return
c
 2000 call lnkerr('m1990 cannot open the plot file.')
      return
      end
