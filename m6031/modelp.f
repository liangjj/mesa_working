*deck modelp
c***begin prologue     modelp 
c***date written       921701   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6005, link 6005, model potential
c***author             schneider, barry (nsf)
c***source             m6005
c***purpose            put model potentials on numerical integration grid 
c***description        current subroutine will handle square-well,
c***                   exponential and yukawa potentials. they are assumed
c***                   to be attractive.
c***references         none
c
c***routines called    none
c***end prologue       modelp 
      subroutine modelp(v,grid,grdtyp,type,nchan,nreg,pntbuf,nolst)
      implicit integer (a-z)
      parameter (ncmx=9)
      character *1600 card
      character *8 cpass
      character *(*) type, grdtyp
      character *1 itoc
      real *8 grid, v, boxsiz, fpkey, vcple
      dimension v(*), grid(4,*), vcple(ncmx,ncmx)
      ntri=nchan*(nchan+1)/2
      call posinp('$modpot',cpass)
      call cardin(card)
      if (type.eq.'square-well') then
          boxsiz=fpkey(card,'well-size',1.d+00,' ')
	  do 10 i=1,nchan
	     do 20 j=1,i
	        call fparr(card,'v'//itoc(i)//itoc(j),vcple(i,j),1,' ')
   20        continue
   10     continue
	  do 30 i=1,nchan
	     do 40 j=1,i
		vcple(j,i)=vcple(i,j)
   40        continue
   30     continue
      elseif(type.eq.'exponential'.or.type.eq.'yukawa') then
	  do 50 i=1,nchan
	     do 60 j=1,i
                call fparr(card,'alpha'//itoc(i)//itoc(j),
     1			   vcple(i,j),1,' ')
   60        continue
   50     continue
	  do 70 i=1,nchan
	     do 80 j=1,i
		vcple(j,i)=vcple(i,j)
   80        continue
   70     continue
      else
	  call lnkerr('potential type '//type// 'not currently handled')
      endif
      npnts=pntbuf
      do 90 ireg=1,nreg
	 if(ireg.eq.nreg) then
	    npnts=nolst
         endif
	 call iosys('read real '//grdtyp//' from grid without rewinding',
     1               4*npnts,grid,0,' ')
	 if (type.eq.'square-well') then
	     call well(v,vcple,grid,boxsiz,nchan,npnts,ntri)
         elseif (type.eq.'exponential') then
	     call expot(v,vcple,grid,nchan,npnts,ntri)
         elseif (type.eq.'yukawa') then
	     call yukawa(v,vcple,grid,nchan,npnts,ntri)
         endif
	 words=ntri*npnts
	 call iosys ('write real "static potential" to vstat without '//
     1               'rewinding',words,v,0,' ')
   90 continue
      call iosys ('rewind all on vstat read-and-write',0,0,0,' ')
      call iosys ('rewind all on grid read-and-write',0,0,0,' ')
      return
      end
