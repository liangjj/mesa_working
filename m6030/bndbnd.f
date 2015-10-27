*deck m6030
c***begin prologue     m6030
c***date written       920209   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6030 ,link 6030
c***author             schneider, barry (nsf)
c***source             m6030
c***purpose            simple bound-bound matrix elements
c***description        calculate bound-bound matrix elements of hamiltonian
c***                   for kohn calculations in simple model problems.
c***references
c
c***routines called    iosys, util and mdutil
c***end prologue       m6030
      program bndbnd
      implicit integer (a-z)
      parameter ( dimc=30 , dimbf=200 )
      real*8 z
      common a(1)
      dimension z(1)
      dimension ngauss(dimc), ngch(dimbf,dimc), nmoc(dimc)
      dimension nmoch(dimbf,dimc)
      equivalence (z,a)
      common /memory/ ioff
      common /io / inp, iout
      logical logkey, posinp, bsym
      logical prntcn, prntmo, unitmx, modpot
      character *4096 ops
      character *1600 card
      character *32 xform
      character *3 sym
      character *16 type 
      character *13 chrkey
      character *128 filgrd, filorb, filkne, filpot
      character *10 cpass
      character *13 grdtyp
      character *3 itoc
      character *80 title
      call drum
      call iosys ('read character options from rwf',-1,0,0,ops)
      prntcn=logkey(ops,'print=m6030=contracted',.false.,' ')
      prntmo=logkey(ops,'print=m6030=molecular',.false.,' ')
      unitmx=logkey(ops,'unit-matrix',.false.,' ')
      modpot=logkey(ops,'model-potential',.false.,' ')
      if (modpot) then
          type=chrkey(ops,'model-potential','square-well',' ')
      endif
c----------------------------------------------------------------------c
c             open the grid, orbital and kohn data files               c
c----------------------------------------------------------------------c
      call iosys ('read character "grid filename" from rwf',-1,0,0,
     1             filgrd)
      call iosys ('open grid as old',0,0,0,filgrd)
      call iosys ('read character "orbital filename" from rwf',-1,0,0,
     1             filorb)
      call iosys ('open orbs as old',0,0,0,filorb)
      call iosys ('read character "grid type" from orbs',0,0,0,grdtyp)
      call iosys ('read integer "no. grid pts" from orbs',1,npts,0,
     1            ' ')
      call iosys ('read integer "point buffer" from orbs',1,pntbuf,
     1            0,' ')
      call iosys ('read integer "no. cont" from orbs',1,ncon,0,0)
      call iosys ('read integer "no. regions" from orbs',1,nreg,0,' ')
      call iosys ('read integer "final pts" from orbs',1,nolst,0,' ')
      call iosys ('read character "kohn data filename" from rwf',-1,
     1              0,0,filkne)
      call iosys ('open kohndt as old',0,0,0,filkne)
      call iosys ('read character symmetry from kohndt',-1,0,0,sym)
      call iosys ('read integer "no. channels" from kohndt',1,nchan,
     1              0,' ')
      call iosys ('read integer "total no. mos" from kohndt',1,nmotot,
     1             0,' ')
      bsym=.false.
      if (sym.eq.'on') then
	  bsym=.true.
      endif
      maxbf=0
      maxmo=0
      do 10 ch1=1,nchan
         if( posinp('$chan-'//itoc(ch1),cpass) ) then
             call cardin(card)
	 endif
         if (bsym) then
             ngauss(ch1)=intkey(card,'no-aos',1,' ')
	     maxbf=max(maxbf,ngauss(ch1))
	     call intarr(card,'aos',ngch(1,ch1),ngauss(ch1),' ')
	     nmoc(ch1)=intkey(card,'no-mos',1,' ')
	     maxmo=max(maxmo,nmoc(ch1))
	     call intarr(card,'mos',nmoch(1,ch1),nmoc(ch1),' ')
         else
	     ngauss(ch1)=ncon
             maxbf=max(maxbf,ngauss(ch1))
	     nmoc(ch1)=ncon
	     maxmo=max(maxmo,nmoc(ch1))
	     do 20 orb=1,ncon
	        ngch(orb,ch1)=orb
                nmoch(orb,ch1)=orb
   20        continue
	 endif
   10 continue
c----------------------------------------------------------------------c
c                    get memory                                        c
c----------------------------------------------------------------------c
      ntri=nchan*(nchan+1)/2
      grid=1
      basis=grid+4*pntbuf
      ovbeg=basis+ncon*pntbuf
      v=basis+ncon*ncon
      words=v 
      if (modpot) then
          pobeg=v+ntri*pntbuf
          words=pobeg+ntri*ncon*ncon
      endif
      call getscm(words,z,ngot,m6030,0)
      grid=ioff
      basis=grid+4*pntbuf
      ovbeg=basis+ncon*pntbuf
      v=basis+ncon*ncon
      if (modpot) then
          pobeg=v+ntri*pntbuf
      endif
      write (iout,1000) words
c----------------------------------------------------------------------c
c           open potential file to store channel potentials which are  c
c           computed here instead of vstat for this simple model       c
c           problem.                                                   c
c----------------------------------------------------------------------c
      if (modpot) then
          call iosys ('read character "potential filename" from rwf',-1,
     1                 0,0,filpot)
          call iosys ('open vstat as new on ssd',npts*ntri,0,0,
     1                 filpot)
          call modelp(z(v),z(grid),z(grdtyp),type,nchan,nreg,
     1                pntbuf,nolst)
      endif
c----------------------------------------------------------------------c
c                   zero matrix elements                               c
c----------------------------------------------------------------------c
      call rzero(z(ovbeg),ncon*ncon)
      if (modpot) then
	  call rzero(z(pobeg),ntri*ncon*ncon)
      endif
c----------------------------------------------------------------------c
c                   loop over grid                                     c
c----------------------------------------------------------------------c
      npnts=pntbuf
      do 100 ireg=1,nreg
         if (ireg.eq.nreg) then
             npnts=nolst
         endif
         call iosys ('read real '//grdtyp//' from grid without '//
     1               'rewinding',4*npnts,z(grid),0,' ')
c----------------------------------------------------------------------c
c                read in a block of gaussians                          c
c----------------------------------------------------------------------c
         call iosys ('read real "con array" from orbs '//
     1               'without rewinding',npnts*ncon,z(basis),0,' ')
c----------------------------------------------------------------------c
c          do integration over this grid block and accumlate           c
c          we do only the integrals needed for the two channels        c
c          in question but we fill out the whole matix with zeros      c
c----------------------------------------------------------------------c
         call bbovlp(z(ovbeg),z(grid),z(basis),ncon,npnts)
         if (modpot) then
             call iosys ('read real "static potential" from vstat '//
     1                   'without rewinding',npnts*ntri,0,' ')       
             po=pobeg
             do 200 ch1=1,nchan
                do 300 ch2=1,ch1
                   chnpt=ch1*(ch1-1)/2 + ch2
	           call bbpot(z(po),z(grid),z(basis),z(pot),ncon,
     1                        npnts,ntri,ngauss(ch1),ngch(1,ch1),
     2                        ngauss(ch2),ngch(1,ch2),chnpt)
	           po=po+ncon*ncon
  300           continue
  200        continue
  100 continue
      call filtri(z(ovbeg),ncon)
      call iosys ('rewind all on grid read-and-write',0,0,0,' ')
      call iosys ('close grid',0,0,0,' ')
      call iosys ('rewind all on orbs read-and-write',0,0,0,' ')
      call iosys ('close orbs',0,0,0,' ')
      if (modpot) then
          call iosys ('rewind all on vstat read-and-write',0,0,0,' ')
          call iosys ('close vstat',0,0,0,' ')
c----------------------------------------------------------------------c
c          add the kinetic energy integrals to complete the matrix     c
c----------------------------------------------------------------------c
          call iosys ('read real "kinetic integrals" from kohndt',
     1                 ncon*ncon,z(vec),0,' ')
          do 400 ch1=1,nchan
             po=pobeg + ( ch1*(ch1+1)/2 - 1 )*ncon*ncon
             call mkchn(z(po)),z(vec),ncon,ngauss(ch1),ngch(1,ch1))
  400     continue
      endif
c----------------------------------------------------------------------c
c              transform to molecular orbital basis                    c
c----------------------------------------------------------------------c
      if (.not.unitmx) then
	  call iosys ('read character "transformation vector" from '//
     1  	      'kohndt',0,0,0,xform)
          call iosys ('read real '//xform//' from kohndt',ncon*nmotot,
     1                 z(vec),0,' ')
      else
          call onemat(z(vec),ncon)
      endif
      totov=0
      totpo=0
      ov=ovbeg
      po=pobeg
      mov=movbeg
      mpo=mpobeg
      do 400 ch1=1,nchan
	 do 500 ch2=1,ch1
            call tomobs(z(ov),z(mov),z(ov),z(mov),z(vec),z(scri),
     1                  z(scrj),z(scr),z(scr),ncon,nmotot,ngauss(ch1),
     2                  ngch(1,ch1),nmoc(ch1),nmoch(1,ch1),ngauss(ch2),
     3                  ngch(1,ch2),nmoc(ch2),nmoch(1,ch2),'real')
	    ov=ov+ngauss(ch1)*ngauss(ch2)
	    mov=mov+nmoc(ch1)*nmoc(ch2)
	    totov=totov+nmoc(ch1)*nmoc(ch2)
	    if (modpot) then
                call tomobs(z(po),z(mpo),z(po),z(mpo),z(vec),z(scri),
     1                      z(scrj),z(scr),z(scr),ncon,nmotot,
     2                      ngauss(ch1),ngch(1,ch1),nmoc(ch1),
     3                      nmoch(1,ch1),ngauss(ch2),ngch(1,ch2),
     4                      nmoc(ch2),nmoch(1,ch2),'complex')
		po=po+2*ngauss(ch1)*ngauss(ch2)
		mpo=mpo+2*nmoc(ch1)*nmoc(ch2)
		totpo=totpo+2*nmoc(ch1)*nmoc(ch2)
            endif
  500    continue
  400 continue
c----------------------------------------------------------------------c
c                    output matrices                                   c
c----------------------------------------------------------------------c
      call iosys ('write integer "size of bound bound overlap" to '//
     1            'kohndt',1,totov,0,' ')
      call iosys ('write real "bound bound overlap" to kohndt',totov,
     1             z(movbeg),0,' ')
      if (modpot) then
          call iosys ('write integer "size of bound bound potential"'//
     1                ' to kohndt',1,totpo,0,' ')
          call iosys ('write integer "bound bound potential" to kohndt',
     1                1,z(mpobeg),0,' ')
      endif
 1000 format(/,5x,'need',1x,i8,1x,'words for calculation')
      call chainx(0)
      stop
      end
