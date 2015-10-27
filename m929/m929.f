*deck @(#)pm929.f	1.3  9/3/91
      program m929
      implicit integer(a-z)
      real*8 z,refeng,eroot
      real*8 rep
      dimension split(20)
      integer a
      character*4096 ops
      character*128 nmfile, nmcnfg, ints, namchk
      character*4 itoc,ans
      character*8 scatyp, filtyp, dsk
      character*16 key, chrkey
      character*32 xform
      logical logkey
      logical drop, scat
      pointer(p,z(1)),(p,a(1))
c
      common /io/ inp,iout
      data maxtri/1000000/
c
c
      call drum
      write(iout,99)
  99  format(' m929: direct terms of the scattering hamiltonian ')
c
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     -----  get the drt file name to process the CI information
c
      key=chrkey(ops,'int=drt=key','drt',' ')
      call pakstr(key,lenkey)
      call iosys('read character "checkpoint filename" from rwf',
     $            0,0,0,namchk)
      call iosys('open chk as old',0,0,0,namchk)
      call iosys('read character "drt file name '//key(1:lenkey)
     #            //'" from chk',0,0,0,dsk)
      call iosys('read character "drt unit file name '//dsk
     #               //'" from chk',0,0,0,nmcnfg)
      call iosys('open '//dsk//' as unknown',0,0,0,nmcnfg)
      write(iout,*) '          reading information from '//dsk
c
c     read the number of P and Q space walks from the chk file.
c     we presume that you have run m806 in the mode that enables
c     the user to compute both the target and scattering DRT's in
c     one run.
c
      call iosys('read integer nwksp from chk',1,nwksp,0,' ') 
      call iosys('read integer nwksq from chk',1,nwksq,0,' ') 
      call iosys('close chk',0,0,0,namchk)
      scatyp=chrkey(ops,'scattering','none',' ')
      if(scatyp(1:4).eq.'none') then
         write(iout,1)
         call chainx(0)
         stop
      endif
c
c     -----  read the file name holding the scattering information
c
      if(scatyp(1:4).eq.'kohn') then
         call iosys('read character "kohn filename" from rwf',
     #               -1,0,0,nmfile)
         filtyp='kohn'
         filtyp=filtyp(1:4)
      else if(scatyp(1:8).eq.'r-matrix') then
         call iosys('read character "r-matrix filename" from rwf',
     #               -1,0,0,nmfile)
         filtyp='rmtrx'
         filtyp=filtyp(1:5)
      endif
c
      drop=logkey(ops,'drop',.false.,' ')
      key=chrkey(ops,'m929=drt','drt',' ')
      call pakstr(key,lenkey)
      call iosys('read integer "number of basis functions" from rwf',
     $           1,nbf,0,' ')
      nnp=nbf*(nbf+1)/2
c
c     open integral and scattering files
c
      call iosys('open '//filtyp//' as old',0,0,0,nmfile)
c
      call iosys('read character "integral filename" from rwf',
     $            0,0,0,ints)
      call iosys('open ints as old',0,0,0,ints)
c
c
      call iosys('write integer "number of basis functions" to '
     #            //filtyp,1,nbf,0,' ')
c
c
      call iosys('read real "nuclear repulsion energy" from rwf',
     $     1,rep,0,' ')
      call iosys('write real "nuclear repulsion energy" to '//filtyp,
     $     1,rep,0,' ')
c
      if (drop) then
         call iosys('read integer "old nbf" from rwf',1,oldnbf,0,' ')
         call iosys('write integer "old nbf" to '//filtyp,1,oldnbf,
     #               0,' ')
c
         call getmem(oldnbf,p,ngot,'index',0)
         call iosys('read integer "packing index vector" from rwf',
     $               oldnbf,a,0,' ')
         call iosys('write integer "packing index vector" to '//filtyp,
     $               oldnbf,a,0,' ')
         call getmem(-ngot,p,idum,'index',0)
      endif
c
      need=wpadti(max(nnp,nbf*nbf))
      call getmem(need,p,ngot,'vec',0)
      call iosys('read real "overlap integrals" from rwf',
     $           nnp,z,0,' ')
      call iosys('write real "overlap integrals" to '//filtyp,
     $           nnp,z,0,' ')
c
      call iosys('read character "transformation vector" from rwf',
     $           -1,0,0,xform)
      call iosys('read real '//xform//' from rwf',nbf*nbf,z,0,' ')
      call iosys('write character "transformation vector" to '//filtyp,
     $           0,0,0,xform)
      call iosys('write real '//xform//' to '//filtyp,nbf*nbf,z,0,' ')
      call getmem(-ngot,p,idum,'vec',idum)
c
c     the nwks read here is that of the target states.
c
      call iosys('read integer nwks from '//dsk,1,nwks,0,' ')
      call iosys('read integer "target roots" from '//dsk,
     1            1,nroots,0,' ')
      call iosys('read integer "target orbitals" from '//dsk,1,
     1            nsmall,0,' ')
c
c     from the input compute the total number of spin eigenfunctions
c     in the scattering space, the number of L**2 scattering orbitals
c     and the number of P-space vectors.  remember that each primitive
c     target spin eigenfunction is multiplied by all of the L**2 orbitals
c     to create the P-space of the scattering calculations.  these are
c     then contracted into the physical P-space.
c
      ncsfs=nwksp+nwksq
      nl2=nwksp/nwks
      npvec=nroots*nl2
      nden=nroots*(nroots+1)/2
c
      call iosys('does "direct ao 2e scattering hamiltonian" exist'//
     $           ' on '//filtyp,0,0,0,ans)
      if(ans.eq.'no') then
         call iosys('create real file "direct ao 2e scattering '//
     #               'hamiltonian" on '//filtyp,nden*nnp,0,0,' ')
         call iosys('create real file "direct ao 1e scattering '//
     #              'hamiltonian" on '//filtyp,nnp,0,0,' ')
      end if
c
c----------------------------------------------------------------------c
c         nsmall = no. orbitals used in target ci.                     c
c         ncsfs = total no. of sefs in full calculation.               c
c                 need to run m806 before this step to                 c  
c                 determine this value.                                c 
c         nl2 = no. virtual p-space orbitals                           c
c----------------------------------------------------------------------c 
c
c
      if(nsmall.eq.0) then
         write(iout,*)' nsmall must be input '
         call lnkerr(' m929: input error ')
      end if
      if(ncsfs.eq.0) then
         write(iout,*)' scattering=ncsfs=xxx must be input '
         call lnkerr(' m929: input error ')
      end if
c
c
      call iosys('read real "ci energy 1" from rwf',1,refeng,0,' ')
      call iosys('write real "reference energy" to '//filtyp,
     #           1,refeng,0,' ')
      call iosys('write real "reference energy" to '//dsk,
     #           1,refeng,0,' ')
      write(iout,*)'         reference energy = ',refeng
c
c
      do 7001 i=1,nroots
         call iosys('read real "ci energy '//itoc(i)//'" from rwf',
     #              1,eroot,0,' ')
         call iosys('write real "ci energy '//itoc(i)//'" to '//filtyp,
     #              1,eroot,0,' ')
         call iosys('write real "ci energy '//itoc(i)//'" to '//dsk,
     #              1,eroot,0,' ')
 7001 continue
c
c
      npass=nden
c
      scr=1
      h=scr+nbf*nbf
      v=h+nnp
      tden=v+nnp
      dirct=tden+nden*nnp
      int=dirct+nden*nnp
c
      call iosys('read integer mxcore from rwf',1,maxcor,0,' ')
      left=iadtwp(maxcor)-int
      ntriang=left/nnp
      ntriang=min(nnp,ntriang)
      if(ntriang.lt.1) then
         ntriang=maxtri/nnp
         ntriang=max(1,ntriang)
         need=ntriang*nnp
         left=iadtwp(maxcor)-need
         npass=left/(2*nnp)
         npass=min(npass,nden)
         if(npass.lt.1) then
            write(iout,*)' maxcor left nnp  ',iadtwp(maxcor),left,nnp
            call lnkerr(' m929: insufficient core ')
         end if
      end if
c
      dirct=tden+npass*nnp
      int=dirct+npass*nnp
      need=wpadti(int+ntriang*nnp)
      call getmem(need,p,ngot,'ints',0)
c
c
c     ----- read in  t and v one-electron integrals -----
c
      call iosys('read real "kinetic integrals" from rwf',
     $     nnp,z(h),0,' ')
      call iosys('write real "kinetic integrals" to '//filtyp,
     $     nnp,z(h),0,' ')
      call iosys('read real "potential integrals" from rwf',
     $     nnp,z(v),0,' ')
      call iosys('write real "potential integrals" to '//filtyp,
     $     nnp,z(v),0,' ')
      call vadd(z(h),z(h),z(v),nnp)
c
      call iosys('write real "direct ao 1e scattering hamiltonian" to '
     #            //filtyp,nnp,z(h),0,' ')
c
      call iosys('rewind "direct ao 2e scattering hamiltonian" on '
     #            //filtyp,0,0,0,' ')
c
      ix=0
      ip=0
      do 3 i=1,nroots
         do 2 j=1,i
            ip=ip+1
            call iosys('read real "ao t1pdm:'//itoc(i)//
     $                  itoc(j)//'" from '//filtyp,nbf*nbf,z(scr),0,' ')
            call fold(z(tden+ix),z(scr),nbf)
            ix=ix+nnp
            if(ip.eq.npass) then
               call direct(z(tden),z(int),z(dirct),nnp,npass,ntriang)
               ip=0
               ix=0
               call iosys('write real "direct ao 2e scattering '//
     #                    'hamiltonian" to '//filtyp//
     #                    ' without rewinding',npass*nnp,z(dirct),0,' ')
            end if
   2     continue
   3  continue
c
      if(ip.ne.0) then
         call direct(z(tden),z(int),z(dirct),nnp,ip,ntriang)
         call iosys('write real "direct ao 2e scattering hamiltonian"'
     $               //' to '//filtyp//' without rewinding',
     #               ip*nnp,z(dirct),0,' ')
      end if
      call getmem(-ngot,p,idum,'ints',idum)
c
      npdim=nwksq
      call iosys('write integer "q space walks" to '//dsk,1,
     1            nwksq,0,' ')
      call iosys('write integer npdim to '//filtyp,1,npdim,0,' ')
      call iosys('write integer "target walks" to '//filtyp,1,nwksp,
     #            0,' ')
      write(iout,*)'           target  walks  = ',nwksp
      write(iout,*)'                   npdim  = ',npdim
c
      npvec=nl2*nroots
      civec=1
      pvec=civec+nroots*nwks
      need=wpadti(pvec+nroots*ncsfs*nl2)
      call getmem(need,p,ngot,'ci',0)
c
c     ----- read ci vectors -----
c
      ix=0
      do 100 i=1,nroots
         call iosys('read real "ci root '//itoc(i)//'" from rwf',
     $               nwks,z(civec+ix),0,' ')
         call iosys('write real "ci root '//itoc(i)//'" to '//dsk,
     $               nwks,z(civec+ix),0,' ')
         ix=ix+nwks
  100 continue
c
      ntot=ncsfs*npvec
      call iosys('write integer npvec to '//filtyp,1,npvec,0,' ')
      call iosys('write integer "p space vectors" to '//dsk,1,
     1            npvec,0,' ')
      call iosys('write integer "p space vectors" to '//filtyp,1,
     1            npvec,0,' ')
      call iosys('does "p-space vectors" exist on '//filtyp,
     $            0,0,0,ans)
      if(ans.eq.'no') then
         call iosys('create real file "p-space vectors"'//
     $              ' on '//filtyp,ntot,0,0,' ')
      end if
      call iosys('does "p-space vectors" exist on '//dsk,
     $            0,0,0,ans)
      if(ans.eq.'no') then
         call iosys('create real file "p-space vectors"'//
     $              ' on '//dsk,ntot,0,0,' ')
      end if
c
c
      if(logkey(ops,'scattering=split',.false.,' ')) then
         nsplit=intkey(ops,'scattering=nsplit',0,' ')
         if(nsplit.eq.0 .or. nsplit.gt.20) then
            write(iout,*)'m929:  nsplit = 0 or nsplit gt 20 ',nsplit  
            call lnkerr(' m929: nsplit ')
         end if
         call intarr(ops,'scattering=split',split,nsplit,' ')
         write(iout,*)' expanding ci vector with partitioned'
     $              //' virtual space'
         call svectr(z(civec),z(pvec),nroots,nwks,ncsfs,nl2,
     #               nsplit,split)
      else
         call pvectr(z(civec),z(pvec),nroots,nwks,ncsfs,nl2)
      end if
c
c
      call iosys('write real "p-space vectors" to '//filtyp,
     $           ntot,z(pvec),0,' ')
      call iosys('write real "p-space vectors" to '//dsk,
     $           ntot,z(pvec),0,' ')
      call iosys('rewind all on '//filtyp//' read-and-write',0,0,0,' ')
      call iosys('rewind all on '//dsk//' read-and-write',0,0,0,' ')
c
      if(logkey(ops,'print=scattering=p-space',.false.,' ')) then
      call prvec(z(pvec),nroots,nl2,ncsfs)
      end if
      call getmem(-ngot,p,idum,'ci',idum)
c
c
 1    format(/,1x,'This link is only for scattering calculations. Quit')
      call chainx(0)
c
c
      stop
      end



