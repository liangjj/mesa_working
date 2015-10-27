*deck  @(#)pm930.f	5.1 11/6/94
      program m930
      implicit integer(a-z)
c
      real*8 r
      integer a
      integer ihdr(10)
      logical logkey, scat
      character*128 namham,nmfile,nmcnfg,namchk,junk,keystr
      character*4096 ops
      character*8 scatyp, filtyp, dsk
      character*16 chrkey, key
      character*3 ans
      common /io/inp,iout
      pointer (p,r(1)), (p,a(1))
c

c      namelist /iput/mdim,lsize,iprt,itest,ivec,npvec,memory,memsrt,
c     $               energy,stopj,tstopt,npdim,npmax,iout,memopt
c
      data itest/0/,lsize/128/,iprt/0/,ivec/0/,npvec/1/
      save itest,lsize,iprt,ivec,npvec
      call drum
      call iosys('read character options from rwf',-1,0,0,ops)
      call iosys('read character "checkpoint filename" from rwf',
     $            0,0,0,namchk)
      call iosys('open chk as old',0,0,0,namchk)
      key=chrkey(ops,'int=drt=key','drt',' ')
      call pakstr(key,lenkey)
      call iosys('read character "drt file name '//key(1:lenkey)
     #            //'" from chk',0,0,0,dsk)
      call iosys('read character "drt unit file name '//dsk
     #               //'" from chk',0,0,0,nmcnfg)
      call iosys('open '//dsk//' as unknown',0,0,0,nmcnfg)
      write(iout,*) '          reading information from '//dsk
      call iosys('close chk',0,0,0,namchk)
      scatyp=chrkey(ops,'scattering','none',' ')
      if(scatyp(1:4).eq.'none') then
         write(iout,1)
         call chainx(0)
         stop
      endif
      if(scatyp(1:4).eq.'kohn') then
         call iosys('read character "kohn filename" from rwf',
     #               -1,0,0,nmfile)
         filtyp='kohn'
         filtyp=filtyp(1:4)
         len=4
      else if(scatyp(1:8).eq.'r-matrix') then
         call iosys('read character "r-matrix filename" from rwf',
     #               -1,0,0,nmfile)
         filtyp='rmtrx'
         filtyp=filtyp(1:5)
         len=8
      endif    
      call iosys('read character "hamiltonian filename" from rwf',
     $           0,0,0,namham)
      call iosys('open hamiltonian as old',0,0,0,namham)
      call iosys('open '//filtyp//' as old',0,0,0,nmfile)
c
c
c     convert maxcor to a working precision value
      call getmem(0,p,ngot,'m930',0)
      call iosys('read integer mxcore from rwf',1,maxcor,0,' ')
      rcore=iadtwp(maxcor)
      memory=65536
c
c
      write(iout,101)
  101 format(' m930: hamiltonian sort for the optical potential')
c
      call iosys('read integer nwksp from '//dsk,1,nwksp,0,' ')
      call iosys('read integer nwksq from '//dsk,1,nwksq,0,' ')    
      mdim=nwksp+nwksq
      npdim=nwksq
c
      if(logkey(ops,'scattering='//scatyp(1:len)//'=test-hamiltonian',
     #               .true.,' ')) then
         itest=1
      endif
      if(logkey(ops,'scattering='//scatyp(1:len)//'=print-hamiltonian',
     #                 .false.,' ')) then
         iprt=1
      endif
      if(logkey(ops,'scattering='//scatyp(1:len)//'=test-hopt',
     #               .false.,' ')) then
         ivec=1
      endif
      lsize=intkey(ops,'scattering='//scatyp(1:len)//'=lsize',lsize,' ')
      call iosys('read integer "p space vectors" from '//filtyp,1,
     1            npvec,0,' ')
      npvec=intkey(ops,'scattering='//scatyp(1:len)//'=npvec',npvec,' ')
c
c

c
      if(logkey(ops,'scattering='//scatyp(1:len)//'=hqq',
     #               .false.,' ')) then
         formqq=1
         call iosys('read integer npdim from '//filtyp,1,npdim,0,' ')
         write(iout,*) ' npdim = ', npdim
         twalks=mdim-npdim
         if(twalks.le.0) then
            call lnkerr(' twalks le 0 ')
         endif
         maxnpdim=intkey(ops,'scattering='//scatyp(1:len)//'=maxnpdim',
     #                        400,' ')
         if(npdim.gt.maxnpdim) then
            write(iout,*)' the max. npdim (dim of hqq) is ',maxnpdim
            write(iout,*)' npdim is                       ',npdim
            call lnkerr(' npdim gt max. npdim ')
         end if
      else
         formqq=0
         npdim=0
         twalks=0
      end if
c
c
      ncore=rcore-2048
      iget=ncore+1024
      need=wptoin(iget)
      call getmem(need,p,ngot,'m930',0)
c
      if(formqq.ne.0) then
         ncore=ncore-npdim*npdim-1
         if(ncore.lt.memory) then
            write(iout,*)' hqq larger than memory'
            call lnkerr(' m930: memory problem ')
         end if
      end if
c
      call iosys('read integer "ci buffer length" from rwf',
     $            1,lnbuf,0,' ')
c
      lxbuf=max(lnbuf,mdim)
      lnblk=lsize*lsize
      memory=max(lnblk,memory)
      nblks=memory/lnblk
      blksiz=nblks*lnblk
      mxblk=(mdim-1)/lsize+1
      tblks=(mxblk+1)*mxblk/2
      tsize=tblks*lnblk
      nbins=(tblks-1)/nblks+1
c
      ncore=ncore-lxbuf-iadtwp(4*nbins-2*lnbuf)
      binsiz=ncore/nbins
      binsiz=(binsiz-2)/2
      binsiz=min(binsiz,blksiz)
c
c
c -- open the new hamiltonian file
c
      nbpb=(blksiz-1)/binsiz+1
      lnfile=4*nbpb*nbins*(2*binsiz+2)+nbins+200000
      lnham=lnfile-150000
c
c     ----- allocate core -----
      hqq=1
      hbuf=hqq+npdim*npdim
      ibuf=wpadti(hbuf+lxbuf)
      bptr=ibuf+2*lnbuf
      bcntr=bptr+nbins
      xbin=iadtwp(bcntr+nbins)
      ibin=wpadti(xbin+binsiz*nbins)
      need=iadtwp(ibin+(binsiz+2)*nbins)
c
      if(need.gt.iget) then
         write(iout,*)' need ncore ',need,iget
         call lnkerr(' core allocation problem ')
      end if
c
      write(iout,98)
      write(iout,99) lsize,lnblk,nblks,blksiz,mxblk,tblks,tsize,
     1               nbins,binsiz,lnfile,nbpb,need,formqq,npdim,twalks
c
  98  format(5x,' sort parameters ')
  99  format(5x,' lsize lnblk nblks  blksiz mxblk  ',5i8,/,5x,
     1         ' tblks tsize nbins  binsiz lnfile ',5i8,/,5x,
     2         ' nbpb  need  formqq npdim  twalks ',5i8)
c
c
      ihdr(1)=mdim
      ihdr(2)=lsize
      ihdr(3)=lnblk
      ihdr(4)=nbins
      ihdr(5)=binsiz
      ihdr(6)=nblks
      ihdr(7)=blksiz
      ihdr(8)=tsize
      ihdr(9)=mxblk
c
      call iosys('write integer header to '//filtyp,10,ihdr,0,' ')
c
      call iosys('does "sorted hamiltonian" exist on '//filtyp,0,0,0,
     #            ans)
      icreat=0
      if(ans.eq.'no') then
         icreat=1
         call iosys('create integer file "sorted hamiltonian" on '
     #               //filtyp,-1,0,0,' ')
      endif
c
c
c
      binsz2=binsiz+2
      start=nbins
c
      if(formqq.eq.1) then
         call hqqsort(a(bptr),a(bcntr),r(xbin),a(ibin),r(hbuf),
     $                a(ibuf),nbins,binsiz,lnbuf,lnblk,blksiz,
     $                lsize,start,iu,binsz2,mdim,r(hqq),npdim,twalks,
     #                filtyp)
      else
         call bsort(a(bptr),a(bcntr),r(xbin),a(ibin),r(hbuf),
     $              a(ibuf),nbins,binsiz,lnbuf,lnblk,blksiz,
     $              lsize,start,iu,binsz2,mdim,filtyp)
      end if
c
      if(icreat.ne.0) then
         call iosys('endfile "sorted hamiltonian" on '//filtyp,0,0,0,
     #              ' ')
      end if
      if(formqq.eq.1) then
         call iosys('write real hqq to '//filtyp,npdim*npdim,r(hqq),
     #               0,' ')
      end if
c
c
      if(ivec.ne.0) then
         write(iout,202) npvec, filtyp
 202     format(//,' npvec = ',i6,'  unit vectors written to '
     #             'file = ',a8,
     #          /,' used for optical potential test ')
         ncsfs=intkey(ops,'scattering=ncsfs',mdim,' ')
         vec=xbin
         need=vec+npvec*ncsfs
         if(need.gt.ncore) then
            write(iout,*)' pvcout: need ncore ',need,ncore
            call lnkerr(' insufficient memory for pvcout ')
         end if
c
         call iosys('does "test vectors" exist on '//filtyp,0,0,0,ans)
         if(ans.eq.'no') then
            call iosys('create real file "test vectors" on '//filtyp,
     $                  npvec*ncsfs,0,0,' ')
         end if
         filsiz=npvec*ncsfs
         write(iout,*) 'filsiz = ',filsiz
c
c  
         call pvcout(r(vec),ncsfs,npvec,11,filtyp)
      end if
      if(itest.gt.0) then
         ntvec=1
c        ----- allocate some more core -----
         xb=xbin
         ib=wpadti(xb+binsiz)
         xm=iadtwp(ib+binsiz+2)
         vec=xm+blksiz
         svec=vec+mdim*ntvec
         tvec=svec+mdim*ntvec
         need=tvec+lsize*ntvec
         if(need.gt.iget) then
            write(iout,*)' test: need ncore ',need,iget
            call lnkerr(' insufficient memory for test ')
         end if
c
c
         call iosys('read real "ci root 1" from rwf',mdim,r(vec),0,' ')
         call hamtest(a(bptr),nbins,r(xb),a(ib),r(xm),binsiz,binsz2,
     1           tblks,nblks,lsize,blksiz,mxblk,iu,iprt,
     2           r(vec),r(svec),r(tvec),mdim,ntvec,filtyp)
      end if
c
c
      call getmem(-ngot,p,idum,'m930',idum)
      call chainx(0)
c
c
      stop
 1    format(/,1x,'This link is used only for scattering '
     #            'calculations. Exit')
      end
