*deck @(#)m935.f	5.1  11/6/94
      program m935
c***begin prologue     m935
c***date written       890701  yymmdd  
c***revision date      910725  yymmdd
c
c  13 july   1989      bhl at llnl
c     adding a loop over scattering energies
c  15 july   1989      bhl at llnl
c     adding an option to directly invert hqq
c     also an option to iteratively invert hqq
c  25  july   1991     rlm at lanl
c     working on 32-bit version. primarily adding appropriate
c     iadtwp functions. changing the 'sorted hamiltonian' file to
c     be type integer.
c
c***keywords           
c***author             lengsfield, byron 
c***source             @(#)m935.f	5.1   11/6/94
c***purpose            to form the optical potential          
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       m935.f
      implicit integer(a-z)
      parameter (maxeng=200)
      real*8 z
      integer a
      real*8 stopj,fpkey,refeng,escatt
      real*8 energy(maxeng),esave(maxeng)
      integer ihdr(10)
      character*4096 ops
      character*128 nmfile
      character*8 scatyp, filtyp, chrkey
      character*3 ans
      logical logkey
c
c
      data npmax/10/,npdim/0/,iprt/0/
      data lsize/128/
      data energy/maxeng*-9999.d00/, stopj/0.0001d0/
      save npmax,npdim,iprt,lsize,energy,stopj
c
      common /io/inp,iout
      pointer(p,z(1)), (p,a(1))
c
c
      call drum
      write(iout,*) ' m935: Optical Potential'
      call getmem(0,p,ngot,'m935',0)
      call iosys('read integer mxcore from rwf',1,maxcor,0,' ')
c     convert to working precision core available
      call getmem(maxcor,p,ngot,'m935',0)
      ncore=iadtwp(maxcor)
c
      call iosys('read character options from rwf',-1,0,0,ops)
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
      call iosys('open '//filtyp//' as old',0,0,0,nmfile)
c
      call iosys('read integer header from '//filtyp,10,ihdr,0,' ')
c
c     ----- look for options in force -----
      if(logkey(ops,'scattering='//scatyp(1:len)//'=freeze',
     #               .false.,' ')) then
         freeze=1
      else
         freeze=0
      end if
      if(logkey(ops,'scattering='//scatyp(1:len)//'=hqq',
     #               .false.,' ')) then
         hqq=1
      else
         hqq=0
      end if
c
c     ----- static exchange -----
      if (logkey(ops,'scattering='//scatyp(1:len)//'=no-optical-'//
     #               'potential',.false.,' ')) then
          noopt=1
          write (iout,*) ' static exchange run '
          write (iout,*) ' no optical potential constructed '
      else
          noopt=0
      endif      
c
      stopj=fpkey(ops,'scattering='//scatyp(1:len)//'=stopj',stopj,' ')
c
      tstopt=0
      if(logkey(ops,'scattering='//scatyp(1:len)//'=test-hopt',
     #                             .false.,' ')) then
         tstopt=1
      endif
c
      call iosys('read integer npvec from '//filtyp,1,npvec,0,' ')
      call iosys('read integer npdim from '//filtyp,1,npdim,0,' ')
      call iosys('read integer "target walks" from '//filtyp,
     #            1,twalks,0,' ')
c
      npmax=intkey(ops,'scattering='//scatyp(1:len)//'=npmax',npmax,' ')
      npdim=intkey(ops,'scattering='//scatyp(1:len)//'=npdim',npdim,' ')
      npvec=intkey(ops,'scattering='//scatyp(1:len)//'=npvec',npvec,' ')
      lsize=intkey(ops,'scattering='//scatyp(1:len)//'=lsize',lsize,' ')
      mxmax=intkey(ops,'scattering='//scatyp(1:len)//'=mxiter',
     #                  mxmax,' ')
      if (mxmax.gt.200) then
          mxmax=200
          write (iout,*) ' mxiter set to 200 '
      endif
c
      call fparr(ops,'scattering='//scatyp(1:len)//'=energy',energy,
     #                maxeng,' ')
      call iosys('read real "reference energy" from '//filtyp,
     #            1,refeng,0,' ')
      numpts=0
      do 10 i=1,maxeng
         if(energy(i).ne.-9999.d0) then
            numpts=numpts+1
            esave(i)=energy(i)
            energy(i)=energy(i)+refeng
         else
            go to 11
         end if
  10  continue
  11  continue
c
      call iosys('read integer nwks from rwf',1,mdim,0,' ')
      npmax=min(npmax,npvec)
      if (hqq.eq.0) then
          if (freeze.ne.0) then
              msize=mdim-twalks
          else
              msize=mdim
          endif
          if (10*npmax.gt.msize) then
              npmax=msize/10
              npmax=max(npmax,1)
          endif
      endif
      if(numpts.eq.0) then
         write(iout,*)' no scattering energies read from input '
         call lnkerr(' numpts = 0 in m935 ')
      end if
c
      call iosys('write integer nenergy to '//filtyp,1,numpts,
     #            0,' ')
      call iosys('write real "scattering energies" to '//filtyp,
     #            numpts,esave,0,' ')
c
      write(iout,99001) npmax,npvec,lsize,numpts,freeze,hqq,
     #                  stopj,refeng,(esave(k),k=1,numpts)
      write(iout,99002)(energy(k),k=1,numpts)
99001 format('      optical potential input parameters ',/,
     #       '        npmax   = ',i8,/,
     #       '        npvec   = ',i8,/,
     #       '        lsize   = ',i8,/,
     #       '        numpts  = ',i8,/,
     #       '        freeze  = ',i8,/,
     #       '        hqq     = ',i8,/,
     #       '        stopj   = ',f14.8,/,
     #       '        reference  energy    ',f14.8,/,
     #       '        rel. scattering energies (au) ',
     #                40(/,5(8x,f14.8)))
99002 format('      abs. scattering energies (au) ',
     #              40(/,5(8x,f14.8)))
c
c     check to see if hopt need be created on the file scattering
c
      call iosys('does hopt exist on '//filtyp,0,0,0,ans)
      if(ans.eq.'no') then
         call iosys('create real hopt on '//filtyp,numpts*npvec*npvec,
     #               0,0,' ')
         call iosys('create real hpp on '//filtyp,numpts*npvec*npvec,
     #               0,0,' ')
      else
         call iosys('get maximum length of hpp on '//filtyp,
     #               lenhpp,0,0,' ')
         tstpts=lenhpp/(npvec*npvec)
         if(tstpts.eq.0) then
            call lnkerr(' cant overwrite hpp on '//filtyp)
         endif
         if(tstpts.lt.numpts) then
            write(iout,*) ' ** warning **'
            write(iout,*)'   the file hpp already exists on scattering '
            write(iout,*)'   the number of scattering energies is',
     #                       ' being reduced to ',tstpts
            numpts=tstpts
            write(iout,*) ' ** warning **'
        end if
      end if
      nbins=ihdr(4)
      binsiz=ihdr(5)
      nblks=ihdr(6)
      blksiz=ihdr(7)
      tsize=ihdr(8)
      mxblk=ihdr(9)
      minsiz=blksiz+2*binsiz+2
c
c
      npass=(npvec-1)/npmax+1
c
      if(abs(stopj).lt.1.d-9) stopj=1.d-04
      if(hqq.eq.0) then
         mxiter=15*npmax
         mxiter=max(mxiter,mxmax)
         nrhs=npmax
         bblock=(mxiter+2*nrhs)*mdim
         tblock=mdim*npvec 
      else
         nrhs=npvec
         mxiter=1
         bblock=1
         tblock=1
      end if
c
c     ----- allocate core -----
c
      bptr=1
      hopt=iadtwp(bptr+nbins)
      hpp=hopt+npvec*npvec
      b=hpp+npvec*npvec
      t=b+bblock
      c=t+tblock
      tm=c+mxiter*nrhs
      scrtch=tm+mxiter*(mxiter+nrhs)
      c0=scrtch+max(2*mxiter,npvec*npvec,2*npdim)
      hess=c0+mdim*npvec
      if(hqq.eq.0) then
         thess=hess+mxiter*(mxiter+1)/2
         ogl=thess+mxiter*(mxiter+1)/2
         gl=ogl+mxiter*nrhs+nrhs
         tgl=gl+mxiter*nrhs+nrhs
         td=tgl+mxiter*nrhs+nrhs
         grad=td+mdim
      else
         grad=hess+npdim*(npdim+npvec)
      end if
      diag=grad+mdim*npvec
      rmsd=diag+mdim
      r=rmsd+nrhs
      tvec=r+nrhs*mxiter
      hss=tvec+lsize*npvec
c
      left=ncore-hss+1
      if(left.gt.(tsize+iadtwp(2*binsiz+2))) then
         incore=1
         xbin=hss+tsize
         ibin=wpadti(xbin+binsiz)
         need=iadtwp(ibin+binsiz+2)
c
         ihave=need
         write(iout,*)'     incore inversion'
      else
         incore=0
         if(left.lt.minsiz) then
            need=hss+minsiz
            write(iout,*)' left need ',left,need
            write(iout,*)' largest memory block is ',
     $                   ' big = (mxiter+2*nrhs)*mdim '
            ibig=(mxiter+2*nrhs)*mdim
            write(iout,7891) mxiter,nrhs,mdim,ibig
 7891       format(/,' mxiter = ',i8,/,
     $              ' nrhs   = ',i8,/,
     $              ' mdim   = ',i8,/,
     $              ' big    = ',i8)
            call lnkerr(' insufficient core to proceed ')
         else
            xbin=hss+blksiz
            ibin=wpadti(xbin+binsiz)
            need=iadtwp(ibin+binsiz+2)
         endif
c
         ihave=need
         write(iout,*)'      out of core inversion '
      endif
c
c     read the pointers
      call getham(idah,a(bptr),nbins,0,filtyp)
c
c     ----- loop over scattering energies -----
c
      call iosys('rewind hopt on '//filtyp,0,0,0,' ')
      call iosys('rewind hpp on '//filtyp,0,0,0,' ')
      do 60 js=1,numpts
         escatt=energy(js)
         call sethpp(z(hpp),z(hss),z(grad),z(diag),z(t),mdim,npvec,
     $               z(c0),npvec,npdim,z(scrtch),incore,lnbuf,escatt,
     $               nbins,binsiz,nblks,blksiz,tsize,mxblk,lsize,
     $               a(bptr),z(xbin),a(ibin),z(tvec),idah,idav,
     $               iprt,tstopt,freeze,hqq,filtyp)
c
c        ----- static exchange if test -----
         if (noopt.eq.0) then
            if(hqq.eq.0) then
               j=0
               do 50 i=1,npvec,npmax
                  nvec=min(npmax,npvec-i+1)
                  call solver(z(hss),z(b),z(c),z(t+j),z(tm),
     $                        z(scrtch),z(c0),z(hess),z(thess),
     $                        z(ogl),z(gl),z(tgl),z(grad+j),z(td),
     $                        z(diag),z(rmsd),z(r),mdim,nvec,nrow,
     $                        incore,npvec,mxiter,stopj,lnbuf,
     $                        escatt,nbins,binsiz,nblks,blksiz,
     $                        tsize,mxblk,lsize,a(bptr),z(xbin),
     $                        a(ibin),z(tvec),idah,freeze,npdim)
                  j=j+nvec*mdim
  50           continue
c
               call ebtc(z(hopt),z(t),z(grad),npvec,mdim,npvec)
            else
               call invhqq(z(hess),z(c0),z(grad),
     $                     npvec,npdim,mdim,escatt,filtyp)
               call ebtc(z(hopt),z(c0),z(grad),npvec,npdim,npvec)
            end if
            call iosys('write real hopt to '//filtyp//
     #                                     ' without rewinding',
     #                  npvec*npvec,z(hopt),0,' ')
            call iosys('write real hpp to '//filtyp//
     #                                    ' without rewinding',
     #                  npvec*npvec,z(hpp),0,' ')
         else
c           ..static exchange code
c
            call iosys('write real hopt to '//filtyp//
     #                                     ' without rewinding',
     #                  npvec*npvec,z(hpp),0,' ')
            call iosys('write real hpp to '//filtyp//
     #                                    ' without rewinding',
     #                  npvec*npvec,z(hpp),0,' ')
c
         end if
c
c
         if(logkey(ops,'scattering=print=hopt',.false.,' ')) then
            write(iout,101) escatt
 101        format(//,'  scattering energy ',f20.12)
            write(iout,*)' optical potential '
            call matout(z(hopt),npvec,npvec,npvec,npvec,iout)
            write(iout,*)' '
            write(iout,*)' hpp - e '
            call matout(z(hpp),npvec,npvec,npvec,npvec,iout)
         end if
c
         if(tstopt.ne.0) then
            call testh(z(hpp),z(hopt),z(c0),z(t),npvec,iu,
     $                  mdim,iout)
         end if
  60  continue
c
c     ----- end loop over energies -----
c           exit gracefully
      call getmem(-ngot,p,idum,'m935',idum)
      call chainx(0)
      stop
 1    format(/,1x,'This link only run for a scattering calculation.'
     #            ' Exit')
      end
