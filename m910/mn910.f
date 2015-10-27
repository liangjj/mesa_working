*deck @(#)mn910.f	1.4  8/3/91
      subroutine mn910(ops)
c
c***begin prologue     mn910
c***date written       000811   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           ci, hamiltonian matrix
c***author             saxe, paul (lanl)
c***source             @(#)m910.f	1.4   8/3/91
c
c***purpose            to construct the hamiltonian matrix,
c                      or portions thereof.
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       mn910
c
c
      implicit integer (a-z)
c
      character*(*) ops
      character*2 mcorci
      character*32 headr
      character*8 citype
      character*128 gints
      logical logkey
      character*8 prtflg
      character*4 itoc
      character*128 namham
      character*3 answer
      character*16 unit
      integer a
      logical dvd, dagtyp, pfile
      real*8 z
      real*8 rep
      real*8 fzcore
      real*8 cutoff
      real*8 fpkey
      dimension headr(3)
      pointer(p,a(1)),(p,z(1))
c
      common /io/ inp,iout
c
c
      mcorci='ci'
      call iosys('write character mcorci to rwf',0,0,0,mcorci)
      citype='m910'
      call iosys('write character "ci used" to rwf',0,0,0,citype)
c
c     ----- open  the integral file ----
c
      call iosys ('read character "guga integral filename" '//
     $            'from rwf',-1,0,0,gints)
      call iosys('open gints as old',0,0,0,gints)
      unit='gints'
c
      call iosys('read character "print flag" from rwf',-1,0,0,prtflg)
c
c
c     ----- get the constants for dividing up core -----
c
      call iosys('read integer ngroup from rwf',1,ngroup,0,' ')
      call iosys('read integer "symmetries in ci" from rwf',
     $            1,nsym,0,' ')
      call iosys('read integer norbs from rwf',1,norbs,0,' ')
      call iosys('read integer nrows from rwf',1,nrows,0,' ')
      call iosys('read integer nlevs from rwf',1,nlevs,0,' ')
      call iosys('read integer orbfrm from rwf',1,orbfrm,0,' ')
      call iosys('read integer numij from rwf',1,nnp,0,' ')
      call iosys('read integer nmax from rwf',1,nmax,0,' ')
      call iosys('read integer nijvir from rwf',1,nnpvir,0,' ')
      call iosys('read real "nuclear repulsion energy" from rwf',
     $            1,rep,0,' ')
c
c
      pfile=.true.
      call iosys('read integer nwks from rwf',1,nwks,0,' ')
      if (prtflg.ne.'minimum') then
          write (iout,103) nwks
 103      format (5x,'forming the entire hamiltonian of ',i5,
     1               ' configurations')
      end if
c
c        ----- buffer length for storing hamiltonian elements -----
c
      lnbuf=intkey(ops,'ci=buffer-length',20000,' ')
c
      call iosys('write integer "ci buffer length" to rwf',
     $            1,lnbuf,0,' ')
c
      cutoff=fpkey(ops,'ci=cutoff',1.0d-9,' ')
c
c     determine how many roots are desired.
      nroots=intkey(ops,'ci=nroots',1,' ')
      nroots=min(nwks,nroots)
c
c     ----- allocate core for drt arrays, integrals, vectors, etc. -----
c
      imngrp=1
      imxgrp=imngrp+ngroup
      jmngrp=imxgrp+ngroup
      jmxgrp=jmngrp+ngroup
      arc=jmxgrp+ngroup
      wt=arc+4*nrows
      nlwks=wt+4*nrows
      ijadd=nlwks+nrows
      ijgrp=ijadd+nnp
      kadd=ijgrp+nnp
      ladd=kadd+norbs*nsym
      orbsym=ladd+norbs*nsym
      b=orbsym+norbs
      refwt=b+nrows
      refarc=refwt+nlevs
      refb=refarc+nlevs
      irowsv=refb+nlevs
      jrowsv=irowsv+nlevs
      segsv=jrowsv+nlevs
      pagesv=segsv+nlevs
      iwtsv=pagesv+nlevs
      jwtsv=iwtsv+nlevs
      traksv=jwtsv+nlevs
c
      ints=iadtwp(traksv+nlevs)
      acoef=ints+nmax
      bcoef=acoef+nlevs
      c=bcoef+nlevs
      rbuf=c+nwks
      ibuf=wpadti(rbuf+lnbuf)
      need=ibuf+2*lnbuf
c
c     ----- get core for main arrays and then read them in -----
c
      call memory(need,p,ngot,'m910',0)
c
      call iosys('read integer arc from rwf',4*nrows,a(arc),0,' ')
      call iosys('read integer weight from rwf',4*nrows,a(wt),0,' ')
      call iosys('read integer nlwks from rwf',nrows,a(nlwks),0,' ')
      call iosys('read integer ijadd from rwf',nnp,a(ijadd),0,' ')
      call iosys('read integer ijgrp from rwf',nnp,a(ijgrp),0,' ')
      call iosys('read integer kadd from rwf',
     $            norbs*nsym,a(kadd),0,' ')
      call iosys('read integer ladd from rwf',norbs*nsym,
     $            a(ladd),0,' ')
      call iosys('read integer orbsym from rwf',norbs,a(orbsym),0,' ')
      call iosys('read integer b from rwf',nrows,a(b),0,' ')
c
c     ----- fix up the symmetry array and add group mins and maxs -----
c
      call fixdrt(a(orbsym),norbs,a(ijgrp),nnp,a(imngrp),a(imxgrp),
     $            a(jmngrp),a(jmxgrp),ngroup)
c
      call iosys('does "frozen core energy" exist on rwf',0,0,0,
     $            answer)
      if (answer.eq.'no') then
         fzcore=0.0d+00
      else
         call iosys('read real "frozen core energy" from rwf',1,
     #               fzcore,0,' ')
      end if
c
c
c        ----- open either a scratch or permanent file for the 
c        -----               hamiltonian             -----
c
      if(pfile) then
         call iosys('read character "hamiltonian filename" '//
     $              'from rwf',0,0,0,namham)
         call iosys('open hamiltonian as new',0,0,0,namham)
      else
         call iosys('open hamiltonian as scratch',0,0,0,' ')
      end if
      call iosys('write integer "buffer size" to hamiltonian',
     $            1,lnbuf,0,' ')
      prnt=logkey(ops,'print=ci=hmatrix',.false.,' ')
      write(iout,*) 'calculating h'
      headr(1)='buffers'
      headr(2)='"number of elements"'
      headr(3)='diagonals'
      sirow=1
      sjrow=1
      call hamilt(a(arc),a(wt),a(nlwks),a(ijadd),a(kadd),a(ladd),
     #            a(orbsym),a(refwt),a(refarc),a(b),a(refb),
     #            z(ints),z(c),nrows,norbs,nlevs,orbfrm,nsym,
     #            nmax,nwks,nnp,a(irowsv),a(jrowsv),
     #            a(segsv),a(pagesv),a(iwtsv),a(jwtsv),
     #            a(traksv),z(acoef),z(bcoef),lnbuf,a(ibuf),
     #            z(rbuf),cutoff,ngroup,a(imngrp),
     #            a(imxgrp),a(jmngrp),a(jmxgrp),sirow,sjrow,
     #            ops,unit,prtflg,headr)
      call memory(-ngot,p,idum,'m910',idum)
      call diagh('ci',rep,fzcore,lnbuf,nroots,nwks,ops,headr)
      call iosys('close hamiltonian',0,0,0,' ')    
c 
      return
      end
