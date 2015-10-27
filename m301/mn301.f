*deck  @(#)mn301.f	1.3 7/30/91
      subroutine mn301(ops)
c***begin prologue     m301
c***date written       840717   (yymmdd)
c***revision date      880318   (yymmdd)
c
c     13 sept 89 bhl at llnl
c               new routines written to drop functions from the
c               integrals list ( see makinx and fix64 )
c
c    28 august   89  bhl at llnl
c       code to create sparate integral files and
c       to exclude integrals from the ordered list
c
c    19 march    88  bhl at llnl
c         maximum maxpos set to 2,000,000 for cos
c
c     1 december 1986  pws at lanl
c         making 'namint' and iosys open character
c
c    13 june 1985   pws  at lanl
c                      modified to handle integral list with duplicate
c                      integrals with non-canonical labels (pws,lanl).
c***keywords           m301, link 301, sort, integrals
c***author             saxe, paul (lanl)
c***source              @(#)mn301.f	1.3 7/30/91
c***purpose            sorts the integral output from l312 into a
c                      triangular matrix oftriangular matrices.
c                      that is, i>j and k>l, but no restriction between
c                      pairs.
c***description
c     m301 currently recognizes the option strings:
c       timing         print timing statistics for this link.
c
c***references
c
c***routines called
c     m301
c       drum(mdutil)
c       iosys(io)
c       traktm(mdutil)
c       getcon(io)
c       sorter(math)
c       getscm(mdutil)
c       unqfil(io)
c       sort1(local)
c       sorter(math)
c       chainx(mdutil)
c
c***end prologue       m301
      implicit integer (a-z)
c
      parameter (maxnbf=2000)
      character*(*) ops
      character*8 prtflg
      character*16 bflabl(maxnbf),nulabl(maxnbf)
      logical prnt
      logical logkey
      real*8 z
      integer a
      pointer(p,z(1)), (p,a(1))
c
      common /io/     inp,iout
c
c
      prnt=logkey(ops,'print=drop',.false.,' ')
c     has printing been turned off externally?
      call iosys('read character "print flag" from rwf',-1,0,0,prtflg)
      if(prtflg.eq.'minimum') prnt=.false.
c
c     ----- dimensions etc. -----
c
      call iosys('read integer "number of basis functions" from rwf',
     $            1,nbasis,0,' ')
c
c     ----- divide core -----
c
c      call manmem(0,p,idum,'mn301',idum)
      pkndx=1
c     these will be passed as the beginning of available common.
      labels=pkndx+nbasis
      zlabls=iadtwp(labels)
c
c     ----- find how much space the sort would like -----
c
c      call iosys('read integer mxcore from rwf',1,mxcore,0,' ')
c      write(iout,*) mxcore, labels
c      if(labels.ge.mxcore) then
c         call lnkerr('mn301: not enough core')
c      endif
c      call manmem(labels,p,ngot,'mn301',0)
      call getmem(labels,p,ngot,'mn301',0)
c
c
c
      if(logkey(ops,'drop',.false.,' ')) then
         call makinx(a(pkndx),nbasis,newnbf,newnnp,bflabl,nulabl,prnt)
         call iosys ('write integer "packing index vector" to rwf',
     $                nbasis,a(pkndx),0,' ')
         call iosys ('write integer "truncated number of basis'
     $             //' functions" to rwf',1,newnbf,0,' ')
         call iosys ('write character "new basis function labels"'
     $             //' to rwf',len(nulabl(1))*newnbf,0,0,nulabl)
c         call manmem(-ngot,p,idum,'mn301',idum)
         call getmem(-ngot,p,idum,'mn301',idum)
      end if
c
      write(iout,*) '      newnbf:',newnbf
c
c     ----- and exit gracefully ------
c
      return
      end
