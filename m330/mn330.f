*deck @(#)mn330.f	5.1  11/6/94
      subroutine mn330(killr)
c***begin prologue     mn330.f
c***date written       840717   (yymmdd)
c***revision date      11/6/94
c
c     10 july 1991 bis at lanl
c                  new routine written for dropping functions on
c                  32 bit machines. information needed for
c                  dropping functions picked up from the rwf
c                  where it is written in m310.
c
c     6 oct  89 bhl at llnl
c               putting "old nbf" and "packing vector" on rwf
c               even no packing occurs
c
c     13 sept 89 bhl at llnl
c               new routines written to drop functions from the
c               integrals list ( see makinx and fix64 )
c
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
c***keywords           m330, link 330, sort, integrals
c***author             saxe, paul (lanl)
c***source             @(#)mn330.f	5.1 11/6/94
c***purpose            sorts the integral output from l312 into a
c                      triangular matrix oftriangular matrices.
c                      that is, i>j and k>l, but no restriction between
c                      pairs.
c***description
c     m330 currently recognizes the option strings:
c       timing         print timing statistics for this link.
c
c***references
c
c***routines called
c
c***end prologue       mn330.f
      implicit none
c     --- input variables -----
      integer maxr
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
      integer a
      real*8 z
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer maxnbf
      parameter (maxnbf=2000)
c
      integer inp,iout
      integer nbasis,nnp,lenbuf,nintgr,labels
      integer iork,labs,ijv,bins,scr,lenscr,values,acore
      integer pkindx,newnbf,newnnp,vals,klv,jorl
      integer wpadti,wptoin,iadtwp,top,need,maxcor,maxpos,idum
      character*4096 ops
      character*8 prtflg
      logical prnt
      logical logkey
      logical drop, killr
      character*16 bflabl(maxnbf)
c
      data prnt/.true./
      save prnt
c
      common /io/     inp,iout
      pointer (p,z(1)), (p,a(1))
c
c     --- recover the options string ---
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     --- has printing been turned off externally?
      call iosys('read character "print flag" from rwf',-1,0,0,prtflg)
      if(prtflg.eq.'minimum') prnt=.false.
c
c     --- dimensions etc. ---
      call iosys('read integer "number of basis functions" from rwf',
     $            1,nbasis,0,' ')
      nnp=(nbasis+1)*nbasis/2
c
      call iosys('rewind all on rints',0,0,0,' ')
      call iosys('read integer lenbuf from rints',1,lenbuf,0,' ')
c
c     --- following variable used for 32/64 bit decision ---
      nintgr=wptoin(1)
c
c     --- divide core for sort ---
      lenscr=2*lenbuf
      labels=1
      scr=1
      labs=labels+nintgr*lenbuf
      iork=labs
      jorl=iork+lenbuf
      bins=labs+lenscr
      ijv=bins+lenscr
      klv=ijv+lenbuf
      pkindx=klv+lenbuf
      values=iadtwp(pkindx+nbasis)
      vals=values+lenbuf
      acore=max(scr+nnp,vals+lenscr)
c
c     --- find how much space the sort would like ---
c         an in-core sort requires nnp**2 working precision words
c         beginning at acore
      call iosys('read integer mxcore from rwf',1,maxpos,0,' ')
c     --- try to get just a little more than you need, 
c         just to avoid problems with wpadti, etc.
      top=min(maxpos,wpadti(acore+nnp**2))
      call getmem(top,p,maxcor,'m330: main',0)
c
c     --- number of words left to work with is top-acore
c         read in basis function information needed for drop ---
      drop=logkey(ops,'drop',.false.,' ')
      if(drop) then
         call iosys('read integer "packing index vector" from rwf',
     $               nbasis,a(pkindx),0,' ')
         call iosys('read integer "truncated number of basis'
     $              //' functions" from rwf',1,newnbf,0,' ')
         newnnp=(newnbf+1)*newnbf/2
      endif
c
c     --- perform the sort ---
      if(prnt) write(iout,11)
   11 format(1x,'m330:')
c
c
      if (nintgr.eq.2) then
c
c        --- 32 bit integer machines ---
          if (drop) then
             call fix32(a(labels),z(values),z(acore),lenbuf,newnnp,
     #                  'ints',iout,top-acore,z(vals),
     #                   a(labs),a(bins),lenscr,a(ijv),a(klv),a(iork),
     #                   a(jorl),prnt,ops,a(pkindx),killr)
          else
             call sort32(a(labels),z(values),z(acore),lenbuf,nnp,
     #                   'ints',iout,top-acore,z(vals),
     #                    a(labs),a(bins),lenscr,a(ijv),a(klv),a(iork),
     #                    a(jorl),prnt,ops,killr)
         endif
      else
c
c        --- 64 bit integer machines ---
          if (drop) then
             call fix64(a(labels),z(values),z(acore),lenbuf,newnnp,
     #                  'ints',iout,top-acore,z(vals),
     #                   a(labs),a(bins),lenscr,a(ijv),a(klv),a(iork),
     #                   a(jorl),prnt,ops,a(pkindx),killr)
c
          else
             call sort64(a(labels),z(values),z(acore),lenbuf,nnp,
     #                   'ints',iout,top-acore,z(vals),
     #                    a(labs),a(bins),lenscr,a(ijv),a(klv),a(iork),
     #                    a(jorl),prnt,ops,killr)
          endif
      endif
c
c 
      if (drop) then
c        --- write interesting information to the ints file to allow
c            int=reuse information.
c
c            overwrite the number of basis functions and other info
c            to fool succeeding links.
         call iosys ('write integer "number of basis functions" to rwf',
     $               1,newnbf,0,' ') 
         call iosys ('write integer "number of basis functions"'
     $               //' to ints',1,newnbf,0,' ') 
         call iosys('read character "new basis function labels"'
     $              //' from rwf',len(bflabl(1))*newnbf,0,0,bflabl)
         call iosys('write character "basis function labels" to rwf',
     $              len(bflabl(1))*newnbf,0,0,bflabl)
         call iosys('write character "basis function labels" to ints',
     $              len(bflabl(1))*newnbf,0,0,bflabl)
         call iosys ('write integer "old nbf" to rwf',1,nbasis,0,' ') 
         call iosys ('write integer "old nbf" to ints',1,nbasis,0,' ') 
         call iosys('write integer "truncated number of basis'
     $              //' functions" to ints',1,newnbf,0,' ')
         call iosys('write integer "packing index vector"'
     $              //' to ints',nbasis,a(pkindx),0,' ')
c
      endif
      call getmem(-maxcor,p,idum,'m330: main',idum)
c
c
      return
      end
