*deck @(#)mn333.f	1.1  11/30/90
      subroutine mn333(z,a,maxr)
c***begin prologue     m333
c***date written       840717   (yymmdd)
c***revision date      861205   (yymmdd)
c
c   4 december 1986    pws at lanl
c      changing 'namint' and iosys open to character.
c
c   13 june 1986  pws at lanl
c                      modified to handle integral list with duplicate
c                      integrals with non-canonical labels (pws,lanl).
c***keywords           m333, link 333, sort, integrals
c***author             saxe, paul (lanl)
c***source             @(#)mn333.f	1.1   11/30/90
c***purpose            sorts the integral output from m312 into a
c                      triangular matrix oftriangular matrices.
c                      that is, i>j and k>l, but no restriction between
c                      pairs.
c***description
c     m333 currently recognizes the option strings:
c       timing         print timing statistics for this link.
c
c***references
c
c***routines called
c     m333
c       drum(mdutil)
c       iosys(io)
c       traktm(mdutil)
c       getcon(io)
c       sorter(math)
c       getscm(mdutil)
c       unqfil(io)
c       sort32(local)
c       sorter(math)
c       chainx(mdutil)
c
c***end prologue       m333
c
      implicit integer (a-z)
c
      character*128 namint
      character*4096 ops
      character*8 file,prtflg
      logical prnt
      real*8 z(*)
      integer a(*)
c
      common /io/     inp,iout
c
      data prnt/.true./
c
c
c     ----- recover the options string -----
c
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     ----- start timing -----
c
c
c     has printing been turned off externally?
      call iosys('read character "print flag" from rwf',-1,0,0,prtflg)
      if(prtflg.eq.'minimum') prnt=.false.
c
cdir$ fastmd
c
c     ----- dimensions etc. -----
c
      call iosys('read integer "number of basis functions" from rwf',
     $     1,num,0,' ')
      call iosys('read character "integral filename" from rwf',
     $     0,0,0,namint)
      call iosys('open ints as old',0,0,0,namint)
      call iosys('rewind all on ints',0,0,0,' ')
      call iosys('read integer "derivative buffer length" from ints',
     $     1,lenbuf,0,' ')
      nnp=(num+1)*num/2
c
c     ----- divide core -----
c
      lenscr=2*lenbuf
      labels=1
c
c     ----- on 32/64 bit machines, labels is 2 x lenbuf array -----
c
      if (wptoin(1).eq.2) then
         values=iadtwp(labels+2*lenbuf)
      else
         values=iadtwp(labels+lenbuf)
      end if
c
      vals=values+lenbuf
      labs=wpadti(vals+lenscr)
      iork=labs
      jorl=iork+lenbuf
      bins=iadtwp(labs+lenscr)
      ijv=wpadti(bins+lenscr)
      klv=ijv+lenbuf
      rcore=iadtwp(klv+lenbuf)
      acore=wpadti(rcore)
c
c     ----- find out how many derivatives there are -----
c
      call iosys('read integer "number of atoms" from ints',
     $     1,natoms,0,' ')
      nder=natoms*3
c
c     ----- find how much space the sort would like -----
c
      call getscm(0,z,maxpos,'m333 size',0)
      top=min(maxpos,wpadti(rcore+nder*nnp**2))
c
c     ----- get core -----
c
      call getscm(top,z,maxcor,'m333: main',1)
c
c     ----- perform the sort -----
c
      write(iout,11) maxcor
   11 format(1x,'m333:  '/,'     maxcor ',i8)
c
      if(wptoin(1).eq.2) then
       call sort32(a(labels),z(values),a(acore),lenbuf,num,nnp,'ints',
     #           iout,maxcor-acore,z(vals),a(labs),a(bins),
     #           lenscr,a(ijv),a(klv),a(iork),a(jorl),prnt,nder)
      else
       call sort64(a(labels),z(values),a(acore),lenbuf,num,nnp,'ints',
     #           iout,maxcor-acore,z(vals),a(labs),a(bins),
     #           lenscr,a(ijv),a(klv),a(iork),a(jorl),prnt,nder)
      endif
c
c     ----- print timing -----
c
c
c     ----- and exit gracefully ------
c
      return
c
c
      end
