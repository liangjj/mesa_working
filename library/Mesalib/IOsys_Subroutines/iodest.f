*deck @(#)iodest.f	5.1  11/6/94
      subroutine iodest(string,pos,i1)
c
c***begin prologue     iodest
c***date written       850125   (yymmdd)
c***revision date      860112   (yymmdd)
c***keywords           iosys dependent routines
c
c***author             saxe, paul,    (lanl)
c***source             @(#)iodest.f	5.1   11/6/94
c***purpose            to destroy an iosys physical unit, or any unit.
c
c***description        #
c
c
c***references
c
c***routines called    gettok (io)
c                      iounit (io)
c                      ioclos (io)
c                      drabsf (cftlib)
c                      lnkerr (mdutil)
c
c   common blocks:     (none)
c
c***end prologue       iodest
c
      implicit integer (a-z)
c
      character*(*) string
c
      character*80 tmplin
      character*32 keylis,key
      character*16 cconst,unlist,unityp,type,unit
      character*8  filtyp
      integer base,readpt,writpt,eof,end,iconst,unitpt,nfile
      integer un,file,nunits,unitno
      real*8 rconst
      logical align,locked
c
      parameter (maxfil=2000,maxun=20)
c
      common /ioqqq1/ keylis(maxfil),filtyp(maxfil),cconst(maxfil),
     #                unlist(maxun),type,key,unit,tmplin,unityp(maxun)
      common /ioqqq2/ base(maxfil),readpt(maxfil),writpt(maxfil),
     #                eof(maxfil),end(maxfil),iconst(maxfil),
     #                rconst(maxfil),unitpt(maxun),nfile(maxun),
     #                align(maxun),locked(maxun),unitno(maxun),
     #                un,file,nunits
c
c
c     ----- find the unit we are dealing with -----
c
      call gettok(unit,string,pos)
c
      un=iounit(unit,unlist,nunits,error)
      if (error.eq.0) then
         unitnm=unitno(un)
         junk=0
         unityp(un)='scratch'
         call ioclos(unit,junk)
      else
cps         call drabsf(0,i1,-1,ierr)
cps         if (ierr.eq.-1) then
cps            call lnkerr('io system: error destroying named unit')
cps         end if
      end if
c
      return
      end
