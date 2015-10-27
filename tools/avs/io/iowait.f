*deck @(#)iowait.f	4.1  7/7/93
      subroutine iowait(string,pos)
c
c***begin prologue     iowait
c***date written       850125   (yymmdd)
c***revision date      860112   (yymmdd)
c***keywords           iosys dependent routines
c
c***author             saxe, paul,    (lanl)
c***source             @(#)iowait.f	4.1   7/7/93
c***purpose            to wait for asynchronous i/o completion on
c                         a named iosys unit.
c***description        #
c
c
c***references
c
c***routines called    gettok (io)
c                      iounit (io)
c                      lnkerr (mdutil)
c                      iowt   (io)
c
c   common blocks:     (none)
c
c***end prologue       iowait
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
c     ----- parse the string for the unit name (or 'all') -----
c
      call gettok(unit,string,pos)
      if (unit.eq.'for') call gettok(unit,string,pos)
      if (unit.eq.'all') then
         unmin=1
         unmax=nunits
      else
         un=iounit(unit,unlist,nunits,error)
         if (error.ne.0) then
            call lnkerr('io system: unable to find unit '//unit//
     #                  ' in iowait')
         end if
         unmin=un
         unmax=un
      end if
c
c     ----- and wait -----
c
      do 4 un=unmin,unmax
         call iowt(unitno(un))
    4 continue
c
      return
      end
