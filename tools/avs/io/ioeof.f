*deck @(#)ioeof.f	4.1  7/7/93
      subroutine ioeof(string,pos)
c
c***begin prologue     ioeof
c***date written       850125   (yymmdd)
c***revision date      860112   (yymmdd)
c***keywords           iosys dependent routines
c
c***author             saxe, paul,    (lanl)
c***source             @(#)ioeof.f	4.1   7/7/93
c***purpose            to end an 'endless' file.
c
c***description
c    ioeof is called to end an 'endless' file, ie. one that was created
c    with a length of -1. only one such file can be active at a time on
c    a given unit, and while it is open, no other files can be created.
c    ioeof is used to assign a length to such a file, effectively
c    meaning that it can no longer be extended.
c
c
c***references
c
c***routines called    gettok (io)
c                      iounit (io)
c                      lnkerr (mdutil)
c                      ioflno (io)
c
c   common blocks:     (none)
c
c***end prologue       ioeof
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
c     ----- parse the rest of the string  expecting:
c           'end(file) <key> (on) <unit>'
c
      call gettok(key,string,pos)
      call gettok(unit,string,pos)
      if (unit.eq.'on') call gettok(unit,string,pos)
c
c     ----- find unit and file number -----
c
      un=iounit(unit,unlist,nunits,error)
      if (error.ne.0) then
         call lnkerr('io system: error endfiling endless file '//
     #                key//' on unit '//unit//'. unit does not exist')
      end if
      file=ioflno(key,nfile(un),unitpt(un),keylis,error)
      if (error.ne.0) then
         call lnkerr('io system: error endfiling endless file '//
     #               key//' on unit '//unit//' file does not exist')
      end if
c
c     ----- check that the unit is indeed locked -----
c
      if (.not.locked(un)) then
         call lnkerr('io system: error ending endless file '//key//
     #               ' on unit '//unit)
      end if
      locked(un)=.false.
c
      if (align(un)) then
         end(file)=(eof(file)+4096)/4096*4096-1
      else
         end(file)=eof(file)
      end if
c
c     ----- now update empty space file, which must be next... -----
c
      file=file+1
      if (keylis(file).ne.'empty space file') then
         call lnkerr('io system: error endfiling endless file. problems'
     #               //' with the empty space file')
      end if
c
      base(file)=(end(file-1)+itobyt(1)-1)/itobyt(1)*itobyt(1)
      eof(file)=base(file)
      end(file)=base(file)
c
      return
      end
