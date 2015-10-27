*deck @(#)ioexst.f	4.1  7/7/93
      subroutine ioexst(string,pos,ians)
c
c***begin prologue     ioexst
c***date written       850125   (yymmdd)
c***revision date      860112   (yymmdd)
c***keywords           iosys dependent routines
c
c***author             saxe, paul,    (lanl)
c***source             @(#)ioexst.f	4.1   7/7/93
c***purpose            to return 2hno or 3hyes if a file exists or not.
c
c***description        #
c
c
c***references
c
c***routines called    gettok (io)
c                      iounit (io)
c
c   common blocks:     (none)
c
c***end prologue       ioexst
c
      implicit integer (a-z)
c
      character*(*) string,ians
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
c     ----- read the rest of the string. expect:
c       'does <key> (exist) (on) <unit>'
c
      call gettok(key,string,pos)
    1 continue
      call gettok(unit,string,pos)
      if (unit.eq.'exist'.or.unit.eq.'on') go to 1
c
c     ----- find the unit -----
c
      un=iounit(unit,unlist,nunits,error)
      ians='no'
      if (error.ne.0) return
c
c     ----- and find if the file exists -----
c
      do 2 file=unitpt(un)+1,unitpt(un)+nfile(un)
         if (key.eq.keylis(file)) then
            ians='yes'
            return
         end if
    2 continue
c
      return
      end
