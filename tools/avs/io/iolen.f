*deck @(#)iolen.f	4.1  7/7/93
      subroutine iolen(string,pos,len)
c
c***begin prologue     iolen
c***date written       850125   (yymmdd)
c***revision date      860112   (yymmdd)
c***keywords           iosys dependent routines
c
c***author             saxe, paul,    (lanl)
c***source             @(#)iolen.f	4.1   7/7/93
c***purpose            to find the length, space available or
c                         read-write pointers of a file.
c***description        #
c
c
c***references         #
c
c***routines called    gettok (io)
c                      ffnext (chr)
c                      lnkerr (util)
c
c   common blocks:     (none)
c
c***end prologue       iolen
c
      implicit integer (a-z)
c
      character*(*) string
      character*16 junk,ffnext
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
c     ----- parse the string...looking for:
c           'length (of) <key> (on) <unit>'
c    or 'get [length,read,write (pointer)] of <key> on <unit>'
c              maximum length
    1 continue
         call gettok(key,string,pos)
         pos1=0
         junk=ffnext(key,pos1,s,e)
      if (index('maximum length read write pointer of',key(s:e)).ne.0)
     #             go to 1
    2 continue
         call gettok(unit,string,pos)
         pos1=0
         junk=ffnext(unit,pos1,s,e)
      if (index('on',unit(s:e)).ne.0) go to 2
c
c     ----- find local unit number -----
c
      un=iounit(unit,unlist,nunits,error)
      if (error.ne.0) then
         len=-1
         return
      end if
c
c     ----- find file and length -----
c
      file=ioflno(key,nfile(un),unitpt(un),keylis,error)
c
      if (error.ne.0) then
         len=-1
         return
      end if
c
c        ----- arrays or files kept on disc -----
c
      size=iosize(filtyp(file))
      if (index(string,'maximum length').ne.0) then
         len=(end(file)-base(file))/size
      else if (index(string,'length').ne.0) then
         len=max(0,eof(file)-base(file))/size
      else if (index(string,'write').ne.0) then
         len=writpt(file)/size
      else if (index(string,'read').ne.0) then
         len=readpt(file)/size
      else
         call lnkerr('io system: do not understand what to get')
      end if
c
      return
      end
