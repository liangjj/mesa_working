*deck @(#)iodump.f	5.1  11/6/94
      subroutine iodump
c
c***begin prologue     iodump
c***date written       850125   (yymmdd)
c***revision date      860112   (yymmdd)
c***keywords           iosys dependent routines
c
c***author             saxe, paul,    (lanl)
c***source             @(#)iodump.f	5.1   11/6/94
c***purpose            to list the iosys units and files open.
c
c***description        #
c
c
c***references
c
c***routines called    (none)
c
c   common blocks:     io
c
c***end prologue       iodump
c
      implicit integer (a-z)
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
      common /io/     inp,iout
c
      nfiles=0
      do 1 i=1,nunits
         nfiles=nfiles+nfile(i)
    1 continue
      write (iout,2) maxun,maxfil,nunits,nfiles
    2 format (//,' ------------- i/o system table dump -------------',/,
     #           '   maximum number of units:',i5,/,
     #           '   maximum number of files:',i5,/,
     #           '    actual number of units:',i5,/,
     #           '    actual number of files:',i5)
c
c     ----- dump files for each unit -----
c
      do 10 un=1,nunits
         write (iout,3) un,unlist(un),unitno(un),nfile(un)
    3    format (//,'   dump of unit:',i3,' named: ',a16,
     #           ' physical name: ',a6,/,'   number of files:',i4,//,
     #           '     name      base   readpt writpt eof    end',/)
         do 5 file=unitpt(un)+1,unitpt(un)+nfile(un)
            write (iout,4) keylis(file),base(file),
     #                    readpt(file),writpt(file),eof(file),end(file)
    4       format (1x,a32,1x,i7,i7,i7,i7,i7)
    5    continue
   10 continue
c
      return
      end
