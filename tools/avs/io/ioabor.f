*deck @(#)ioabor.f	4.1  7/7/93
      subroutine ioabor
c
c***begin prologue     ioabor
c***date written       850125   (yymmdd)
c***revision date      860112   (yymmdd)
c***keywords           io, i/o, closing files, abort, error
c
c***author             saxe, paul,    (lanl)
c***source             @(#)ioabor.f	4.1   7/7/93
c***purpose            to try to save as many iosys files as possible
c                         before aborting a job.
c***description        ioabor attempts to close all the open iosys
c units. if it cannot close any unit, it prints an informative message
c to the output, and continues to the next file. primarily ioabor is
c called from lnkerr (mdutil). this must be a separate routine from the
c iosys routines to avoid recursive calls is lnkerr is called because
c of an iosys error. there are no arguments.
c
c
c***references
c
c***routines called    ioput
c                      iowtab (io)
c
c   common blocks:     io, ioqqq1, ioqqq2
c
c***end prologue       ioabor
c
      implicit integer (a-z)
c
      character*16 idum(3),itoc
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
      common/io/inp,iout
c
c
c
      do 100 un=1,nunits
c
c        ----- if this is a scratch unit, destroy the file -----
c
         if (unityp(un).eq.'scratch') then
            if (iocls(unitno(un),'delete').ne.0) then
               write (iout,2) unlist(un),unitno(un)
    2          format (//' ######## cannot destroy scratch file ',a16,
     #                 ' which is unit:',i4)
            end if
            go to 100
         end if
c
         nfiles=nfile(un)
         file=unitpt(un)+nfiles
         if (keylis(file).ne.'empty space file') then
            write (iout,1) unlist(un),unitno(un)
    1       format (//' ######## cannot close ',a16,' which is unit:'
     #              ,i4)
            go to 100
         end if
c
         pt=base(file)
         nounit=unitno(un)
         idum(1)='iosys'
         idum(2)=itoc(nfiles)
         idum(3)=itoc(pt)
c
         file=unitpt(un)+1
         call ioputc(nounit,keylis(file),len(keylis(1))*nfiles,pt,pt)
         call ioput(nounit,base(file),itobyt(nfiles),pt,pt)
         call ioput(nounit,readpt(file),itobyt(nfiles),pt,pt)
         call ioput(nounit,writpt(file),itobyt(nfiles),pt,pt)
         call ioput(nounit,eof(file),itobyt(nfiles),pt,pt)
         call ioput(nounit,end(file),itobyt(nfiles),pt,pt)
         call ioputc(nounit,filtyp(file),len(filtyp(1))*nfiles,pt,pt)
         call ioputc(nounit,cconst(file),len(cconst(1))*nfiles,pt,pt)
         call ioput(nounit,iconst(file),itobyt(nfiles),pt,pt)
         call ioput(nounit,rconst(file),wptbyt(nfiles),pt,pt)
c
c        ----- write the pointers back to the beginning of the file
c
         call ioputc(nounit,idum,3*len(idum(1)),0,junk)
c
c     ----- now delete the unit -----
c
         junk=iocls(nounit,'keep')
c
  100 continue
c
      nunits=0
c
      return
      end
