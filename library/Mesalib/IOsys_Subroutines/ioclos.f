*deck @(#)ioclos.f	5.1  11/6/94
      subroutine ioclos(string,pos)
c
c***begin prologue     ioclos
c***date written       850125   (yymmdd)
c***revision date      860112   (yymmdd)
c***keywords           iosys dependent routines
c
c***author             saxe, paul,    (lanl)
c***source             @(#)ioclos.f	5.1   11/6/94
c***purpose            to close an iosys unit in an orderly fashion.
c
c***description        #
c
c
c***references
c
c***routines called    gettok (io)
c                      iounit (io)
c                      lnkerr (mdutil)
c                      ioput (cftlib)
c                      iowt   (io)
c                      close  (cftlib)
c
c   common blocks:     (none)
c
c***end prologue       ioclos
c
      implicit integer (a-z)
c
      character*(*) string
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
c
c
c     ----- find the unit we are dealing with -----
c
      call gettok(unit,string,pos)
      if (unit.eq.'all') then
         unmin=1
         unmax=nunits
      else
         un=iounit(unit,unlist,nunits,error)
         if (error.ne.0) then
            call lnkerr('io system: cannot find unit '//unit//
     #                  ' to close it')
         end if
         unmin=un
         unmax=un
      end if
c
c     ----- find the end of information on the unit, and write the
c           directory info out.
c
      do 100 un=unmin,unmax
c
         if (unityp(un).eq.'scratch') go to 90
c
         nfiles=nfile(un)
         file=unitpt(un)+nfiles
         if (keylis(file).ne.'empty space file') then
            unit=unlist(un)
            call lnkerr('io system: cannot find empty-space file '
     #                  //'while closing unit '//unit)
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
c     ----- write the pointers back to the beginning of the file ----
c
         call ioputc(nounit,idum,3*len(idum(1)),0,junk)
c
c     ----- now delete the unit -----
c
   90    continue
c
c        ----- if this is a scratch unit, destroy the file -----
c
         if (unityp(un).eq.'scratch') then
            junk=iocls(unitno(un),'delete')
         else
            junk=iocls(unitno(un),'keep')
         end if
c
c        ----- and delete it from the directory tables -----
c
         if (unit.ne.'all') then
            do 5 i=un+1,nunits
               unlist(i-1)=unlist(i)
               unityp(i-1)=unityp(i)
               unitpt(i-1)=unitpt(i)
               nfile(i-1)=nfile(i)
               align(i-1)=align(i)
               locked(i-1)=locked(i)
               unitno(i-1)=unitno(i)
    5       continue
            nunits=nunits-1
         end if
  100 continue
c
      if (unit.eq.'all') then
         nunits=0
      end if
c
      return
      end
