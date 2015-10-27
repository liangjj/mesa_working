*deck @(#)iofile.f	5.1  11/6/94
      subroutine iofile(string,pos,length)
c
c***begin prologue     iofile
c***date written       850125   (yymmdd)
c***revision date      860112   (yymmdd)
c***keywords           iosys dependent routines
c
c***author             saxe, paul,    (lanl)
c***source             @(#)iofile.f	5.1   11/6/94
c***purpose            to create an iosys file.
c
c***description        #
c
c
c***references
c
c***routines called    gettok (io)
c                      iounit (io)
c                      lnkerr (mdutil)
c                      dirsrt (io)
c
c   common blocks:     (none)
c
c***end prologue       iofile
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
c     ----- parse the rest of 'string'. expecting:
c           'create character (file) <key> (on) <unit>'
c                   integer
c                   real
c
      call gettok(type,string,pos)
c
c     ----- check for an appropriate file type -----
c
      if (index('char inte real',type(1:4)).eq.0) then
         call lnkerr('trying to create an illegal type of file...'//
     #               'should be character, integer or real: '//tmplin)
      end if
c
      call gettok(key,string,pos)
      if (key.eq.'file') call gettok(key,string,pos)
      call gettok(unit,string,pos)
      if (unit.eq.'on') call gettok(unit,string,pos)
c
c     ----- find the unit we are to create the file on -----
c
      un=iounit(unit,unlist,nunits,error)
      if (error.ne.0) then
         call lnkerr('io system: trying to create '//key//' on a '
     #               //'non-existant unit '//unit)
      end if
c
c     ----- check if the unit is locked, ie. cannot create a file -----
c
      if (locked(un)) then
         call lnkerr(' io system trying to create file '//key//
     #               ' on unit '//unit//' which is locked')
      end if
c
c     ----- check that the file does not already exist -----
c
      do 4 file=unitpt(un)+1,unitpt(un)+nfile(un)
         if (keylis(file).eq.key) then
            call lnkerr('io system: trying to create a file that '
     #                  //'already exists:'//key//' on '//unit)
         end if
    4 continue
c
c     ----- all ok, so let's create a new file -----
c           checking if the directory has space, of course ....
c
      if (un.lt.nunits.and.unitpt(un)+nfile(un).ge.unitpt(un+1).or.
     #    un.eq.nunits.and.unitpt(un)+nfile(un).ge.maxfil) then
c
c        ----- if got in here, we must move the file lists to
c              create some space to add a new one.
c
         call dirsrt
      end if
c
c     ----- check that the last file is open-space -----
c
      file=unitpt(un)+nfile(un)
      if (keylis(file).ne.'empty space file') then
         call lnkerr('io system: error adding '//key//' to unit '//
     #                unit)
      end if
c
c     ----- set up pointers, and new free-space file -----
c
      keylis(file)=key
      filtyp(file)=type(1:1)
      if (align(un)) base(file)=(base(file)+4096)/4096*4096
      base(file)=(base(file)+itobyt(1)-1)/itobyt(1)*itobyt(1)
      if (length.lt.0) then
         end(file)=-1
         locked(un)=.true.
      else
         end(file)=base(file)+length*iosize(type)
         locked(un)=.false.
         if (align(un)) end(file)=(end(file)+4096)/4096*4096
      end if
c
      file=file+1
      nfile(un)=nfile(un)+1
      keylis(file)='empty space file'
      filtyp(file)='e'
      base(file)=(end(file-1)+itobyt(1)-1)/itobyt(1)*itobyt(1)
      readpt(file)=0
      writpt(file)=0
      eof(file)=base(file)
      end(file)=base(file)
c
      return
      end
