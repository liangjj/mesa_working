*deck @(#)iocopy.f	4.1  7/7/93
      subroutine iocopy(string,pos,buffer,lenbuf,char)
c
c***begin prologue     iocopy
c***date written       860123   (yymmdd)
c***revision date      870706   (yymmdd)
c     6 july 1987      pws at lanl
c         modifying to use ioget, ioput and byte sizes for i/o
c
c***keywords           iosys dependent routines
c
c***author             saxe, paul,    (lanl)
c***source             @(#)iocopy.f	4.1   7/7/93
c***purpose            to copy an iosys file from one unit to another,
c                      possibly changing the name in the process.
c
c***description        #
c
c
c***references
c
c***routines called    gettok (io)
c                      iounit (io)
c                      lnkerr (mdutil)
c                      ioflno (io)
c                      ioget (cftlib)
c                      unitwt (io)
c
c   common blocks:     (none)
c
c***end prologue       iocopy
c
      implicit integer (a-z)
c
      character*(*) string,char
      integer buffer(lenbuf)
      character*32 key1
      character*16 unit1
      character*32 key2
      character*16 unit2
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
c     ----- parse the rest of 'string', expecting:
c           'copy <key> (from) <unit> (to) (<key>) (on) <unit>'
c
      call gettok(key1,string,pos)
      if (index('character integer   real        ',key1(1:9))
     #     .ne.0) call gettok(key1,string,pos)
c
      call gettok(unit1,string,pos)
      if (unit1.eq.'from') call gettok(unit1,string,pos)
c
      call gettok(key2,string,pos)
      if (key2.eq.'to') call gettok(key2,string,pos)
c
      call gettok(unit2,string,pos)
      if (unit2.eq.'on') call gettok(unit2,string,pos)
c
c     ----- if there is no key specified, we should get eos from
c           gettok, and make key name the same.
c
      if (unit2.eq.'eos') then
         unit2=key2
         key2=key1
      end if
c
c     ----- find the unit from its name -----
c
      un1=iounit(unit1,unlist,nunits,error)
      if (error.ne.0) then
         if (char.eq.'noerror') return
         call lnkerr('io system: error trying to copy '//key1//' from'
     #               //' non-existant unit: '//unit1)
      end if
c
      un2=iounit(unit2,unlist,nunits,error)
      if (error.ne.0) then
         if (char.eq.'noerror') return
         call lnkerr('io system: error trying to copy '//key2//' to'
     #               //' non-existant unit: '//unit2)
      end if
c
c     ----- find if the file exists -----
c
      file1=ioflno(key1,nfile(un1),unitpt(un1),keylis,error)
      if (error.ne.0) then
         if (char.eq.'noerror') return
         call lnkerr('io system: error trying to copy non-existant '
     #               //'file '//key1//' from unit '//unit1)
      end if
      type=filtyp(file1)
      if (type(1:1).eq.'c') then
         type='character'
      else if (type(1:1).eq.'i') then
         type='integer'
      else if (type(1:1).eq.'r') then
         type='real'
      end if
      nb=eof(file1)-base(file1)
      size=iosize(type)
      length=nb/size
c
c     ----- find if the destination file exists -----
c
      file2=ioflno(key2,nfile(un2),unitpt(un2),keylis,error)
      if (error.eq.0) then
c
c        ----- check the length of the file we are copying to. if too
c              short, create a new file.
c
         nb2=end(file2)-base(file2)
         if (nb2.lt.nb) then
            error=1
            keylis(file2)='dead file'
         end if
      end if
c
      if (error.ne.0) then
c
c        ----- the destination file does not exist, so create it...
c
         newpos=0
         call iofile(type//' "'//key2//'" '//unit2,newpos,length)
c
c        ----- creating a new file may have changed the numbering -----
c
         file1=ioflno(key1,nfile(un1),unitpt(un1),keylis,error)
         file2=ioflno(key2,nfile(un2),unitpt(un2),keylis,error)
      end if
c
c     ----- the files exist -----
c
c
c     ----- check for constants -----
c
      if (filtyp(file1).eq.'c'.and.nb.le.len(cconst(1))) then
         cconst(file2)=cconst(file1)
         eof(file2)=base(file2)+nb
         return
      else if (filtyp(file1).eq.'i'.and.nb.eq.itobyt(1)) then
         iconst(file2)=iconst(file1)
         eof(file2)=base(file2)+nb
         return
      else if (filtyp(file1).eq.'r'.and.nb.eq.wptbyt(1)) then
         rconst(file2)=rconst(file1)
         eof(file2)=base(file2)+nb
         return
      end if
c
      ida1=base(file1)
      ida2=base(file2)
      ntogo=nb
c
c     ----- round of the buffer length to a multiple of 8192 bytes -----
c
      if (itobyt(lenbuf).lt.8192) then
         bite=itobyt(lenbuf)
      else
         bite=itobyt(lenbuf)/8192*8192
      end if
c
c     ----- read and write buffers till done -----
c
  400 continue
         nbytes=min(ntogo,bite)
         call ioget(unitno(un1),buffer,nbytes,ida1,junk)
         call ioput(unitno(un2),buffer,nbytes,ida2,junk)
c
         ida1=ida1+nbytes
         ida2=ida2+nbytes
c
         ntogo=ntogo-nbytes
      if (ntogo.gt.0) go to 400
c
c     ----- sort out pointers, etc -----
c
      eof(file2)=base(file2)+(eof(file1)-base(file1))
c
      return
c
c
      end
