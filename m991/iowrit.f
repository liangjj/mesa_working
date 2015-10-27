*deck @(#)iowrit.f	5.1  11/6/94
      subroutine iowrit(string,pos,nwords,array,rarray,offset,c1)
c***begin prologue     iowrit
c***date written       850125   (yymmdd)
c***revision date      919818   (yymmdd)
c   18 august   1991   bhl at llnl
c      fixing this routine so that when the phrase "without rewinding'
c      is in the string one word writes actually go to disk, as opposed to
c      being stored in the constants table. see also ioread.f
c***keywords           iosys dependent routines
c
c***author             saxe, paul,    (lanl)
c***source             @(#)iowrit.f	5.1   11/6/94
c***purpose            to write data to an iosys file.
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
c                      ioput (cftlib)
c                      unitwt (io)
c
c   common blocks:     (none)
c
c***end prologue       iowrit
c
      implicit integer (a-z)
c
      character*(*) string,c1
      real*8 rarray(*)
      integer array(*)
      logical asynch,rewind
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
      common/io/ inp, iout         
c
c
c     ----- parse the rest of string, expecting:
c           'write character <key> (to) <unit> (asynchronously)'
c                  integer         (on)      (after rewinding)
c                  real                       (without rewinding)
c
      call gettok(type,string,pos)
c
c     ----- check that the type is legitimate -----
c
      if (index('char inte real',type(1:4)).eq.0) then
         call lnkerr('illegal or missing type (character, integer '//
     #               'or real) for write: '//tmplin)
      end if
c
      call gettok(key,string,pos)
      call gettok(unit,string,pos)
      if (unit.eq.'to'.or.unit.eq.'on') call gettok(unit,string,pos)
      asynch=index(string,'asynch').ne.0
      rewind=index(string,'without rewind').eq.0
c
c     ----- find the unit from its name -----
c
      un=iounit(unit,unlist,nunits,error)
      if (error.ne.0) then
         call lnkerr('io system: trying to write '//key//' on non-'//
     #               'existant unit:'//unit)
      end if
c
c     ----- find if the file exists -----
c
      file=ioflno(key,nfile(un),unitpt(un),keylis,error)
c
      if (error.ne.0) then
c
c        ----- it doesn't exist, so create the file ... -----
c
         newpos=0
         nw=nwords
         if (type.eq.'character'.and.nw.eq.0) nw=len(c1)
         call iofile(type//' "'//key//'" '//unit,
     #               newpos,nw)
c
c     #              ,maxfil,maxun,nunits,
c     #               unlist,unitpt,nfile,keylis,base,readpt,writpt,
c     #               eof,end,align,locked)
         file=ioflno(key,nfile(un),unitpt(un),keylis,error)
         if (error.ne.0) then
            call lnkerr(' io system: cannot create a file in iowrit')
         end if
      end if
c
c     ----- check that the type of the file is the same as that
c           in the operation string .....
c
      if (type(1:1).ne.filtyp(file)) then
         call lnkerr('a request to write '//key//' on '//unit//
     #               ' has a type ('//type(1:9)//') inconsistent '//
     #               'with a previous write: ('//filtyp(file))
      end if
c
c     ----- now the file exists, one way or another, so write....
c
      size=iosize(filtyp(file))
      if (nwords.eq.-1) then
         nb=eof(file)-base(file)
      else
         nb=nwords*size
      end if
      if (type.eq.'character') then
         if (nwords.eq.0) then
            nb=len(c1)
         end if
      end if
c
c     ----- look for a constant to slip into the directory -----
c
c     ----- note that this is not done if the phrase 'without rewinding'
c           appears in the iosys call.
      maxbyt=eof(file)-base(file)
      if(rewind) then
         writpt(file)=0
         if (type.eq.'character'.and.maxbyt.le.len(cconst(1)).and.
     $       nb.le.len(cconst(1))) then
            cconst(file)=c1
            eof(file)=base(file)+nb
            return
         else if (type.eq.'integer'.and.maxbyt.le.itobyt(1).and.
     $           nb.eq.itobyt(1)) then
            iconst(file)=array(1)
            eof(file)=base(file)+nb
            return
         else if (type.eq.'real'.and.maxbyt.le.wptbyt(1).and.
     $           nb.eq.wptbyt(1)) then
            rconst(file)=rarray(1)
            eof(file)=base(file)+nb
            return
         end if
      end if
c
c     ----- and handle the offset, if it exists -----
c
      writpt(file)=writpt(file)+offset*size
      write(iout,*) 'size= ',size
      write(iout,*) 'offset = ',offset
      write(iout,*) 'writpt = ',writpt(file)
      if (writpt(file).lt.0) then
         call lnkerr('io system: writing '//key//' on unit '//unit
     #                //' with bad offset')
      end if
c
      ida=base(file)+writpt(file)
      write(iout,*) 'base = ',base(file)
      write(iout,*) 'writpt = ', writpt(file)
      write(iout,*)  'ida= ',ida
      if (ida+nb.gt.end(file).and.end(file).gt.0) then
         call lnkerr('io system: writing '//key//' on '//unit//
     #               ' exceeds file length')
      end if
c
      if ((ida+3)/4*4.ne.ida) then
         write(iout,*) 'ida = ',ida
         call lnkerr('iowrite is trying to write '//key//' on '//unit//
     #               ' starting at a non-integral number of words'//
     #               '......fatal')
      end if
c
      if (type.eq.'character') then
         call ioputc(unitno(un),c1,nb,ida,junk)
      else
         call ioput(unitno(un),array,nb,ida,junk)
      end if
c
c     ----- check if asynchronous io requested -----
c
      if (asynch) go to 7
      if (unitwt(unitno(un))) 7,900,910
    7 continue
      writpt(file)=writpt(file)+nb
      if (base(file)+writpt(file).gt.eof(file)) eof(file)=
     #                                    base(file)+writpt(file)
c
      return
c
  900 continue
  910 continue
      call lnkerr('io system: i/o error while writing '//key//' on '
     #            //unit)
c
c
      end
