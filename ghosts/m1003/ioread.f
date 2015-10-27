*deck %W%  %G%
      subroutine ioread(string,pos,nwords,array,rarray,offset,c1)
c
c***begin prologue     ioread
c***date written       850125   (yymmdd)
c***revision date      860112   (yymmdd)
c***keywords           iosys dependent routines
c
c***author             saxe, paul,    (lanl)
c***source             %W%   %G%
c***purpose            to read data from an iosys file.
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
c***end prologue       ioread
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
c
c
c     ----- parse the rest of 'string', expecting:
c           'read (character) <key> (from) <unit> (after rewinding)
c                 (integer)                     (without rewinding)
c                 (real)                        (asynchronously)'
c
      call gettok(key,string,pos)
      if (index('char inte real',key(1:4)).eq.0) then
         type='unknown'
      else
         type=key
         call gettok(key,string,pos)
      end if
c
      call gettok(unit,string,pos)
      if (unit.eq.'from') call gettok(unit,string,pos)
      asynch=index(string,'asynch').ne.0
      rewind=index(string,'without rewind').eq.0
c
c     ----- find the unit from its name -----
c
      un=iounit(unit,unlist,nunits,error)
      if (error.ne.0) then
         call lnkerr('io system: error trying to read '//key//' from'
     #               //' non-existant unit: '//unit)
      end if
c
c     ----- find if the file exists -----
c
      file=ioflno(key,nfile(un),unitpt(un),keylis,error)
      if (error.ne.0) then
         call lnkerr('io system: error trying to read non-existant '
     #               //'file '//key//' from unit '//unit)
      end if
c
c     ----- if the operation string gives a type, check against file
c
      if (type.ne.'unknown'.and.type(1:1).ne.filtyp(file)) then
         call lnkerr('read request for '//key//' on '//unit//
     #               ' has an incompatible type ('//type(1:9)//
     #               '). should be '//filtyp(file))
      end if
c
c     ----- the file exists -----
c
      size=iosize(filtyp(file))
      if (nwords.eq.-1) then
         nb=eof(file)-base(file)
      else
         nb=nwords*size
      end if
      if (filtyp(file).eq.'c') then
         if (nwords.eq.0) then
            nb=len(c1)
         end if
      end if
c
      maxbyt=eof(file)-base(file)
c
c     ----- check for constants -----
c
      if(rewind) then
      if (filtyp(file).eq.'c'.and.maxbyt.le.len(cconst(1))) then
         if (nb.gt.maxbyt) then
            call lnkerr('character constant '//key//' on '//unit//
     #                  ' is shorter than you the read you have '//
     #                  'requested')
         end if
         c1=cconst(file)(1:nb)
         return
      else if (filtyp(file).eq.'i'.and.maxbyt.eq.itobyt(1)) then
         if (nb.gt.maxbyt) then
            call lnkerr('integer constant '//key//' on '//unit//
     #                  ' is shorter than you the read you have '//
     #                  'requested')
         end if
         array(1)=iconst(file)
         return
      else if (filtyp(file).eq.'r'.and.maxbyt.eq.wptbyt(1)) then
         if (nb.gt.maxbyt) then
            call lnkerr('real constant '//key//' on '//unit//
     #                  ' is shorter than you the read you have '//
     #                  'requested')
         end if
         rarray(1)=rconst(file)
         return
      end if
      end if
c
c     ----- check for a possible rewind of the file -----
c
      if (rewind) readpt(file)=0
c
c     ----- and handle the offset, if it exsits -----
c
      readpt(file)=readpt(file)+offset*size
      if (readpt(file).lt.0) then
         call lnkerr('io system: reading '//key//' on '//unit//
     #               ' with bad offset')
      end if
c
      ida=base(file)+readpt(file)
      if (ida+nb.gt.eof(file)) then
         call lnkerr('io system: reading '//key//' on '//unit//
     #               'exceeds file length')
      end if
c
      if ((ida+itobyt(1)-1)/itobyt(1)*itobyt(1).ne.ida) then
         call lnkerr('ioread is trying to read '//key//' on '//unit//
     #               ' starting at a non-integral number of words'//
     #               '......fatal')
      end if
c
      if (type.eq.'character') then
         call iogetc(unitno(un),c1,nb,ida,junk)
      else
         call ioget(unitno(un),array,nb,ida,junk)
      end if
c
c     ----- check if asynchronous io requested -----
c
      if (asynch) go to 7
      if (unitwt(unitno(un))) 7,900,910
    7 continue
      readpt(file)=readpt(file)+nb
c
      return
c
  900 continue
  910 continue
      call lnkerr('io system: i/o error reading '//key//' from '//
     #             unit)
c
c
      end
