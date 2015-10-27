*deck @(#)iorew.f	4.1  7/7/93
      subroutine iorew(string,pos)
c
c***begin prologue     iorew
c***date written       850125   (yymmdd)
c***revision date      860112   (yymmdd)
c***keywords           iosys dependent routines
c
c***author             saxe, paul,    (lanl)
c***source             @(#)iorew.f	4.1   7/7/93
c***purpose            to rewind the read and/or write pointers of an
c                        iosys file.
c***description        #
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
c***end prologue       iorew
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
c     ----- parse the string, which is of form:
c           'rewind <key> (on) <unit> (read)'
c                   (all)      (all)  (write)
c                                    (read-and-write)
c
      call gettok(key,string,pos)
      call gettok(unit,string,pos)
      if (unit.eq.'on') call gettok(unit,string,pos)
      call gettok(type,string,pos)
      if (type.eq.'eos') type='read-and-write'
c
c     ----- set up limits on units -----
c
      if (unit.eq.'all') then
         unmin=1
         unmax=nunits
      else
         un=iounit(unit,unlist,nunits,error)
         if (error.ne.0) then
            call lnkerr('io system: error in unit '//unit//
     #                  ' on rewind call')
         end if
         unmin=un
         unmax=un
      end if
c
c     ----- loop through the unit(s), rewinding as needed -----
c
      do 100 un=unmin,unmax
         pt=unitpt(un)
         nfiles=nfile(un)
         if (key.eq.'all') then
            do 4 file=pt+1,pt+nfiles
               if (type.eq.'read'.or.type.eq.'read-and-write')
     #             readpt(file)=0
               if (type.eq.'write'.or.type.eq.'read-and-write')
     #             writpt(file)=0
    4       continue
         else
c
c           ----- find key in list -----
c
            file=ioflno(key,nfiles,pt,keylis,error)
            if (error.ne.0) then
               call lnkerr('io system: cannot find file '//key//
     #                     ' on '//unit//' to rewind')
            end if
c
            if (type.eq.'read'.or.type.eq.'read-and-write')
     #          readpt(file)=0
            if (type.eq.'write'.or.type.eq.'read-and-write')
     #          writpt(file)=0
         end if
  100 continue
c
c
      return
      end
