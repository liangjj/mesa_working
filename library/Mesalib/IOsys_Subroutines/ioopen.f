*deck @(#)ioopen.f	5.1  11/6/94
      subroutine ioopen(string,pos,numunt,filesz)
c
c***begin prologue     ioopen
c***date written       850125   (yymmdd)
c***revision date      910625   (yymmdd)
c
c  25 june 1991        rlm at lanl
c      on the sun, the ioinq function gives undependable results when
c      called with the 'open' query.  it is apparently ok with the 'exist'
c      query. the section of code below uses this function only to see
c      if some gopher is trying to open an iosys file with the same name
c      as some non-iosys file like input or output.
c      the ioinq function in mdutil/sun has been modified to always return
c      .false. if asked if a file is open, so beware of the section of
c      code here where it is called.  this is apparently the only call
c      with the 'open' query in all of mesa, it protects against
c      a fairly unlikely event, and it will be fixed in the next version
c      of the compiler so i decided not to modify this routine. 
c  10 july 1990        rlm at lanl
c      giving scratch files a name scr[unit number]
c      i.e. scr11, etc.
c   6 july 1987        pws at lanl
c      changing to a machine independent form -- or at least more
c      so. using the routines ioget, iogetc, ioput, ioputc, ioinq
c      iorm and ioopn for this purpose.
c
c   1 february 1987  pws at lanl
c      changing 'unitnm' and 'number'to character*256 variables to
c      allow for long file names to accomodate unix directory system.
c
c***keywords           iosys dependent routines
c
c***author             saxe, paul,    (lanl)
c***source             @(#)ioopen.f	5.1   11/6/94
c***purpose            to open an existing, or create a new, iosys unit
c
c***description        #
c
c
c***references
c
c***routines called    gettok (io)
c                      lnkerr (mdutil)
c                      ioget (mdutil)
c                      iogetc(mdutil)
c                      ioopn(mdutil)
c                      ioputc(mdutil)
c                      iorm(mdutil)
c                      dirsrt (io)
c
c   common blocks:     (none)
c
c***end prologue       ioopen
c
      implicit integer (a-z)
c
      parameter (maxfil=2000,maxun=20)
c
c
      character*(*) string,numunt
      character*16 idum(3)
      character*16 itoc
c
      character*80 extra
      character*80 tmplin
      character*256 unitnm(maxun),number
      character*32 keylis,key
      character*16 cconst,unlist,unityp,type,unit
      character*8  filtyp
      character*16 temp
      integer base,readpt,writpt,eof,end,iconst,unitpt,nfile
      integer un,file,nunits,unitno
      real*8 rconst
      logical align,locked,open
      logical ioinq
c
c
      common /ioqqq1/ keylis(maxfil),filtyp(maxfil),cconst(maxfil),
     #                unlist(maxun),type,key,unit,tmplin,unityp(maxun)
      common /ioqqq2/ base(maxfil),readpt(maxfil),writpt(maxfil),
     #                eof(maxfil),end(maxfil),iconst(maxfil),
     #                rconst(maxfil),unitpt(maxun),nfile(maxun),
     #                align(maxun),locked(maxun),unitno(maxun),
     #                un,file,nunits
      common /pointr/ wptr(128),tptr(128)
      common/io/inp,iout
c
c
      number=numunt
c
c     ----- parse the rest of 'string'. expecting:
c           'open <unit> (as [new,old,unknown,scratch]) (on) (ssd)'
c
      call gettok(unit,string,pos)
c
c     ----- find out if opening as new, old or unknown -----
c
      call gettok(type,string,pos)
      if (type.eq.'as') call gettok(type,string,pos)
c
c     ----- decide if this unit is already open -----
c
      do 1 i=1,nunits
         if (unlist(i).eq.unit) go to 9
    1 continue
c
c     ----- the following checks whether the file is assigned in
c           any way to the program, e.g. as output
c
      if (type.ne.'scratch') then
         open=ioinq(number,'open')
      else
         open=.false.
      end if
      if (open) go to 9
c
c     ----- the unit is not open, check the unit number -----
c
      if (type.ne.'scratch') then
         do 2 i=1,nunits
            if (unitnm(i).eq.number) go to 19
    2    continue
      end if
c
c     ----- find a unit number for this unit -----
c
      do 101 i=10,90
         do 100 un=1,nunits
            if (i.eq.unitno(un)) go to 101
  100    continue
         nounit=i
         go to 102
  101 continue
c
      call lnkerr('cant find a unit number from 10 to 90 available')
c
  102 continue
c
c     ----- all ok, so open this new unit -----
c
      nunits=nunits+1
      if (nunits.gt.maxun) then
         call lnkerr('io system: exceeded the maximum number of units')
      end if
      unlist(nunits)=unit
      if (type.eq.'scratch') then
         unityp(nunits)='scratch'
         number='scr'//itoc(nounit)
      else
         unityp(nunits)='permanent'
      end if
      unitnm(nunits)=number
      unitno(nunits)=nounit
      nfile(nunits)=0
c
c     ----- check if the file exists, and act accordingly -----
c
      open=ioinq(unitnm(nunits),'exist')
c
      if (type.eq.'new') then
         if (open) then
            call iorm(unitnm(nunits))
            open=ioinq(unitnm(nunits),'exist')
         end if
c
         if (open) then
            unit=unlist(nunits)
            call lnkerr('io system: error opening unit '//unit//
     #                  ' as new. it can''t be destroyed')
         end if
      else if (type.eq.'old') then
         if (.not.open) then
            unit=unlist(nunits)
            call lnkerr('io system: error opening unit '//unit//
     #                  ' as old...it does not exist')
         end if
      else if (type.eq.'unknown') then
         if (.not.open) then
            type='new'
         else
            type='old'
         end if
      else if (type.eq.'scratch') then
      else
         unit=unlist(nunits)
         call lnkerr('io system: error with type opening unit '//
     #                unit)
      end if
c
c        ----- check if this unit is to be on the ssd -----
c
      extra=' '
      call gettok(unit,string,pos)
      if (unit.eq.'on') call gettok(unit,string,pos)
      if (unit.eq.'ssd') extra='ssd'
      if (filesz.gt.0) then
         temp=itoc(filesz)
         do 345 i=1,len(temp)
            if (temp(i:i).ne.' ') go to 346
 345     continue
         i=1
 346     continue
         start=i
         extra(5:)='size='//temp(start:)
      end if
c
c
      call ioopn(nounit,unitnm(nunits),type,ierr,extra)
c
c
      if(ierr.ne.0) then
         call lnkerr('error opening '//unlist(nunits)//' (file: '//
     #        unitnm(nunits)//') as '//type)
      end if
c
      if (type.eq.'old') then
c
c     ----- unit is now opened and the first family member exists.
c           now need to set up tables, etc.
c
         call iogetc(nounit,idum,3*len(idum(1)),0,junk)
c
         if (idum(1).ne.'iosys') then
            unit=unlist(nunits)
            call lnkerr('io system: trying to open unsuitable unit '
     #                  //unit//'  directory corrupted')
         end if
c
         nfiles=ctoi(idum(2))
         nfile(nunits)=nfiles
         if (nunits.gt.1) then
            unitpt(nunits)=unitpt(nunits-1)+nfile(nunits-1)
         else
            unitpt(nunits)=0
         end if
c
         call dirsrt
c
c        ----- read in the directory info, bit by bit -----
c
         pt=ctoi(idum(3))
         file=unitpt(nunits)+1
         call iogetc(nounit,keylis(file),len(keylis(1))*nfiles,pt,pt)
         call ioget(nounit,base(file),itobyt(nfiles),pt,pt)
         call ioget(nounit,readpt(file),itobyt(nfiles),pt,pt)
         call ioget(nounit,writpt(file),itobyt(nfiles),pt,pt)
         call ioget(nounit,eof(file),itobyt(nfiles),pt,pt)
         call ioget(nounit,end(file),itobyt(nfiles),pt,pt)
         call iogetc(nounit,filtyp(file),len(filtyp(1))*nfiles,pt,pt)
         call iogetc(nounit,cconst(file),len(cconst(1))*nfiles,pt,pt)
         call ioget(nounit,iconst(file),itobyt(nfiles),pt,pt)
         call ioget(nounit,rconst(file),wptbyt(nfiles),pt,pt)
c
      else if (type.eq.'new'.or.type.eq.'scratch') then
c
c        ----- if a new unit, simply set up the empty file -----
c
         call dirsrt
c
         nfile(nunits)=1
         file=unitpt(nunits)+1
         keylis(file)='empty space file'
         base(file)=4096
         readpt(file)=0
         writpt(file)=0
         eof(file)=4096
         end(file)=4096
      else
         call lnkerr('illegal type opening unit: '//type)
      end if
c
c     ----- change the header word on the file so that if the
c           job crashes without writing the directory back,
c           the file cannot be reopened.
c
      idum(1)='changed'
      call ioputc(nounit,idum,3*len(idum(1)),0,junk)
c
c
      return
c
    9 continue
   19 continue
      call lnkerr('io system: duplicate open of unit '//unit)
c
c
      end
