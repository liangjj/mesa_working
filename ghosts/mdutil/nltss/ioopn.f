*deck %W%  %G%
      subroutine ioopn(unit,file,status,iostat,extra)
c
c***begin prologue     ioopn
c***date written       870706   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           iosys dependent routine
c***author             saxe, paul (lanl)
c***source             %W%   %G%
c
c***purpose            to open a unit to a disc file
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       ioopn
c
      implicit integer (a-z)
c..bhl
      common /bhl/ maxwrd(100)
      common /blname/ bhlnam(100)
      character*8 bhlnam
c..bhl
c
      character*(*) file
      logical query
      character*(*) status
      character*(*) extra
      character*8 iounq
      character*8 scrnam
      dimension list(15)
      integer unit
      integer intkey
      logical logkey
c
c set list for ssd on nltss using fortlib calls
c 41b is the physical unit for the ssd
c
      do 2 i=1,15
      list(i)=0
  2   continue
      list(3)=4
      list(4)=3
      list(6)=61440
      do 3 j=7,10
      list(j)=-1
  3   continue
      list(15)=41b
c
c     ----- open the unit as a cftlib random, familied unit -----
c
      size=intkey(extra,'size',262144,' ')
c
      if (status.eq.'scratch') then
c
c        ----- get a unique file-name -----
c
         scrnam=iounq()
c        call assign(unit,scrnam,0)
c
       if (logkey(extra,'ssd',.false.,' ')) then
        call eqvfix(list(1),scrnam)
        if(list(2).ne.0) call lnkerr(' ioopn: ssd error ')
        list(2)=size
        call createl(unit,list,15)
       else
        call create(unit,scrnam,4,size)
       end if
c
      bhlnam(unit)=scrnam
      maxwrd(unit)=size
c
      else
c
c        call assign(unit,file,0)
      inquire(file=file,exist=query)
c
      if(query) then
       call open(unit,file,4,len)
       if(len.eq.-1) call lnkerr(' cant open file ')
       bhlnam(unit)=file
       maxwrd(unit)=len
      else
       if (logkey(extra,'ssd',.false.,' ')) then
        call eqvfix(list(1),file)
        if(list(2).ne.0) call lnkerr(' ioopn: ssd error ')
        list(2)=size
        call createl(unit,list,15)
       else
        call create(unit,file,4,size)
       end if
      bhlnam(unit)=file
      maxwrd(unit)=size
      end if
c
      end if
c
c     ----- put the file on the ssd if requested -----
c
c     if (logkey(extra,'ssd',.false.,' ')) then
c        call setssd(unit)
c     end if
c
c     ----- set the family size and wait parameters -----
c
c     call famsiz(unit,size)
c     call famwait(unit,1)
c
c
      return
      end
