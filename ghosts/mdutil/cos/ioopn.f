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
c
      character*(*) file
      character*(*) status
      character*(*) extra
      character*(8) iounq
      character*(8) scrnam
      integer unit
      integer intkey
      logical logkey
      character*8 lcname
      character*8 cname
      logical iiii
      logical ifdnt
c
c     ----- open the unit as a cftlib random, familied unit -----
c
      if (status.eq.'scratch') then
c
c        ----- get a unique file-name -----
c
         scrnam=iounq()
         cname=scrnam
      else
         cname=file
      end if
      call locase(cname,cname)
      lcname=cname
      do 58 i=1,len(lcname)
         if(lcname(i:i).eq.' ') lcname(i:i)=char(0)
   58 continue
c
c     ----- set up an alias name of the form 'ftnn' -----
c
      write (scrnam,60) unit
 60   format ('ft',i2)
      if (scrnam(3:3).eq.' ') scrnam(3:3)='0'
c      do 65 i=5,8
c         scrnam(i:i)=char(0)
c 65   continue
c
c     ----- assign the cos file with the maximum size limit  
c           but first make sure it wasn't assigned in the jcl ----
c
c
      iiii=ifdnt(lcname)
      iiii=.false.
c 
      if (.not.iiii) then
c
c        ----- put the file on the ssd if requested -----
c
         if (logkey(extra,'ssd',.false.,' ')) then
            iass=1
            call assign(iass,'dn'l,cname,'dv'l,'ssd-0-20'l,
     #               'u'l,'lm'l,'320000'l,'a'l,scrnam) 
         else
            call assign(iasss,'dn'l,cname,'u'l,'lm'l,'320000'l,
     $                  'a'l,scrnam)
         end if
         if(iasss.ne.0) then
            call lnkerr('io system: unsuccessful cos assign') 
         endif
      endif
      call wopen(unit,20,0)
c
c
      return
      end
