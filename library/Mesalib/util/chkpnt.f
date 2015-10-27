*deck @(#)chkpnt.f	5.1  11/6/94
      subroutine chkpnt(action)
c***begin prologue     chkpnt
c***date written       850601  (yymmdd)
c***revision date      930120  (yymmdd)
c
c   20 january  1993   rlm at lanl
c       units 0 and 7 are now reserved for standard error.
c       in order to avoid a conflict, the data file is opened as unit 1.
c    1 december 1986   pws at lanl
c       changing iosys open to character.
c
c***keywords           checkpoint
c***author             martin, richard (lanl)
c***source
c***purpose            checkpoints the calculation.
c***description
c                      call chkpnt(action)
c                        action   what to do, character*(*).
c                                 if 'save', move from rwf to chk.
c                                 if 'restore', move from chk to rwf.
c
c***references
c***routines called    iosys(io)
c***end prologue       chkpnt
c
      implicit integer(a-z)
c
      parameter (lenbuf=4096)
c
      character*(*) action
      character sink*3,source*3,actdum*8
      character*128 namchk,namdat
      character*32 file
      character*16 cjunk
      integer a(lenbuf)
      logical positn
c
c
      if(action.eq.'save') then
         sink='chk'
         source='rwf'
      else if(action.eq.'restore') then
         sink='rwf'
         source='chk'
      else
         actdum=action
         call lnkerr('unrecognized command in chkpnt: '//actdum)
      endif
c
c     ----- open the data file for the file-names to transfer -----
c
      numdat=1
      call iosys('read character "data filename" from rwf',
     $     0,0,0,namdat)
      open (unit=numdat,file=namdat,status='old')
c
c     ----- position the data to the file names -----
c
      if (.not.positn('$checkpoint',cjunk,numdat)) then
         call lnkerr('can''t find the checkpoint data on mesadat')
      end if
c
c     open the checkpoint file.
c
      call iosys('read character "checkpoint filename" from rwf',
     $     0,0,0,namchk)
      call iosys('open chk as unknown',1000000,0,0,namchk)
c
c     checkpoint molecular parameters.
c
      do 10 i=1,100000
         read (numdat,1,end=20) file
 1       format (a32)
         if (file(1:1).eq.'$') go to 20
         if (file.ne.' ') then
            call iosys('copy "'//file//'" from '//source//' to '//
     $           sink,lenbuf,a,0,'noerror')
         end if
 10   continue
 20   continue
c
c
c
      close (unit=numdat)
      call iosys('close chk',namchk,0,0,' ')
c
c
      return
      end
