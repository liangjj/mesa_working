*deck @(#)dattim.f	5.1  11/6/94
      subroutine dattim(dt)
c***begin prologue     dattim
c***date written       850601  (yymmdd)
c***revision date      870207  (yymmdd)
c
c   7 february 1987  pws at lanl
c      rewriting functionality for bsd 4.2 unix on sun 3/50 and 3/160
c      workstations.
c
c***keywords           date, time
c***author             martin, richard (lanl)
c***source             @(#)dattim.f	5.1   11/6/94
c***purpose            returns the current date and time of day.
c***description
c                      call dattim(dt)
c                        dt      the date and time (character*20)
c                                  dd-mmm-yyyy hh mm ss
c
c***references
c***routines called    
c***end prologue       dattim
c
      implicit integer(a-z)
      character dt*(*)
      integer timarr(9)
      integer dum,clock
      character*1 colon
      character*2 itoc
      character*3 month(12)
c
      data month/'jan','feb','mar','apr','may','jun','jul','aug','sep',
     $           'oct','nov','dec'/
      data colon/':'/
      save month,colon
c
c     get the date.  the format is 'day mon dd hh:mm:ss yyyy'
c                                   123456789012345678901234
c
      dum=time(clock)
      call movdt(%val(localtime(clock)),timarr)
c
c     put it all together.
c
      dt=itoc(timarr(4))//'-'//month(timarr(5)+1)//'-'
     $   //'19'//itoc(timarr(6))//' '//
     $   itoc(timarr(3))//':'//itoc(timarr(2))//':'//itoc(timarr(1))
c
c
      return
      end
