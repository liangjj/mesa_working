*deck %W%  %G%
      subroutine dattim(dt)
c***begin prologue     dattim
c***date written       850601  (yymmdd)
c***revision date      860123  (yymmdd)
c***cos version
c***keywords           date, time
c***author             martin, richard (lanl)
c***source             mdutil
c***purpose            returns the current date and time of day.
c***description
c                      call dattim(dt)
c                        dt      the date and time (character*20)
c                                  dd-mmm-yyyy hh mm ss
c
c***references
c***routines called    date(cos), clock(cos)
c***end prologue       dattim
      implicit integer(a-z)
      character dt*(*)
      character today*8, time*8, nmonth(12)*2, amonth(12)*3
      character year*4, month*3, day*2, nummon*2
      data amonth/'jan','feb','mar','apr','may','jun',
     $            'jul','aug','sep','oct','nov','dec'/
      data nmonth/'01','02','03','04','05','06',
     $            '07','08','09','10','11','12'/
c
c
c     get the date.  the format is 'mm/dd/yy'
      call date(today)
      year='19'//today(7:8)
      nummon=today(1:2)
c
c     change from a numerical month to something more readable.
      do 100 i=1,12
         if(nummon.eq.nmonth(i)) then
            month=amonth(i)
         endif
  100 continue
      day=today(4:5)
c
c     get the time.  the format is 'hh mm ss'.
      call clock(time)
c
c     put it all together.
      dt=day//'-'//month//'-'//year//' '//time
c
c
      return
      end
