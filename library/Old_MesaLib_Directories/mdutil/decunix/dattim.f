*deck @(#)dattim.f	1.1  11/30/90
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
c***source             @(#)dattim.f	1.1   11/30/90
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
      character*9 dat
      character*8 tim  
c
c     get the date.  the format is 'day mon yy hh:mm:ss'
c                                   123456789012345678901234
c
      call date(dat)
      call time(tim)
c
c     put it all together.
c
      call locase(dat,dat)
      dt=dat(1:3)//dat(4:7)//'19'//dat(8:9)//' '//tim(1:8)
c
c
      return
      end
