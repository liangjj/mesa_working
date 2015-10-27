*deck @(#)usrnam.f	5.1   11/6/94
      subroutine usrnam(user)
c
c***begin prologue     usrnam
c***date written       850601  yymmdd
c***revision date      870207  yymmdd
c
c   7 february 1987   pws at lanl
c       creating bsd 4.2 unix version for sun 3/50 and 3/160 workstations
c
c***keywords           user, name
c***author             martin, richard (lanl)
c***source             @(#)usrnam.f	5.1   11/6/94
c***purpose            retrieves the name of the person running the job.
c***description
c     call usrnam(user)
c       user    the users' name.  defaults to 'martin' if a production run.
c***references
c***routines called    
c***end prologue       usrnam
c
      implicit integer(a-z)
      character*(*) user
c
      call getlog(user)
c     call locase(user(1:1),user(1:1))
c
c
      return
      end
