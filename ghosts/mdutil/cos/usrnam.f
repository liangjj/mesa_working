*deck %W%  %G%
      subroutine usrnam(user)
c***begin prologue     usrnam
c***date written       850601  yymmdd
c***revision date      yymmdd  yymmdd
c***keywords           user, name
c***author             martin, richard (lanl)
c***source             %W%   %G%
c***purpose            retrieves the name of the person running the job.
c***description
c     call usrnam(user)
c       user    the users' name.  defaults to 'martin' if a production run.
c***references
c***routines called    usernoo(ctss),who(ctss),captlz(chr)
c***end prologue       usrnam
      implicit integer(a-z)
      character*(*) user
      character blank*1
      dimension ibuf(2000b)
      data blank/' '/
c
c
      user=blank
c      call userno(iuser)
c      call who(iuser,ibank,user,ibuf,len)
c      if(user.eq.'producti') user='martin'
      user='bhl'
      call locase(user,user)
c
c
      return
      end
