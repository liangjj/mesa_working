*deck %W%  %G%
      subroutine iorm(file)
c
c***begin prologue     iorm
c***date written       870706   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           removal of a file, rm
c***author             saxe, paul (lanl)
c***source             %W%   %G%
c
c***purpose            to destroy or remove a file .
c
c***description        
c
c***references         
c
c***routines called    (none)
c
c***end prologue       iorm
c
      implicit integer (a-z)
c
      character*(*) file
c
      call drabsf(0,file,-1,ierr)
c
c
      return
      end
