*deck @(#)iorm.f	5.1  11/6/94
      subroutine iorm(file)
c
c***begin prologue     iorm
c***date written       870706   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           removal of a file, rm
c***author             saxe, paul (lanl)
c***source             @(#)iorm.f	5.1   11/6/94
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
      character*130 tfile
c
c     junk=unlink(file)
      tfile='rm '//file(1:128)
      call system(tfile)
c
c
      return
      end
