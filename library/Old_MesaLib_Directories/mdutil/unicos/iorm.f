*deck @(#)iorm.f	5.1  11/6/94
      subroutine iorm(file)
c
c***begin prologue     iorm
c***date written       870706   (yymmdd)
c***revision date      910617   (yymmdd)
c
c   17 june   1991     rlm at lanl
c      unicos version
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
c
      open(1,file=file)
      close(1,status='delete')
c
      return
      end
