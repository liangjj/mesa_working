*deck %W%  %G%
      subroutine unqfil(file)
c***begin prologue     unqfil
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           files, input/output
c***author             saxe, paul (lanl)
c***source             %W%   %G%
c***purpose            generates a unique file name.
c***description
c                      call unqfil(file)
c                        file    on input, the file name to be made unique.
c                                on output, it has the form [suffix]file,
c                                where [suffix] is the suffix the job is
c                                running under.  used primarily to create
c                                scratch files which are hidden from the user.
c***references
c***routines called    suffix(ctss)
c***end prologue       unqfil
      implicit integer(a-z)
      character*(*) file
      character*8 unique,userno,account,dropfile
c
c
c     call suffix(unique,junk)
      call userinfo(userno,account,dropfile,unique)
      unique(2:8)=file
      file=unique
c
c
      return
      end
