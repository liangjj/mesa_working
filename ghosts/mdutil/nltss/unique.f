*deck %W%  %G%
      subroutine unique(name)
      character*8 name,account,userno,dropfil,suffix
c
      call dropfile(0)
c
      call userinfo(userno,account,dropfil,suffix)
      name(1:1)=suffix(1:1)
c
       return
       end
