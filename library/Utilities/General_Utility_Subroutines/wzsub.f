*deck @(#)wzsub.f	5.1  11/6/94
      subroutine wzsub(file,nvar,names,values,intvec,fpvec)
c***begin prologue     wzsub
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           rwf, chk, z-matrix
c***author             martin, richard (lanl)
c***source
c***purpose            writes the arrays describing the z-matrix
c                      substitution information to an iosys file.
c***description
c                      call wzsub(file,nvar,names,values,intvec,fpvec)
c                        file     character string giving the file to search
c                                 (e.g., 'rwf').
c                        nvar     the number of variables in the z-matrix.
c                        names    the variable names.
c                        values   the varaible values
c                        intvec
c                        fpvec
c
c***references
c
c***iosys i/o                            unit unknown
c                      znames        character    written
c                      zvalues       real         written
c                      zintvec       integer      written
c                      zfpvec        real         written
c
c***routines called    iosys(io)
c***end prologue       wzsub
      implicit integer(a-z)
      real*8 values(nvar), fpvec(nvar)
      character*(*) names(*), file
      character filnam*8
      dimension intvec(nvar)
c
c     write /zsubst/
      filnam=file
      if(nvar.ne.0) then
         lennms=nvar*16
         call iosys('write character "variable names" on '//filnam,
     $               lennms,0,0,names)
         call iosys('write real zvalues on '//filnam,nvar,values,0,' ')
         call iosys('write integer zintvec on '//filnam,nvar,intvec,
     #               0,' ')
         call iosys('write real zfpvec on '//filnam,nvar,fpvec,0,' ')
      endif
c
c
      return
      end
