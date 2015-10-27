*deck @(#)rzsub.f	5.1  11/6/94
      subroutine rzsub(file,nvar,names,values,intvec,fpvec)
c***begin prologue     rzsub
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           rwf, chk, z-matrix
c***author             martin, richard (lanl)
c***source
c***purpose            retrieves the arrays describing the z-matrix
c                      substitution information from an iosys file.
c***description
c                      call rzsub(file,nvar,names,values,intvec,fpvec)
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
c***iosys i/o                        unit unknown
c                      znames      character     read
c                      zvalues     real          read
c                      zintvec     integer       read
c                      zfpvec      real          read
c
c***routines called    iosys(io)
c***end prologue       rzsub
      implicit integer(a-z)
      real*8 values(nvar),fpvec(nvar)
      character*(*) names(*), file
      character*8 filnam
      dimension intvec(nvar)
c
c
c     read /zsubst/
      filnam=file
      if(nvar.ne.0) then
         call iosys('read character "variable names" from '//filnam,
     $               -1,0,0,names)
         call iosys('read real zvalues from '//filnam,-1,values,0,' ')
         call iosys('read integer zintvec from '//filnam,-1,intvec,0,
     $        ' ')
         call iosys('read real zfpvec from '//filnam,-1,fpvec,0,' ')
      endif
c
c
      return
      end
