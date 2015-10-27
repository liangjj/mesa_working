*deck @(#)rzmat.f	5.1  11/6/94
      subroutine rzmat(file,nz,nvar,ianz,iz,bl,alpha,beta,lbl,lalpha,
     $                 lbeta)
c***begin prologue     rzmat
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           rwf, chk, z-matrix
c***author             martin, richard (lanl)
c***source
c***purpose            retrieves the arrays describing the z-matrix from an
c                      iosys file.
c***description
c                      call rzmat(file,nz,nvar,ianz,iz,bl,alpha,beta,lbl,
c                                 lalpha,lbeta)
c                        file     character string giving the file to search
c                                 (e.g., 'rwf').
c                        nz       the number of entries in the z-matrix.
c                        nvar     the number of variables in the z-matrix.
c                        ianz     the z-matrix atomic charges.
c                        iz       an array describing the z-matrix entries.
c                        bl       the bond lengths.
c                        alpha    the bond angles.
c                        beta     the dihedral angles.
c                        lbl      symbolic bond lengths.
c                        lalpha   symbolic bond angles.
c                        lbeta    symbolic dihedral angles.
c
c***references
c
c***iosys i/o                   unit unknown
c                      zian      integer       read
c                      ziz       integer       read
c                      zbl       real          read
c                      zalpha    real          read
c                      zbeta     real          read
c                      zlbl      integer       read
c                      zlalpha   integer       read
c                      zlbeta    integer       read
c
c***routines called    iosys(io)
c***end prologue       rzmat
      implicit integer(a-z)
      character*(*) file
      character*8 filnam
      real*8 bl(nz),alpha(nz),beta(nz)
      integer ianz(nz),iz(4,nz),lbl(nz),lalpha(nz),lbeta(nz)
c
 
c
c     read the /zmat/.
      filnam=file
      call iosys('read integer zian from '//filnam,-1,ianz,0,' ')
      call iosys('read integer ziz from '//filnam,-1,iz,0,' ')
      call iosys('read real zbl from '//filnam,-1,bl,0,' ')
      call iosys('read real zalpha from '//filnam,-1,alpha,0,' ')
      call iosys('read real zbeta from '//filnam,-1,beta,0,' ')
      call iosys('read integer zlbl from '//filnam,-1,lbl,0,' ')
      call iosys('read integer zlalpha from '//filnam,-1,lalpha,0,' ')
      call iosys('read integer zlbeta from '//filnam,-1,lbeta,0,' ')
c
c
      return
      end
