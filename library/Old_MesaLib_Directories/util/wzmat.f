*deck @(#)wzmat.f	5.1  11/6/94
      subroutine wzmat(file,nz,nvar,ianz,iz,bl,alpha,beta,lbl,lalpha,
     $                 lbeta)
c***begin prologue     wzmat
c***date written       850601    (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           rwf, chk, z-matrix
c***author             martin, richard (lanl)
c***source
c***purpose            writes the arrays describing the z-matrix to an
c                      iosys file.
c***description
c                      call wzmat(file,nz,nvar,ianz,iz,bl,alpha,beta,lbl,
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
c***iosys i/o                           unit unknown
c                      zian      integer       written
c                      ziz       integer       written
c                      zbl       real          written
c                      zalpha    real          written
c                      zbeta     real          written
c                      zlbl      integer       written
c                      zlalpha   integer       written
c                      zlbeta    integer       written
c
c***routines called    iosys(io)
c***end prologue       wzmat
      implicit integer(a-z)
      character*(*) file
      character*8 filnam
      real*8 bl(nz),alpha(nz),beta(nz)
      integer ianz(nz),iz(4,nz),lbl(nz),lalpha(nz),lbeta(nz)
c
c
c     write the /zmat/.
      filnam=file
      call iosys('write integer zian on '//filnam,nz,ianz,0,' ')
      call iosys('write integer ziz on '//filnam,4*nz,iz,0,' ')
      call iosys('write real zbl on '//filnam,nz,bl,0,' ')
      call iosys('write real zalpha on '//filnam,nz,alpha,0,' ')
      call iosys('write real zbeta on '//filnam,nz,beta,0,' ')
      call iosys('write integer zlbl on '//filnam,nz,lbl,0,' ')
      call iosys('write integer zlalpha on '//filnam,nz,lalpha,0,' ')
      call iosys('write integer zlbeta on '//filnam,nz,lbeta,0,' ')
c
c
      return
      end
