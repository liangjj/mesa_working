*deck @(#)zconv.f	5.1  11/6/94
      subroutine zconv(action,nz,bl,alpha,beta)
c***begin prologue     zconv
c***date written       850601  yymmdd
c***revision date      yymmdd  yymmdd
c***keywords           z-matrix, bond lengths, bond angles, conversion
c***author             martin, richard (lanl)
c***source
c***purpose            converts bond lengths/angles from angstroms/degrees
c                      to bohrs/radians, and vice versa.
c***description
c     call zconv(action,nz,bl,alpha,beta)
c       action   either 'toborrad' or 'toangdeg'.
c       nz      the number of z-matrix entries.
c       bl      a vector(nz) containing the bond lengths.
c       alpha   a vector(nz) containing the bond angles.
c       beta    a vector(nz) containing the dihedral angles.
c***references
c***routines called    lnkerr(mdutil)
c***end prologue       zconv
      implicit integer(a-z)
      character action*8
      real*8 bl(nz),alpha(nz),beta(nz)
      real*8 toang,one,f45,conlen,angcon
      parameter (one=1.0d+00, f45=45.0d+00)
c
c
      call iosys('read real angstrom/bohr from rwf',-1,toang,0,' ')
      if(action.eq.'toborrad') then
         angcon=atan(one)/f45
         conlen=one/toang
      else if(action.eq.'toangdeg') then
         angcon=f45/atan(one)
         conlen=toang
      else
         call lnkerr(' unrecognized command in zconv: '//action)
      endif
c
c     convert the z-matrix.
      do 20 i=1,nz
         bl(i)    = bl(i)    * conlen
         alpha(i) = alpha(i) * angcon
         beta(i)  = beta(i)  * angcon
   20    continue
c
c
      return
      end
