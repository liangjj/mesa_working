*deck @(#)ordn.f	5.1  11/6/94
      subroutine ordn(maxap3,a,b,aset,atmchg,npop,
     $     nset,natoms,norder,dump)
      implicit real*8(a-h,o-z)
c
c     this routine orients symmetric top molecules in the point
c     groups dn, dnd, dnh, and cnv.  in dn and dnh molecules a
c     c2 axis is coincident with the cartesian y axis.  in dnd and
c     c2v molecules there is a vertical plane coincident with the
c     yz plane.  the molecule is rotated by pi/norder about
c     the z axis until one of the following occurs:
c     1-- the projection of the key atom on the y axis is a maximum.
c     2-- if two orientations give the greatest projection on y, the
c     one where the key atom has a positive x coordinate is chosen.
c     spherical top molecules in the point groups t, td, th, o, oh,
c     and ih are also oriented here.  the proper axis has already been
c     oriented with z.
c
      integer orkey
      logical dump
      common/tol/ toler,tol2
      common/io/inp,iout
      dimension a(maxap3,3), b(maxap3,3), aset(maxap3,3), atmchg(1),
     $     npop(1), nset(1), t(3,3)
      data zero,one,four/0.0d+00,1.0d+00,4.0d+00/
      save zero,one,four
 1000 format(' ordn-- key atom:',i3)
c
c     call rtrace(6hordn  ,1)
      pi = four * atan(one)
      numatm = natoms + 3
      theta = pi / norder
      key = orkey(maxap3,natoms,a,atmchg,nset,npop,aset)
      if(dump) write(iout,1000) key
      curpy = a(key,2)
      curpx = a(key,1)
c
c     rotate the molecule so that the key atom has a positive
c     y coordinate.
c
 20   continue
      if (curpy .gt. zero) goto 40
      call rotate(maxap3,a,b,numatm,t,3,pi)
      call move(maxap3,b,a,numatm)
      curpx = a(key,1)
      curpy = a(key,2)
c
c     search for the maximum curpy.
c
 40   continue
      itry = 0
      direct = - sign(one,curpx)
 60   phi = direct * theta
      call rotate(maxap3,a,b,numatm,t,3,phi)
      itry = itry + 1
      if (abs(curpy-b(key,2)) .lt. toler) goto 100
      if (curpy .gt. b(key,2)) goto 80
      call move(maxap3,b,a,numatm)
      curpy = a(key,2)
      if (abs(a(key,1)) .lt. toler) goto 200
      goto 60
 80   continue
      if (itry .gt. 1) goto 200
      direct = - direct
      goto 60
 100  continue
      if (b(key,1) .lt. zero) goto 200
 120  call move(maxap3,b,a,numatm)
c
c     planar molecules contained in the xy plane have not been
c     completely specified.
c
 200  continue
      call orptst(maxap3,a,natoms,ixyz)
      if (ixyz .eq. 3) call oraxis(maxap3,a,b,natoms,atmchg,1)
      return
      end
