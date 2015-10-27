*deck @(#)findv.f	5.1  11/6/94
      subroutine findv(maxap3,a,b,d,natoms,npop,nset,atmchg,itst)
      implicit real*8(a-h,o-z)
c
c     this routine tests for a set of norder vertical planes.  it
c     is assumed that the principal axis is aligned with the
c     the cartesian z axis.  if a set of planes is found, it leaves
c     one of them coincident with the yz cartesian plane.
c
      common/tol/ toler,tol2
      dimension a(maxap3,3), b(maxap3,3), d(maxap3,3), atmchg(1),
     $     npop(1), nset(1), t(3,3)
      data half,one,two/0.5d+00,1.0d+00,2.0d+00/
      save half,one,two
c
c     call rtrace(6hfindv ,1)
      numatm = natoms + 3
      halfpi = two * atan(one)
c
c     look for a circular-set of atoms.  a vertical mirror must
c     pass through the point midway between any atom in the set
c     and one other atom in the set.
c
      call cirset(maxap3,natoms,a,atmchg,3,nset,npop,d,numset)
c
      iset = 1
      iattop = natoms - 1
      do 20 iat=1,iattop
         if (nset(iat) .ne. iset) goto 20
         goto 40
 20   continue
      itst = 0
      return
 40   continue
      j1 = iat + 1
      do 60 jat=j1,natoms
         if (nset(jat) .ne. iset) goto 60
         x = (a(iat,1)+a(jat,1)) * half
         y = (a(iat,2)+a(jat,2)) * half
         if (abs(x).gt.toler .or. abs(y).gt.toler) goto 50
         x = half * a(jat,1)
         y = half * a(jat,2)
 50      theta = halfpi
         if (abs(y) .gt. toler) theta = - atan(x/y)
         call rotate(maxap3,a,b,numatm,t,3,theta)
         call reflct(maxap3,b,d,natoms,t,1)
         call equiv(maxap3,b,d,atmchg,natoms,itst)
         if (itst .eq. 0) goto 60
         call move(maxap3,b,a,3*numatm)
         itst = 1
         return
 60   continue
      itst = 0
      return
      end
