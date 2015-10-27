*deck @(#)findcn.f	5.1  11/6/94
      subroutine findcn(maxap3,natoms,a,b,d,atmchg,npop,nset,
     $     ixyz,norder)
      implicit real*8(a-h,o-z)
c
c     this routine finds the highest order proper axis which is
c     coincident with cartesian axis ixyz.
c
c     an axis of order n will produce a number of "circular-sets" of
c     equivalent atoms (circular-set is more fully defined in
c     subroutine cirset).  furthermore, the population of each of
c     these sets must be an integer multiple of n.
c
      dimension a(maxap3,3), b(maxap3,3), d(maxap3,3), npop(1),
     $     nset(1),
     $     atmchg(1), t(3,3)
      data one,eight/1.0d+00,8.0d+00/
      save one,eight
c
c     call rtrace(6hfindcn,1)
      twopi = eight * atan(one)
c
c     cirset determines the populations of all circular sets which
c     are present.
c
      call cirset(maxap3,natoms,a,atmchg,ixyz,nset,npop,d,numset)
c
c     test the common multiples of the elements of npop in descending
c     order as possible orders for an axis of symetry.
c
      maxmul = 1
      do 20 i=1,numset
 20      maxmul = max(maxmul,npop(i))
         do 60 i=1,maxmul
            multst = maxmul - i + 1
            do 40 j=1,numset
               if (mod(npop(j),multst) .ne. 0) goto 60
 40         continue
            theta = twopi / float(multst)
            call rotate(maxap3,a,b,natoms,t,ixyz,theta)
            call equiv(maxap3,a,b,atmchg,natoms,itst)
            if (itst .eq. 0) goto 60
            norder = multst
            return
 60      continue
         norder = 1
         return
         end
