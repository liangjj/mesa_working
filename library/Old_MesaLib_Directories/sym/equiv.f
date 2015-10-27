*deck @(#)equiv.f	5.1  11/6/94
      subroutine equiv(maxap3,a,b,atmchg,natoms,itst)
      implicit real*8(a-h,o-z)
c
c     itst is set to 1 if the two molecular orientations in a and
c     in b are equivalent.  otherwise itst is set to 0.
c
      common/tol/ toler,tol2
      dimension a(maxap3,3), b(maxap3,3), atmchg(1)
c
c     call rtrace(6hequiv ,1)
      do 100 iat=1,natoms
         do  80 jat=1,natoms
            if (abs(atmchg(iat)-atmchg(jat)) .gt. tol2) goto 80
            do  60 ixyz=1,3
               test = a(iat,ixyz) - b(jat,ixyz)
               if (abs(test) .gt. toler) goto 80
 60         continue
            goto 100
 80      continue
         itst = 0
         return
 100  continue
      itst = 1
      return
      end
