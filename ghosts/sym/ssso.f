*deck ssso
      subroutine ssso(maxap3,natoms,iatflg,a,itst)
      implicit real*8(a-h,o-z)
c
c     flag an atom at the origin.
c
      common/tol/ toler,tol2
      dimension iatflg(1), a(maxap3,3)
c
c     call rtrace(6hssso  ,1)
      itst = 0
      do 40 iat=1,natoms
         if (iatflg(iat) .ne. 0) goto 40
         do 20 ixyz=1,3
            if (abs(a(iat,ixyz)) .gt. toler)  goto 40
 20      continue
         iatflg(iat) = 2
         itst = 1
         return
 40   continue
      return
      end
