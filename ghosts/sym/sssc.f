*deck sssc
      subroutine sssc(maxap3,natoms,iatflg,a,ixyz,itst)
      implicit real*8(a-h,o-z)
c
c     flag those atoms which are on cartesian axis ixyz.
c
      common/tol/ toler,tol2
      dimension iatflg(1), a(maxap3,3)
c
c     call rtrace(6hsssc  ,1)
      i1 = 1 + mod(ixyz,3)
      i2 = 1 + mod(i1  ,3)
      itst = 0
      do 20 iat=1,natoms
         if (abs(a(iat,i1)) .gt. toler  .or.
     $        abs(a(iat,i2)) .gt. toler  .or.
     $        iatflg(iat) .ne. 0)  goto 20
         itst = 1
         iatflg(iat) = 2
 20   continue
      return
      end
