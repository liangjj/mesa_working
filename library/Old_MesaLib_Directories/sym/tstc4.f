*deck @(#)tstc4.f	5.1  11/6/94
      subroutine tstc4(maxap3,a,b,natoms,atmchg,iat,jat,kat,
     $     centr,itst)
      implicit real*8(a-h,o-z)
c
c     are the three atoms iat, jat, and kat interchangeable via a
c     4-fold rotation?
c     if no  itst=0, return
c     if yes itst=1, align c4 with z, return
c     centr is the point of intersection of the c4 axis with the
c     plane defined by the 3 atoms.
c
      common/tol/ toler,tol2
      dimension a(maxap3,3), b(maxap3,3), atmchg(1), centr(3), t(3,3)
      data half,one,two/0.5d+00,1.0d+00,2.0d+00/
      save half,one,two
c
c     call rtrace(6htstc4 ,1)
      numatm = natoms + 3
      itst   = 0
      halfpi = two * atan(one)
c
c     get the angles and sides of the triangle defined by the three
c     atoms.
c
      call triang(maxap3,a,iat,jat,kat,alpha,beta,gamma,dij,dik,djk)
c
c     are any of these angles equal to 90 degrees (i.e. the internal
c     angle of a square) and thus possibly equivalent by a
c     4-fold axis of symmetry?
c     are two of the sides of the triangle of equal length?
c
      if (abs(alpha-halfpi) .gt. tol2  .or.
     $     abs(dij-dik)      .gt. toler) goto 20
      centr(1) = half * (a(jat,1)+a(kat,1))
      centr(2) = half * (a(jat,2)+a(kat,2))
      centr(3) = half * (a(jat,3)+a(kat,3))
      goto 60
c
 20   continue
      if (abs(beta-halfpi) .gt. tol2  .or.
     $     abs(dij-djk)     .gt. toler) goto 40
      centr(1) = half * (a(iat,1)+a(kat,1))
      centr(2) = half * (a(iat,2)+a(kat,2))
      centr(3) = half * (a(iat,3)+a(kat,3))
      goto 60
c
 40   continue
      if (abs(gamma-halfpi) .gt. tol2  .or.
     $     abs(dik-djk)      .gt. toler) return
      centr(1) = half * (a(iat,1)+a(jat,1))
      centr(2) = half * (a(iat,2)+a(jat,2))
      centr(3) = half * (a(iat,3)+a(jat,3))
c
 60   continue
      call putsym(maxap3,a,b,t,centr,numatm,3)
      call rotate(maxap3,a,b,natoms,t,3,halfpi)
      call equiv(maxap3,a,b,atmchg,natoms,itst)
      return
      end
