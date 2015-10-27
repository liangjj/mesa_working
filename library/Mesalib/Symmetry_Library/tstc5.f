*deck @(#)tstc5.f	5.1  11/6/94
      subroutine tstc5(maxap3,a,b,natoms,atmchg,iat,jat,kat,
     $     centr,itst)
      implicit real*8(a-h,o-z)
c
c     are the three atoms iat, jat, and kat interchangeable via a
c     5-fold rotation?
c     if no  itst=0, return
c     if yes itst=1, align c5 with z, return
c     centr is the point of intersection of the c5 axis with the
c     plane defined by the 3 atoms.
c
      common/tol/ toler,tol2
      dimension a(maxap3,3), b(maxap3,3), atmchg(1), centr(3), t(3,3)
      data half,pt8,one,onept6/0.5d+00,0.8d+00,1.0d+00,1.6d+00/,
     $     two,twopt4/2.0d+00,2.4d+00/
      save half,pt8,one,onept6,two,twopt4
c
c     call rtrace(6htstc5 ,1)
      numatm = natoms + 3
      itst   = 0
      piovr4 = atan(one)
      phi5   = onept6 * piovr4
      theta5 = twopt4 * piovr4
      fact5  = one / ( two * sin(pt8*piovr4)**2 )
c
c     get the angles and sides of the triangle defined by the three
c     atoms.
c
      call triang(maxap3,a,iat,jat,kat,alpha,beta,gamma,dij,dik,djk)
c
c     are any of these angles equal to 108 degrees (i.e. the internal
c     angle of a regular pentagon) and thus possibly equivalent by a
c     5-fold axis of symmetry?
c     are two of the sides of the triangle of equal length?
c
      if (abs(alpha-theta5) .gt. tol2  .or.
     $     abs(dij-dik)      .gt. toler) goto 20
      px = half * (a(jat,1)+a(kat,1))
      py = half * (a(jat,2)+a(kat,2))
      pz = half * (a(jat,3)+a(kat,3))
      centr(1) = a(iat,1) + fact5 * (px-a(iat,1))
      centr(2) = a(iat,2) + fact5 * (py-a(iat,2))
      centr(3) = a(iat,3) + fact5 * (pz-a(iat,3))
      goto 60
c
 20   continue
      if (abs(beta-theta5) .gt. tol2  .or.
     $     abs(dij-djk)     .gt. toler) goto 40
      px = half * (a(iat,1)+a(kat,1))
      py = half * (a(iat,2)+a(kat,2))
      pz = half * (a(iat,3)+a(kat,3))
      centr(1) = a(jat,1) + fact5 * (px-a(jat,1))
      centr(2) = a(jat,2) + fact5 * (py-a(jat,2))
      centr(3) = a(jat,3) + fact5 * (pz-a(jat,3))
      goto 60
c
 40   continue
      if (abs(gamma-theta5) .gt. tol2  .or.
     $     abs(dik-djk)      .gt. toler) return
      px = half * (a(iat,1)+a(jat,1))
      py = half * (a(iat,2)+a(jat,2))
      pz = half * (a(iat,3)+a(jat,3))
      centr(1) = a(kat,1) + fact5 * (px-a(kat,1))
      centr(2) = a(kat,2) + fact5 * (py-a(kat,2))
      centr(3) = a(kat,3) + fact5 * (pz-a(kat,3))
c
 60   continue
      call putsym(maxap3,a,b,t,centr,numatm,3)
      call rotate(maxap3,a,b,natoms,t,3,phi5)
      call equiv(maxap3,a,b,atmchg,natoms,itst)
      return
      end
