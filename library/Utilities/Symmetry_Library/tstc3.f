*deck @(#)tstc3.f	5.1  11/6/94
      subroutine tstc3(maxap3,a,b,natoms,atmchg,iat,jat,kat,
     $     centr,itst)
      implicit real*8(a-h,o-z)
c
c     are the three atoms iat, jat, and kat interchangeable via a
c     3-fold rotation?
c     if no  itst=0, return
c     if yes itst=1, align c3 with z, return
c     centr is the point of intersection of the c3 axis with the
c     plane defined by the 3 atoms.
c
      common/tol/ toler,tol2
      dimension a(maxap3,3), b(maxap3,3), atmchg(1), centr(3), t(3,3)
      data half,one,two,three/0.5d+00,1.0d+00,2.0d+00,3.0d+00/
      data eight/8.0d+00/
      save half,one,two,three,eight
c
c     call rtrace(6htstc3 ,1)
      numatm = natoms + 3
      itst   = 0
      phi3   = (eight/three) * atan(one)
      theta3 = half * phi3
      fact3  = two / three
c
c     get the angles and sides of the triangle defined by the three
c     atoms.
c
      call triang(maxap3,a,iat,jat,kat,alpha,beta,gamma,dij,dik,djk)
c
c     do the three points form an equilateral triangle?  it is only
c     necessary to check to see if one angle is 60 degrees and that
c     two sides are equal.
c
      if (abs(alpha-theta3) .gt. tol2  .or.
     $     abs(dij-dik) .gt. toler) return
      px = half * (a(jat,1)+a(kat,1))
      py = half * (a(jat,2)+a(kat,2))
      pz = half * (a(jat,3)+a(kat,3))
      centr(1) = a(iat,1) + fact3 * (px-a(iat,1))
      centr(2) = a(iat,2) + fact3 * (py-a(iat,2))
      centr(3) = a(iat,3) + fact3 * (pz-a(iat,3))
      call putsym(maxap3,a,b,t,centr,numatm,3)
      call rotate(maxap3,a,b,natoms,t,3,phi3)
      call equiv(maxap3,a,b,atmchg,natoms,itst)
      return
      end
