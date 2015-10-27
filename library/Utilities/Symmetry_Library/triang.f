*deck @(#)triang.f	5.1  11/6/94
      subroutine triang(maxap3,a,i,j,k,alp,bet,gam,dij,dik,djk)
      implicit real*8(a-h,o-z)
c
c     given the atomic coordinates in a and the three atoms i, j, and k
c     find:
c     1-- the distances defined by the locations of the three atoms,
c     dij, dik, and djk.
c     2-- the angles defined by the location of the three atoms,
c     alp    k.i.j
c     bet    i.j.k
c     gam    j.k.i
c
      dimension a(maxap3,3)
c
c     call rtrace(6htriang,1)
      xi  = a(i,1)
      yi  = a(i,2)
      zi  = a(i,3)
      xj  = a(j,1)
      yj  = a(j,2)
      zj  = a(j,3)
      xk  = a(k,1)
      yk  = a(k,2)
      zk  = a(k,3)
      dij    = sqrt( (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2 )
      dik    = sqrt( (xi-xk)**2 + (yi-yk)**2 + (zi-zk)**2 )
      djk    = sqrt( (xj-xk)**2 + (yj-yk)**2 + (zj-zk)**2 )
      dotjk  = (xj-xi)*(xk-xi) + (yj-yi)*(yk-yi) + (zj-zi)*(zk-zi)
      dotik  = (xi-xj)*(xk-xj) + (yi-yj)*(yk-yj) + (zi-zj)*(zk-zj)
      dotij  = (xi-xk)*(xj-xk) + (yi-yk)*(yj-yk) + (zi-zk)*(zj-zk)
      alp    = acos( dotjk / (dij*dik) )
      bet    = acos( dotik / (dij*djk) )
      gam    = acos( dotij / (dik*djk) )
      return
      end
