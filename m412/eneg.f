*deck @(#)eneg.f	5.1  11/6/94
      function eneg(ia,j)
c***begin prologue     eneg
c***date written       850601  yymmdd
c***revision date      yymmdd  yymmdd
c***keywords           huckel, energies
c***author             martin, richard (lanl)
c***source             @(#)eneg.f	5.1   11/6/94
c***purpose            returns the diagonal energy associated with a
c                      specific atom/orbital.
c***description
c     eneg is a real function used as:
c            diagnl=eneg(ia,j)
c            ia     atomic number.
c            j      orbital type.
c
c
c     the values for elements 1-18 are from basch,viste,and gray.
c     values for elements 19-110 are from desclaux's numerical dirac-fock
c     work.
c
c***references         basch,viste,and gray,theor.chim.acta,vol.3,458(1965).
c                      desclaux,atomic data and nuclear data tables,
c                      vol.12,311(1973).
c***routines called    (none)
c***end prologue       eneg
c
      implicit real*8(a-h,o-z)
c
      real*8 eneg
c
      dimension val(110,4), vals(110), valp(110), vald(110), valf(110),
     $          itype(19)
      equivalence (val(1,1),vals(1)), (val(1,2),valp(1)),
     $            (val(1,3),vald(1)), (val(1,4),valf(1))
      data vals/
     $         0.50d0,0.90d0,0.20d0,0.34d0,0.52d0,0.71d0,0.94d0,
     $         1.19d0,1.47d0,1.78d0,0.19d0,0.28d0,0.41d0,0.54d0,0.69d0,
     $         0.76d0,0.93d0,1.07d0,0.15d0,0.20d0,0.21d0,0.22d0,0.23d0,
     $         0.21d0,0.25d0,0.26d0,0.27d0,0.28d0,0.24d0,0.30d0,
     $         0.43d0,0.57d0,0.71d0,0.86d0,1.02d0,1.19d0,0.14d0,
     $         0.18d0,0.20d0,0.22d0,0.21d0,0.21d0,0.25d0,0.22d0,0.23d0,
     $         0.23d0,0.23d0,0.28d0,0.40d0,0.51d0,0.63d0,0.75d0,0.88d0,
     $         1.01d0,0.13d0,0.16d0,0.18d0,0.18d0,0.17d0,0.17d0,0.18d0,
     $         0.18d0,0.18d0,0.20d0,0.19d0,0.19d0,0.19d0,0.19d0,0.19d0,
     $         0.20d0,0.22d0,0.24d0,0.25d0,0.27d0,0.28d0,0.29d0,0.30d0,
     $         0.28d0,0.30d0,0.33d0,0.45d0,0.57d0,0.69d0,0.81d0,0.94d0,
     $         1.07d0,0.13d0,0.17d0,0.19d0,0.21d0,0.20d0,0.20d0,0.21d0,
     $         0.19d0,0.19d0,0.22d0,0.20d0,0.20d0,0.20d0,0.20d0,0.21d0,
     $         0.21d0,0.25d0,0.27d0,0.30d0,0.32d0,0.34d0,0.36d0,0.38d0,
     $         0.41d0/
      data valp/
     $         0.0d0,0.0d0,0.13d0,0.23d0,0.30d0,0.39d0,0.48d0,0.58d0,
     $         0.69d0,0.79d0,0.11d0,0.16d0,0.22d0,0.28d0,0.37d0,
     $         0.43d0,0.50d0,0.58d0,0.07d0,0.10d0,0.10d0,0.11d0,0.11d0,
     $         0.10d0,0.12d0,0.13d0,0.13d0,0.14d0,0.12d0,0.15d0,0.21d0,
     $         0.27d0,0.33d0,0.39d0,0.46d0,0.52d0,0.07d0,0.09d0,0.10d0,
     $         0.11d0,0.10d0,0.11d0,0.12d0,0.11d0,0.11d0,0.11d0,0.12d0,
     $         0.14d0,0.19d0,0.24d0,0.30d0,0.35d0,0.40d0,0.45d0,0.07d0,
     $         0.08d0,0.09d0,0.09d0,0.09d0,0.09d0,0.09d0,0.09d0,0.09d0,
     $         0.10d0,0.09d0,0.09d0,0.09d0,0.09d0,0.09d0,0.10d0,0.11d0,
     $         0.12d0,0.13d0,0.13d0,0.14d0,0.14d0,0.15d0,0.14d0,0.14d0,
     $         0.16d0,0.20d0,0.26d0,0.32d0,0.37d0,0.44d0,0.50d0,0.06d0,
     $         0.08d0,0.09d0,0.10d0,0.10d0,0.10d0,0.10d0,0.10d0,0.10d0,
     $         0.11d0,0.10d0,0.10d0,0.10d0,0.10d0,0.10d0,0.10d0,0.12d0,
     $         0.14d0,0.15d0,0.15d0,0.17d0,0.18d0,0.19d0,0.20d0/
      data vald/
     $         20*0.0d0,
     $         0.34d0,0.40d0,0.45d0,0.32d0,0.55d0,0.60d0,0.64d0,0.69d0,
     $         0.49d0,0.77d0,1.17d0,1.62d0,2.10d0,2.62d0,3.18d0,3.78d0,
     $         0.00d0,0.00d0,0.23d0,0.29d0,0.26d0,0.30d0,0.46d0,0.39d0,
     $         0.43d0,0.34d0,0.53d0,0.74d0,1.03d0,1.34d0,1.66d0,1.99d0,
     $         2.34d0,2.71d0,2*0.0d0,14*0.24d0,0.19d0,0.25d0,0.3d0,
     $         0.35d0,0.40d0,0.45d0,0.50d0,0.45d0,0.49d0,0.65d0,0.89d0,
     $         1.14d0,1.39d0,1.65d0,1.91d0,2.19d0,0.00d0,0.00d0,0.18d0,
     $         0.22d0,14*0.19d0,0.25d0,0.30d0,0.34d0,0.38d0,0.43d0,
     $         0.47d0/
      data valf/
     $         56*0.0d0,0.12d0,
     $         0.54d0,0.33d0,0.37d0,0.40d0,0.42d0,0.45d0,0.74d0,0.48d0,
     $         0.50d0,0.51d0,0.52d0,0.53d0,0.54d0,0.86d0,1.18d0,1.53d0,
     $         1.89d0,2.28d0,2.68d0,3.10d0,3.41d0,3.87d0,4.47d0,5.19d0,
     $         5.93d0,6.70d0,7.50d0,8.34d0,9.19d0,0.00d0,0.00d0,0.10d0,
     $         0.10d0,0.29d0,0.35d0,0.40d0,0.30d0,0.34d0,0.54d0,0.41d0,
     $         0.44d0,0.47d0,0.51d0,0.54d0,0.57d0,0.82d0,1.07d0,1.32d0,
     $         1.58d0,1.85d0,2.13d0,2.42d0,2.72d0/
      data itype/1,1,2,1,2,1,2,3,1,2,3,1,2,3,4,1,2,3,4/
      data zero/0.0d0/, two/2.0d0/, pt32/0.32d0/, pt0045/0.0045d0/
      save vals,valp,vald,valf,itype,zero,two,pt32,pt0045
c
c     valence values for 110 elements spdf.  itype maps j onto spdf.
c
      eneg = zero
      if(j.gt.0) eneg = -val(ia,itype(j))
      if(j.eq.-1) eneg = -pt32 * float(ia)**2
      if(j.le.-2) eneg = -pt0045 * float(ia)**3 / float(-j-1)
      if(ia.gt.17) eneg = eneg * two
      return
      end
