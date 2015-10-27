*deck bsplpp
      subroutine bsplpp ( t, bcoef, n, k, scrtch, break, coef, l )
converts b-representation to pp-representation.
c     dimension t(n+k), break(*+1), coef(k,*)
c  here, * = the final value of the output parameter l .
      implicit real *8 (a-h,o-z)
      save
      dimension t(1),bcoef(n), scrtch(k,k), break(1),coef(k,1)
      dimension biatx(20)
      l = 0
      break(1) = t(k)
      do 50 ileft=k,n
c        find the next nontrivial knot interval.
         if (t(ileft+1) .eq. t(ileft))  go to 50
         l = l + 1
         break(l+1) = t(ileft+1)
         if (k .gt. 1)                  go to 9
         coef(1,l) = bcoef(ileft)
                                        go to 50
c        store the k b-spline coeff.s relevant to current knot interval
c        in  scrtch(.,1) .
    9    do 10 i=1,k
   10       scrtch(i,1) = bcoef(ileft-k+i)
c        for j=1,...,k-1, compute the k-j b-spline coeff.s relevant to
c        current knot interval for the j-th derivative by differencing
c        those for the (j-1)st derivative, and store in scrtch(.,j+1) .
         do 20 jp1=2,k
            j = jp1 - 1
            kmj = k - j
            fkmj = float(kmj)
            do 20 i=1,kmj
               diff = t(ileft+i) - t(ileft+i - kmj)
               if (diff .gt. 0.d0)  scrtch(i,jp1) =
     *                       ((scrtch(i+1,j)-scrtch(i,j))/diff)*fkmj
   20          continue
c        starting with the one b-spline of order 1 not zero at t(ileft),
c        find the values at t(ileft) of the j+1 b-splines of order j+1
c        not identically zero there from those of order j, then combine
c        with the b-spline coeff.s found earlier to compute the (k-j)-
c        th derivative at t(ileft) of the given spline.
         call bsplvb ( t, 1, 1, t(ileft), ileft, biatx )
         coef(k,l) = scrtch(1,k)
         do 30 jp1=2,k
            call bsplvb ( t, jp1, 2, t(ileft), ileft, biatx )
            kmj = k+1 - jp1
            sum = 0.d0
            do 28 i=1,jp1
   28          sum = biatx(i)*scrtch(i,kmj) + sum
   30       coef(kmj,l) = sum
   50    continue
                                        return
      end
