*deck cdpst
      subroutine cdpst (el, f, fa, h, impl, jacobn, matdim, miter, ml,
     8   mu, n, nde, nq, save2, t, users, y, yh, ywt, uround, nfe, nje,
     8   a, dfdy, fac, ier, ipvt, save1, iswflg, bnd, jstate)
c***begin prologue  cdpst
c***subsidiary
c***purpose  subroutine cdpst evaluates the jacobian matrix of the right
c            hand side of the differential equations.
c***library   slatec (sdrive)
c***type      complex (sdpst-s, ddpst-d, cdpst-c)
c***author  kahaner, d. k., (nist)
c             national institute of standards and technology
c             gaithersburg, md  20899
c           sutherland, c. d., (lanl)
c             mail stop d466
c             los alamos national laboratory
c             los alamos, nm  87545
c***description
c
c  if miter is 1, 2, 4, or 5, the matrix
c  p = i - l(0)*h*jacobian is stored in dfdy and subjected to lu
c  decomposition, with the results also stored in dfdy.
c
c***routines called  cgbfa, cgefa, scnrm2
c***revision history  (yymmdd)
c   790601  date written
c   900329  initial submission to slatec.
c***end prologue  cdpst
      integer i, iflag, imax, impl, info, iswflg, j, j2, jstate, k,
     8        matdim, miter, ml, mu, mw, n, nde, nfe, nje, nq
      complex a(matdim,*), cfctr, dfdy(matdim,*), dy, fac(*), save1(*),
     8        save2(*), y(*), yh(n,*), yj, ys, ywt(*)
      real bl, bnd, bp, br, bu, dfdymx, diff, el(13,12), facmax, facmin,
     8     factor, h, scale, scnrm2, t, uround, zmax, zmin
      integer ipvt(*)
      logical ier
      parameter(facmax = .5e0, bu = 0.5e0)
c***first executable statement  cdpst
      nje = nje + 1
      ier = .false.
      if (miter .eq. 1 .or. miter .eq. 2) then
        if (miter .eq. 1) then
          call jacobn (n, t, y, dfdy, matdim, ml, mu)
          if (n .eq. 0) then
            jstate = 8
            return
          end if
          if (iswflg .eq. 3) bnd = scnrm2(n*n, dfdy, 1)
          factor = -el(1,nq)*h
          do 110 j = 1,n
            do 110 i = 1,n
 110          dfdy(i,j) = factor*dfdy(i,j)
        else if (miter .eq. 2) then
          br = uround**(.875e0)
          bl = uround**(.75e0)
          bp = uround**(-.15e0)
          facmin = uround**(.78e0)
          do 170 j = 1,n
            if (abs(y(j)) .gt. abs(ywt(j))) then
              ys = y(j)
            else
              ys = ywt(j)
            end if
 120        dy = fac(j)*ys
            if (dy .eq. 0.e0) then
              if (real(fac(j)) .lt. facmax) then
                fac(j) = min(100.e0*real(fac(j)), facmax)
                go to 120
              else
                dy = ys
              end if
            end if
            dy = (y(j) + dy) - y(j)
            yj = y(j)
            y(j) = y(j) + dy
            call f (n, t, y, save1)
            if (n .eq. 0) then
              jstate = 6
              return
            end if
            y(j) = yj
            cfctr = -el(1,nq)*h/dy
            do 140 i = 1,n
 140          dfdy(i,j) = (save1(i) - save2(i))*cfctr
c                                                                 step 1
            diff = abs(save2(1) - save1(1))
            imax = 1
            do 150 i = 2,n
              if (abs(save2(i) - save1(i)) .gt. diff) then
                imax = i
                diff = abs(save2(i) - save1(i))
              end if
 150          continue
c                                                                 step 2
            if (min(abs(save2(imax)), abs(save1(imax))) .gt. 0.e0) then
              scale = max(abs(save2(imax)), abs(save1(imax)))
c                                                                 step 3
              if (diff .gt. bu*scale) then
                fac(j) = max(facmin, real(fac(j))*.5e0)
              else if (br*scale .le. diff .and. diff .le. bl*scale) then
                fac(j) = min(real(fac(j))*2.e0, facmax)
c                                                                 step 4
              else if (diff .lt. br*scale) then
                fac(j) = min(bp*real(fac(j)), facmax)
              end if
            end if
 170        continue
          if (iswflg .eq. 3) bnd = scnrm2(n*n, dfdy, 1)/(-el(1,nq)*h)
          nfe = nfe + n
        end if
        if (impl .eq. 0) then
          do 190 i = 1,n
 190        dfdy(i,i) = dfdy(i,i) + 1.e0
        else if (impl .eq. 1) then
          call fa (n, t, y, a, matdim, ml, mu, nde)
          if (n .eq. 0) then
            jstate = 9
            return
          end if
          do 210 j = 1,n
            do 210 i = 1,n
 210          dfdy(i,j) = dfdy(i,j) + a(i,j)
        else if (impl .eq. 2) then
          call fa (n, t, y, a, matdim, ml, mu, nde)
          if (n .eq. 0) then
            jstate = 9
            return
          end if
          do 230 i = 1,nde
 230        dfdy(i,i) = dfdy(i,i) + a(i,1)
        else if (impl .eq. 3) then
          call fa (n, t, y, a, matdim, ml, mu, nde)
          if (n .eq. 0) then
            jstate = 9
            return
          end if
          do 220 j = 1,nde
            do 220 i = 1,nde
 220          dfdy(i,j) = dfdy(i,j) + a(i,j)
        end if
        call cgefa (dfdy, matdim, n, ipvt, info)
        if (info .ne. 0) ier = .true.
      else if (miter .eq. 4 .or. miter .eq. 5) then
        if (miter .eq. 4) then
          call jacobn (n, t, y, dfdy(ml+1,1), matdim, ml, mu)
          if (n .eq. 0) then
            jstate = 8
            return
          end if
          factor = -el(1,nq)*h
          mw = ml + mu + 1
          do 260 j = 1,n
            do 260 i = max(ml+1, mw+1-j), min(mw+n-j, mw+ml)
 260          dfdy(i,j) = factor*dfdy(i,j)
        else if (miter .eq. 5) then
          br = uround**(.875e0)
          bl = uround**(.75e0)
          bp = uround**(-.15e0)
          facmin = uround**(.78e0)
          mw = ml + mu + 1
          j2 = min(mw, n)
          do 340 j = 1,j2
            do 290 k = j,n,mw
              if (abs(y(k)) .gt. abs(ywt(k))) then
                ys = y(k)
              else
                ys = ywt(k)
              end if
 280          dy = fac(k)*ys
              if (dy .eq. 0.e0) then
                if (real(fac(k)) .lt. facmax) then
                  fac(k) = min(100.e0*real(fac(k)), facmax)
                  go to 280
                else
                  dy = ys
                end if
              end if
              dy = (y(k) + dy) - y(k)
              dfdy(mw,k) = y(k)
 290          y(k) = y(k) + dy
            call f (n, t, y, save1)
            if (n .eq. 0) then
              jstate = 6
              return
            end if
            do 330 k = j,n,mw
              dy = y(k) - dfdy(mw,k)
              y(k) = dfdy(mw,k)
              cfctr = -el(1,nq)*h/dy
              do 300 i = max(ml+1, mw+1-k), min(mw+n-k, mw+ml)
 300            dfdy(i,k) = cfctr*(save1(i+k-mw) - save2(i+k-mw))
c                                                                 step 1
              imax = max(1, k - mu)
              diff = abs(save2(imax) - save1(imax))
              do 310 i = max(1, k - mu)+1, min(k + ml, n)
                if (abs(save2(i) - save1(i)) .gt. diff) then
                  imax = i
                  diff = abs(save2(i) - save1(i))
                end if
 310            continue
c                                                                 step 2
              if (min(abs(save2(imax)), abs(save1(imax))) .gt.0.e0) then
                scale = max(abs(save2(imax)), abs(save1(imax)))
c                                                                 step 3
                if (diff .gt. bu*scale) then
                  fac(j) = max(facmin, real(fac(j))*.5e0)
                else if (br*scale .le.diff .and. diff .le.bl*scale) then
                  fac(j) = min(real(fac(j))*2.e0, facmax)
c                                                                 step 4
                else if (diff .lt. br*scale) then
                  fac(k) = min(bp*real(fac(k)), facmax)
                end if
              end if
 330          continue
 340        continue
          nfe = nfe + j2
        end if
        if (iswflg .eq. 3) then
          dfdymx = 0.e0
          do 345 j = 1,n
            do 345 i = max(ml+1, mw+1-j), min(mw+n-j, mw+ml)
              zmax = max(abs(real(dfdy(i,j))), abs(aimag(dfdy(i,j))))
              zmin = min(abs(real(dfdy(i,j))), abs(aimag(dfdy(i,j))))
              if (zmax .ne. 0.e0)
     8        dfdymx = max(dfdymx, zmax*sqrt(1.e0+ (zmin/zmax)**2))
 345          continue
          bnd = 0.e0
          if (dfdymx .ne. 0.e0) then
            do 350 j = 1,n
              do 350 i = max(ml+1, mw+1-j), min(mw+n-j, mw+ml)
                bnd = bnd + (real(dfdy(i,j))/dfdymx)**2 +
     8          (aimag(dfdy(i,j))/dfdymx)**2
 350            continue
            bnd = dfdymx*sqrt(bnd)/(-el(1,nq)*h)
          end if
        end if
        if (impl .eq. 0) then
          do 360 j = 1,n
 360        dfdy(mw,j) = dfdy(mw,j) + 1.e0
        else if (impl .eq. 1) then
          call fa (n, t, y, a(ml+1,1), matdim, ml, mu, nde)
          if (n .eq. 0) then
            jstate = 9
            return
          end if
          do 380 j = 1,n
            do 380 i = max(ml+1, mw+1-j), min(mw+n-j, mw+ml)
 380          dfdy(i,j) = dfdy(i,j) + a(i,j)
        else if (impl .eq. 2) then
          call fa (n, t, y, a, matdim, ml, mu, nde)
          if (n .eq. 0) then
            jstate = 9
            return
          end if
          do 400 j = 1,nde
 400        dfdy(mw,j) =  dfdy(mw,j) + a(j,1)
        else if (impl .eq. 3) then
          call fa (n, t, y, a(ml+1,1), matdim, ml, mu, nde)
          if (n .eq. 0) then
            jstate = 9
            return
          end if
          do 390 j = 1,nde
            do 390 i = max(ml+1, mw+1-j), min(mw+nde-j, mw+ml)
 390          dfdy(i,j) = dfdy(i,j) + a(i,j)
        end if
        call cgbfa (dfdy, matdim, n, ml, mu, ipvt, info)
        if (info .ne. 0) ier = .true.
      else if (miter .eq. 3) then
        iflag = 1
        call users (y, yh(1,2), ywt, save1, save2, t, h, el(1,nq), impl,
     8              n, nde, iflag)
        if (iflag .eq. -1) then
          ier = .true.
          return
        end if
        if (n .eq. 0) then
          jstate = 10
          return
        end if
      end if
      return
      end
