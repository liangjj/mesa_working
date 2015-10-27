*deck cdcor
      subroutine cdcor (dfdy, el, fa, h, ierror, impl, ipvt, matdim,
     8   miter, ml, mu, n, nde, nq, t, users, y, yh, ywt, evalfa, save1,
     8   save2, a, d, jstate)
c***begin prologue  cdcor
c***subsidiary
c***purpose  subroutine cdcor computes corrections to the y array.
c***library   slatec (sdrive)
c***type      complex (sdcor-s, ddcor-d, cdcor-c)
c***author  kahaner, d. k., (nist)
c             national institute of standards and technology
c             gaithersburg, md  20899
c           sutherland, c. d., (lanl)
c             mail stop d466
c             los alamos national laboratory
c             los alamos, nm  87545
c***description
c
c  in the case of functional iteration, update y directly from the
c  result of the last call to f.
c  in the case of the chord method, compute the corrector error and
c  solve the linear system with that as right hand side and dfdy as
c  coefficient matrix, using the lu decomposition if miter is 1, 2, 4,
c  or 5.
c
c***routines called  cgbsl, cgesl, scnrm2
c***revision history  (yymmdd)
c   790601  date written
c   900329  initial submission to slatec.
c***end prologue  cdcor
      integer i, ierror, iflag, impl, j, jstate, matdim, miter, ml, mu,
     8        mw, n, nde, nq
      complex a(matdim,*), dfdy(matdim,*), save1(*), save2(*), y(*),
     8        yh(n,*), ywt(*)
      real d, el(13,12), h, scnrm2, t
      integer ipvt(*)
      logical evalfa
c***first executable statement  cdcor
      if (miter .eq. 0) then
        if (ierror .eq. 1 .or. ierror .eq. 5) then
          do 100 i = 1,n
 100        save1(i) = (h*save2(i) - yh(i,2) - save1(i))/ywt(i)
        else
          do 102 i = 1,n
            save1(i) = (h*save2(i) - yh(i,2) - save1(i))/
     8      max(abs(y(i)), abs(ywt(i)))
 102        continue
        end if
        d = scnrm2(n, save1, 1)/sqrt(real(n))
        do 105 i = 1,n
 105      save1(i) = h*save2(i) - yh(i,2)
      else if (miter .eq. 1 .or. miter .eq. 2) then
        if (impl .eq. 0) then
          do 130 i = 1,n
 130        save2(i) = h*save2(i) - yh(i,2) - save1(i)
        else if (impl .eq. 1) then
          if (evalfa) then
            call fa (n, t, y, a, matdim, ml, mu, nde)
            if (n .eq. 0) then
              jstate = 9
              return
            end if
          else
            evalfa = .true.
          end if
          do 150 i = 1,n
 150        save2(i) = h*save2(i)
          do 160 j = 1,n
            do 160 i = 1,n
 160          save2(i) = save2(i) - a(i,j)*(yh(j,2) + save1(j))
        else if (impl .eq. 2) then
          if (evalfa) then
            call fa (n, t, y, a, matdim, ml, mu, nde)
            if (n .eq. 0) then
              jstate = 9
              return
            end if
          else
            evalfa = .true.
          end if
          do 180 i = 1,n
 180        save2(i) = h*save2(i) - a(i,1)*(yh(i,2) + save1(i))
        else if (impl .eq. 3) then
          if (evalfa) then
            call fa (n, t, y, a, matdim, ml, mu, nde)
            if (n .eq. 0) then
              jstate = 9
              return
            end if
          else
            evalfa = .true.
          end if
          do 140 i = 1,n
 140        save2(i) = h*save2(i)
          do 170 j = 1,nde
            do 170 i = 1,nde
 170          save2(i) = save2(i) - a(i,j)*(yh(j,2) + save1(j))
        end if
        call cgesl (dfdy, matdim, n, ipvt, save2, 0)
        if (ierror .eq. 1 .or. ierror .eq. 5) then
          do 200 i = 1,n
            save1(i) = save1(i) + save2(i)
 200        save2(i) = save2(i)/ywt(i)
        else
          do 205 i = 1,n
            save1(i) = save1(i) + save2(i)
 205        save2(i) = save2(i)/max(abs(y(i)), abs(ywt(i)))
        end if
        d = scnrm2(n, save2, 1)/sqrt(real(n))
      else if (miter .eq. 4 .or. miter .eq. 5) then
        if (impl .eq. 0) then
          do 230 i = 1,n
 230        save2(i) = h*save2(i) - yh(i,2) - save1(i)
        else if (impl .eq. 1) then
          if (evalfa) then
            call fa (n, t, y, a(ml+1,1), matdim, ml, mu, nde)
            if (n .eq. 0) then
              jstate = 9
              return
            end if
          else
            evalfa = .true.
          end if
          do 250 i = 1,n
 250        save2(i) = h*save2(i)
          mw = ml + 1 + mu
          do 260 j = 1,n
            do 260 i = max(ml+1, mw+1-j), min(mw+n-j, mw+ml)
              save2(i+j-mw) = save2(i+j-mw)
     8                        - a(i,j)*(yh(j,2) + save1(j))
 260        continue
        else if (impl .eq. 2) then
          if (evalfa) then
            call fa (n, t, y, a, matdim, ml, mu, nde)
            if (n .eq. 0) then
              jstate = 9
              return
            end if
          else
            evalfa = .true.
          end if
          do 280 i = 1,n
 280        save2(i) = h*save2(i) - a(i,1)*(yh(i,2) + save1(i))
        else if (impl .eq. 3) then
          if (evalfa) then
            call fa (n, t, y, a(ml+1,1), matdim, ml, mu, nde)
            if (n .eq. 0) then
              jstate = 9
              return
            end if
          else
            evalfa = .true.
          end if
          do 270 i = 1,n
 270        save2(i) = h*save2(i)
          mw = ml + 1 + mu
          do 290 j = 1,nde
            do 290 i = max(ml+1, mw+1-j), min(mw+nde-j, mw+ml)
              save2(i+j-mw) = save2(i+j-mw)
     8                        - a(i,j)*(yh(j,2) + save1(j))
 290        continue
        end if
        call cgbsl (dfdy, matdim, n, ml, mu, ipvt, save2, 0)
        if (ierror .eq. 1 .or. ierror .eq. 5) then
          do 300 i = 1,n
            save1(i) = save1(i) + save2(i)
 300        save2(i) = save2(i)/ywt(i)
        else
          do 305 i = 1,n
            save1(i) = save1(i) + save2(i)
 305        save2(i) = save2(i)/max(abs(y(i)), abs(ywt(i)))
        end if
        d = scnrm2(n, save2, 1)/sqrt(real(n))
      else if (miter .eq. 3) then
        iflag = 2
        call users (y, yh(1,2), ywt, save1, save2, t, h, el(1,nq), impl,
     8              n, nde, iflag)
        if (n .eq. 0) then
          jstate = 10
          return
        end if
        if (ierror .eq. 1 .or. ierror .eq. 5) then
          do 320 i = 1,n
            save1(i) = save1(i) + save2(i)
 320        save2(i) = save2(i)/ywt(i)
        else
          do 325 i = 1,n
            save1(i) = save1(i) + save2(i)
 325        save2(i) = save2(i)/max(abs(y(i)), abs(ywt(i)))
        end if
        d = scnrm2(n, save2, 1)/sqrt(real(n))
      end if
      return
      end
