*               g r a n g e
*     ------------------------------------------------------------------
*
*       controls the calculation of off-diagonal energy parameters.
*   it searches for all pairs (i,j) which are constrained through  an
*   orthogonality requirement.   eq. (7-10) is used to calculate the
*   PARAMETER.
*
      SUBROUTINE grange
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      PARAMETER(nt=100,kx=8,nwd=10)
      COMMON /coeff/coef(100),ijptr(5,5)
      CHARACTER config*50,el*3,atom*6,term*6
      COMMON /LABEL/config,el(nwd),atom,term
*
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
      COMMON /ORBSCF/sum(nwd), s(nwd), dpm(nwd), acc(nwd), meth(nwd),
     :               iord(nwd), e(nwd,nwd)
*
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
      LOGICAL fail,omit,rel,all,trace
      COMMON /TEST/fail,omit,rel,all,trace
*
*  *****  rotate pairs connected by orthogonality but not when one of
*         the orbitals is simultaneously orthogonal to a non-orthogonal
*         pair
*
      do 1 i = ib,nwf-1
         do 2 j = i+1,nwf
            if (dabs(e(i,j)) .gt. 1.d-10) call rotate(i,j)
2        continue
1     continue
      print *, 'sumi, sumj', sum(i), sum(j)
*
*  *****   compute off-diagonal energy parameters
*
      do 10 i = max0(2,ib),nwf
         do 12 j = 1,i-1
           if (dabs(e(i,j)) .gt. 1.d-10) then
             if ( j .lt. ib) then
               e(i,j) = hl(el,i,j,rel) - ekin(i,j,rel)
               e(j,i) = 0.d0
             else if (sum(i) .eq. sum(j)) then
               c=hl(el,i,j,rel)-(ekin(i,j,rel)+ekin(j,i,rel))/2.d0
               e(i,j) = c
               e(j,i) = c
             else
               res = 0.d0
               do 14 ii = 1,nwf
                 if (ii.eq.i .or. ii.eq.j) then
                   do 22 k = 0,2*l(i),2
                     if(ii.eq.i) then
                       c= a(i,i,k)-a(j,i,k)-b(j,i,k)
                       if (dabs(c).gt.1.d-10)
     :                          res = res + c*rky(i,i,i,j,k,rel)
		       print *, i,j,k,c,res
                     else if (ii.eq.j) then
                       c= a(j,j,k)-a(i,j,k)-b(i,j,k)
                       if (dabs(c).gt.1.d-10)
     :                         res = res - c*rky(j,j,j,i,k,rel)
		       print *, i,j,k,c,res
                     end if
22                continue
                else
                  do 24 k = 0,2*min0(l(i),l(ii)),2
                    c = a(i,ii,k) - a(j,ii,k)
                    if (dabs(c).gt.1.d-10)
     :                          res = res + c*rky(i,ii,j,ii,k,rel)
                    kk = abs(l(i)-l(ii)) + k
                    c = b(i,ii,kk) - b(j,ii,kk)
                    if(dabs(c) .gt. 1.d-10)
     :                            res = res + c*rky(i,ii,ii,j,kk,rel)
		       print *, i,j,k,c,res
24                continue
                end if
14            continue
              e(i,j) = 2.d0*sum(j)*res/(sum(i)-sum(j))
              e(j,i) = sum(i)*e(i,j)/sum(j)
            end if
          end if
          if (dabs(e(i,j)) .gt. 1.d-10) then
            write(6,35) el(i),el(j),e(i,j),el(j),el(i),e(j,i)
35          format(7x,2(3x,'e(',2a3,') =',f12.5))
          end if
 12     continue
 10   continue
      return
      end
