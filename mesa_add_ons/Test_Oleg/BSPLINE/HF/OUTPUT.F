*
*     ------------------------------------------------------------------
*               o u t p u t
*     ------------------------------------------------------------------
*
*       the radial functions and orthogonality integrals are printed,
*   if print is .true.   the  functions  will  also  be  punched  (or
*   stored) on unit ouf, if ouf .ne. 0.
*
      SUBROUTINE output(print)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER(nt=100,kx=8,nwd=10)
      LOGICAL print
      INTEGER ouf
      COMMON /INOUT/ iuf,ouf,iscw
      CHARACTER config*50,el*3,atom*6,term*6
      COMMON /LABEL/config,el(nwd),atom,term
*
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
      COMMON /ORBSCF/sum(nwd), s(nwd), dpm(nwd), acc(nwd), meth(nwd),
     :               iord(nwd), e(nwd,nwd)
*
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
      dimension out(8)
      if ( print ) then
*
*  *****  print radial functions, 7 per page
*
        ml = ib
2       mu = min0(ml+7,nwf)
        i = mu - ml + 1
        mx = 0
        do 1 j = ml,mu
1         mx = max0(mx,max(j))
        write (3,5) atom,term,(el(j),j=ml,mu)
5       format(1h1,9x,19hwave functions for  ,2a6//10x,1hr,8(10x,a3))
        k= 0
        kk = 0
        do 6 j = 1,mx
          do 9 jj = ml,mu
            ij = jj - ml + 1
            out(ij) = p(j,jj)
 9        continue
          k = k+1
          if (k .lt. 10) then
            k = 1
            kk = kk+1
            if (kk .ge. 5) then
              kk = 0
              write (3,23)
23            format(1h1//)
            else
              write (3,10) t(j+ks),(out(jj),jj=1,i)
10            format(/(f13.5,f15.6,7f13.6))
	    end if
          end if
 6      continue
        do 15 j = ml,mu
          ij = j - ml + 1
          out(ij) = dpm(j)
15      continue
        write (3,16) (out(j),j=1,i)
16      format(4x,10hmax. diff. ,f15.7,7f13.7)
        ml = ml+8
        if (ml .le. nwf) go to 2
      end if
30    if ( ouf .ne. 0) then
*
*  *****  output functions on unit ouf for future input
*
        do 3 i = 1,nwf
          mx = max(i)
          write (ouf) atom,term,el(i),mx,z,e(i,i),
     :                (p(j,i),j=1,mx)
3       continue
      end if
*
      end
