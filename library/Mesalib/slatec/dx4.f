*deck dx4
      subroutine dx4 (u, idmn, i, j, uxxx, uxxxx)
c***begin prologue  dx4
c***subsidiary
c***purpose  subsidiary to sepx4
c***library   slatec
c***type      single precision (dx4-s)
c***author  (unknown)
c***description
c
c     this program computes second order finite difference
c     approximations to the third and fourth x
c     partial derivatives of u at the (i,j) mesh point.
c
c***see also  sepx4
c***routines called  (none)
c***common blocks    spl4
c***revision history  (yymmdd)
c   801001  date written
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  dx4
c
      common /spl4/   kswx       ,kswy       ,k          ,l          ,
     1                ait        ,bit        ,cit        ,dit        ,
     2                mit        ,nit        ,is         ,ms         ,
     3                js         ,ns         ,dlx        ,dly        ,
     4                tdlx3      ,tdly3      ,dlx4       ,dly4
      dimension       u(idmn,*)
c***first executable statement  dx4
      if (i.gt.2 .and. i.lt.(k-1)) go to  50
      if (i .eq. 1) go to  10
      if (i .eq. 2) go to  30
      if (i .eq. k-1) go to  60
      if (i .eq. k) go to  80
c
c     compute partial derivative approximations at x=a
c
   10 if (kswx .eq. 1) go to  20
      uxxx = (-5.0*u(1,j)+18.0*u(2,j)-24.0*u(3,j)+14.0*u(4,j)-
     1                                               3.0*u(5,j))/(tdlx3)
      uxxxx = (3.0*u(1,j)-14.0*u(2,j)+26.0*u(3,j)-24.0*u(4,j)+
     1                                      11.0*u(5,j)-2.0*u(6,j))/dlx4
      return
c
c     periodic at x=a
c
   20 uxxx = (-u(k-2,j)+2.0*u(k-1,j)-2.0*u(2,j)+u(3,j))/(tdlx3)
      uxxxx = (u(k-2,j)-4.0*u(k-1,j)+6.0*u(1,j)-4.0*u(2,j)+u(3,j))/dlx4
      return
c
c     compute partial derivative approximations at x=a+dlx
c
   30 if (kswx .eq. 1) go to  40
      uxxx = (-3.0*u(1,j)+10.0*u(2,j)-12.0*u(3,j)+6.0*u(4,j)-u(5,j))/
     1       tdlx3
      uxxxx = (2.0*u(1,j)-9.0*u(2,j)+16.0*u(3,j)-14.0*u(4,j)+6.0*u(5,j)-
     1                                                      u(6,j))/dlx4
      return
c
c     periodic at x=a+dlx
c
   40 uxxx = (-u(k-1,j)+2.0*u(1,j)-2.0*u(3,j)+u(4,j))/(tdlx3)
      uxxxx = (u(k-1,j)-4.0*u(1,j)+6.0*u(2,j)-4.0*u(3,j)+u(4,j))/dlx4
      return
c
c     compute partial derivative approximations on the interior
c
   50 continue
      uxxx = (-u(i-2,j)+2.0*u(i-1,j)-2.0*u(i+1,j)+u(i+2,j))/tdlx3
      uxxxx = (u(i-2,j)-4.0*u(i-1,j)+6.0*u(i,j)-4.0*u(i+1,j)+u(i+2,j))/
     1        dlx4
      return
c
c     compute partial derivative approximations at x=b-dlx
c
   60 if (kswx .eq. 1) go to  70
      uxxx = (u(k-4,j)-6.0*u(k-3,j)+12.0*u(k-2,j)-10.0*u(k-1,j)+
     1                                                 3.0*u(k,j))/tdlx3
      uxxxx = (-u(k-5,j)+6.0*u(k-4,j)-14.0*u(k-3,j)+16.0*u(k-2,j)-
     1                                     9.0*u(k-1,j)+2.0*u(k,j))/dlx4
      return
c
c     periodic at x=b-dlx
c
   70 uxxx = (-u(k-3,j)+2.0*u(k-2,j)-2.0*u(1,j)+u(2,j))/tdlx3
      uxxxx = (u(k-3,j)-4.0*u(k-2,j)+6.0*u(k-1,j)-4.0*u(1,j)+u(2,j))/
     1        dlx4
      return
c
c     compute partial derivative approximations at x=b
c
   80 uxxx = -(3.0*u(k-4,j)-14.0*u(k-3,j)+24.0*u(k-2,j)-18.0*u(k-1,j)+
     1                                                 5.0*u(k,j))/tdlx3
      uxxxx = (-2.0*u(k-5,j)+11.0*u(k-4,j)-24.0*u(k-3,j)+26.0*u(k-2,j)-
     1                                    14.0*u(k-1,j)+3.0*u(k,j))/dlx4
      return
      end
