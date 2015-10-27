*deck dy4
      subroutine dy4 (u, idmn, i, j, uyyy, uyyyy)
c***begin prologue  dy4
c***subsidiary
c***purpose  subsidiary to sepx4
c***library   slatec
c***type      single precision (dy4-s)
c***author  (unknown)
c***description
c
c     this program computes second order finite difference
c     approximations to the third and fourth y
c     partial derivatives of u at the (i,j) mesh point.
c
c***see also  sepx4
c***routines called  (none)
c***common blocks    spl4
c***revision history  (yymmdd)
c   801001  date written
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  dy4
c
      common /spl4/   kswx       ,kswy       ,k          ,l          ,
     1                ait        ,bit        ,cit        ,dit        ,
     2                mit        ,nit        ,is         ,ms         ,
     3                js         ,ns         ,dlx        ,dly        ,
     4                tdlx3      ,tdly3      ,dlx4       ,dly4
      dimension       u(idmn,*)
c***first executable statement  dy4
      if (j.gt.2 .and. j.lt.(l-1)) go to  50
      if (j .eq. 1) go to  10
      if (j .eq. 2) go to  30
      if (j .eq. l-1) go to  60
      if (j .eq. l) go to  80
c
c     compute partial derivative approximations at y=c
c
   10 if (kswy .eq. 1) go to  20
      uyyy = (-5.0*u(i,1)+18.0*u(i,2)-24.0*u(i,3)+14.0*u(i,4)-
     1                                                 3.0*u(i,5))/tdly3
      uyyyy = (3.0*u(i,1)-14.0*u(i,2)+26.0*u(i,3)-24.0*u(i,4)+
     1                                      11.0*u(i,5)-2.0*u(i,6))/dly4
      return
c
c     periodic at x=a
c
   20 uyyy = (-u(i,l-2)+2.0*u(i,l-1)-2.0*u(i,2)+u(i,3))/tdly3
      uyyyy = (u(i,l-2)-4.0*u(i,l-1)+6.0*u(i,1)-4.0*u(i,2)+u(i,3))/dly4
      return
c
c     compute partial derivative approximations at y=c+dly
c
   30 if (kswy .eq. 1) go to  40
      uyyy = (-3.0*u(i,1)+10.0*u(i,2)-12.0*u(i,3)+6.0*u(i,4)-u(i,5))/
     1       tdly3
      uyyyy = (2.0*u(i,1)-9.0*u(i,2)+16.0*u(i,3)-14.0*u(i,4)+6.0*u(i,5)-
     1                                                      u(i,6))/dly4
      return
c
c     periodic at y=c+dly
c
   40 uyyy = (-u(i,l-1)+2.0*u(i,1)-2.0*u(i,3)+u(i,4))/tdly3
      uyyyy = (u(i,l-1)-4.0*u(i,1)+6.0*u(i,2)-4.0*u(i,3)+u(i,4))/dly4
      return
c
c     compute partial derivative approximations on the interior
c
   50 continue
      uyyy = (-u(i,j-2)+2.0*u(i,j-1)-2.0*u(i,j+1)+u(i,j+2))/tdly3
      uyyyy = (u(i,j-2)-4.0*u(i,j-1)+6.0*u(i,j)-4.0*u(i,j+1)+u(i,j+2))/
     1        dly4
      return
c
c     compute partial derivative approximations at y=d-dly
c
   60 if (kswy .eq. 1) go to  70
      uyyy = (u(i,l-4)-6.0*u(i,l-3)+12.0*u(i,l-2)-10.0*u(i,l-1)+
     1                                                 3.0*u(i,l))/tdly3
      uyyyy = (-u(i,l-5)+6.0*u(i,l-4)-14.0*u(i,l-3)+16.0*u(i,l-2)-
     1                                     9.0*u(i,l-1)+2.0*u(i,l))/dly4
      return
c
c     periodic at y=d-dly
c
   70 continue
      uyyy = (-u(i,l-3)+2.0*u(i,l-2)-2.0*u(i,1)+u(i,2))/tdly3
      uyyyy = (u(i,l-3)-4.0*u(i,l-2)+6.0*u(i,l-1)-4.0*u(i,1)+u(i,2))/
     1        dly4
      return
c
c     compute partial derivative approximations at y=d
c
   80 uyyy = -(3.0*u(i,l-4)-14.0*u(i,l-3)+24.0*u(i,l-2)-18.0*u(i,l-1)+
     1                                                 5.0*u(i,l))/tdly3
      uyyyy = (-2.0*u(i,l-5)+11.0*u(i,l-4)-24.0*u(i,l-3)+26.0*u(i,l-2)-
     1                                    14.0*u(i,l-1)+3.0*u(i,l))/dly4
      return
      end
