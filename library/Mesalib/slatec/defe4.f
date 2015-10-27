*deck defe4
      subroutine defe4 (cofx, idmn, usol, grhs)
c***begin prologue  defe4
c***subsidiary
c***purpose  subsidiary to sepx4
c***library   slatec
c***type      single precision (defe4-s)
c***author  (unknown)
c***description
c
c     this subroutine first approximates the truncation error given by
c     trun1(x,y)=dlx**2*tx+dly**2*ty where
c     tx=afun(x)*uxxxx/12.0+bfun(x)*uxxx/6.0 on the interior and
c     at the boundaries if periodic (here uxxx,uxxxx are the third
c     and fourth partial derivatives of u with respect to x).
c     tx is of the form afun(x)/3.0*(uxxxx/4.0+uxxx/dlx)
c     at x=a or x=b if the boundary condition there is mixed.
c     tx=0.0 along specified boundaries.  ty has symmetric form
c     in y with x,afun(x),bfun(x) replaced by y,dfun(y),efun(y).
c     the second order solution in usol is used to approximate
c     (via second order finite differencing) the truncation error
c     and the result is added to the right hand side in grhs
c     and then transferred to usol to be used as a new right
c     hand side when calling blktri for a fourth order solution.
c
c***see also  sepx4
c***routines called  dx4, dy4
c***common blocks    spl4
c***revision history  (yymmdd)
c   801001  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  defe4
c
      common /spl4/   kswx       ,kswy       ,k          ,l          ,
     1                ait        ,bit        ,cit        ,dit        ,
     2                mit        ,nit        ,is         ,ms         ,
     3                js         ,ns         ,dlx        ,dly        ,
     4                tdlx3      ,tdly3      ,dlx4       ,dly4
      dimension       grhs(idmn,*)           ,usol(idmn,*)
      external cofx
c***first executable statement  defe4
         do  30 i=is,ms
            xi = ait+(i-1)*dlx
            call cofx (xi,ai,bi,ci)
         do 30 j=js,ns
c
c     compute partial derivative approximations at (xi,yj)
c
            call dx4(usol,idmn,i,j,uxxx,uxxxx)
            call dy4(usol,idmn,i,j,uyyy,uyyyy)
            tx = ai*uxxxx/12.0+bi*uxxx/6.0
             ty=uyyyy/12.0
c
c     reset form of truncation if at boundary which is non-periodic
c
            if (kswx.eq.1 .or. (i.gt.1 .and. i.lt.k)) go to  10
            tx = ai/3.0*(uxxxx/4.0+uxxx/dlx)
   10       if (kswy.eq.1 .or. (j.gt.1 .and. j.lt.l)) go to  20
            ty = (uyyyy/4.0+uyyy/dly)/3.0
   20 grhs(i,j)=grhs(i,j)+dly**2*(dlx**2*tx+dly**2*ty)
   30    continue
c
c     reset the right hand side in usol
c
      do  60 i=is,ms
         do  50 j=js,ns
            usol(i,j) = grhs(i,j)
   50    continue
   60 continue
      return
      end
