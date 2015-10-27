*deck @(#)surface.f	5.2 5/22/95
      subroutine surface(xyzch,xcnt,ycnt,zcnt,rad,
     $         npoints,ncnt,prnt,nfocmx,mark,nark,
     $         xq,yq,zq,indx,u,kount,
     $         nfocus,ifocus,icindx,ipindx,index1,index2,index3)
c***begin prologue     surface.f
c***date written       yymmdd
c***revision date      5/22/95
c
c***keywords
c***author             tawa,greg(lanl)
c***source             @(#)surface.f	5.2   5/22/95
c***purpose
c***description
c
c
c
c***references
c
c***routines called
c
c***end prologue       surface.f
      implicit none
c     --- input variables -----
      integer npoints,ncnt,nfocmx
      logical prnt
c     --- input arrays (unmodified) ---
      real*8 xcnt(ncnt),ycnt(ncnt),zcnt(ncnt),rad(ncnt)
c     --- input arrays (scratch) ---
      integer mark(ncnt+1),nark(ncnt+1)
      integer ifocus(ncnt)
      real*8 u(2*ncnt)
c     --- output arrays ---
      integer index1(nfocmx)
      integer index2(npoints),index3(nfocmx)
      integer icindx(npoints),ipindx(npoints)
      integer indx(npoints)
      real*8 xyzch(3,nfocmx)
      real*8 xq(npoints),yq(npoints),zq(npoints)
c     --- output variables ---
      integer nfocus
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer i,j,icnt,jcnt,ndim,kount
      integer ncntp1,nfocus2,npo,itest
      integer n,iadd,iadd2,icomp,ido,ioff,ioff1
      integer jmin
      real*8 zero,one,two,four
      real*8 pi,twopi,fourpi
      real*8 f,ux,uy,uz
      real*8 atot,d,old
      real*8 r,hudge,dx,dy,dz
      parameter (hudge=1.0d+10)
      logical buried,timing
      logical debug
c
      data zero/0.0d+00/,one/1.0d+00/,two/2.0d+00/,four/4.0d+00/
      save zero,one,two,four
c
      parameter (timing=.false.)
      parameter (debug=.false.)
c
      common/io/inp,iout
c
c     --- set up some constants and times ---
      pi=four*atan(one)
      twopi=two*pi
      fourpi=four*pi
c
c     --- we want to distribute the points among the sphere surfaces
c         approximately proportionally to the surface
c         area of the spheres -- thus we want the
c         area fraction for each sphere and to get this we
c         first find the total area of the sphere surfaces
      atot =  zero
      do 5 icnt = 1, ncnt
         atot = atot + rad(icnt)**2
  5   continue
      mark(1) = 1
      do 10 icnt = 2, ncnt
         mark(icnt) =  mark(icnt-1)
     &              +int(float(npoints)*rad(icnt-1)**2/atot)
  10  continue
      n = ncnt+1
      mark(n) = npoints+1
c
c     --- get quadrature points ---
      do 15 icnt = 1, ncnt
         ndim = -2
	 call sobseq1(ndim, u)
	 ndim = -ndim
	 kount = 0
	 do 20 i = mark(icnt), mark(icnt+1)-1
            kount = kount + 1
            call sobseq1(ndim,u)
            ux = u(1)
            uz = u(2)
c           cos theta
            uz = two*uz - one
c           sin theta
            f = sqrt(one-uz*uz)
c           phi
            ux = twopi*ux
c           y
            uy = f*sin(ux)		
c           x
            ux = f*cos(ux)
            xq(i) = xcnt(icnt) + rad(icnt)*ux
            yq(i) = ycnt(icnt) + rad(icnt)*uy
            zq(i) = zcnt(icnt) + rad(icnt)*uz
  20     continue
  15  continue
      if(debug) then
         write(iout,22)
  22     format(/4x,'2d SOBOL SERIES SAMPLING OF THE SPHERES'/)
      endif
c
c     --- determine all buried points
      kount = 0
      do 25 icnt = 1, ncnt
         nark(icnt) = kount + 1
         do 30 i = mark(icnt), mark(icnt+1)-1
	    buried = .false.
            do 35 jcnt = 1, ncnt
               if(jcnt.ne.icnt)then
                  d=(xq(i)-xcnt(jcnt))**2
     &             +(yq(i)-ycnt(jcnt))**2
     &             +(zq(i)-zcnt(jcnt))**2
                  if(d.lt.rad(jcnt)**2) then
                     buried =.true.
                     goto 77
                  endif
               endif
  35        continue
  77        continue
            if(.not.buried) then
               kount = kount + 1
               indx(kount) = i
               ipindx(kount)=
     &                     float(mark(icnt+1)-mark(icnt))
               icindx(kount)=icnt
            endif
  30     continue
  25  continue
      ncntp1 = ncnt + 1
      nark(ncntp1) = kount + 1
c
      if(debug) then
         write(iout,88)kount
   88    format(/4x,'THERE ARE ',i9,' POINTS EXPOSED ON THE SURFACE'/)
      endif
c
c     --- contract the list of points so that only non-buried points are
c         included in the xq, yq, zq vectors.
      do 40 i = 1, kount
         xq(i) = xq(indx(i))
         yq(i) = yq(indx(i))
         zq(i) = zq(indx(i))
  40  continue
c
c      --- Determine focus points.  Assign number of focus points
c          proportional to the surface area of atom (iopt = 1).
c          Assign number of focus points proportional to the number
c          of exposed points on the surface area of the atom (iopt=2)
c          Determine how many focus points per sphere.
      nfocus2 = 0
      do 45 icnt=1,ncnt
         npo = nark(icnt+1)-nark(icnt)
         ifocus(icnt)=nfocmx*npo/kount
         itest = nark(icnt+1)-nark(icnt)-ifocus(icnt)
         if(debug) then
            write(iout,'(/4x,''For center = '',t35,i10)') icnt
            write(iout,'(4x,''Trial nfocus = '',t35,i10)') ifocus(icnt)
            write(iout,'(4x,''Number of non-buried points = '',
     &                 t35,i10)') nark(icnt+1)-nark(icnt)
         endif
         if (itest .le. 0) then
            if(debug) then
               write(iout,'(4x,''nfocus too large, scaling back'')')
            endif
            ifocus(icnt)=ifocus(icnt)+itest
         endif
         if(debug) then
            write(iout,'(4x,''Final nfocus = '',t35,i10)') ifocus(icnt)
         endif
         nfocus2= nfocus2+ ifocus(icnt)
  45  continue
      if(debug) then
         write(iout,'(/4x,''Total number of focus points = '',
     &                 t35,i10)') nfocus2
      endif
c****greg,did i do this right
      if (nfocus2 .lt. nfocmx .and. kount .gt. nfocmx) then
         if(debug) then
            write(iout,
     $            '(/4x,''Total number of focus point < input,'
     $              //' adding'')')
         endif
c        number of focus points to be added
         iadd = nfocmx - nfocus2
c        number of focus points added so far
         iadd2 = 0
  54     continue
         do 55 icnt = 1,ncnt
            icomp = nark(icnt+1)-nark(icnt)
            if (icomp .gt. 0 .and. iadd2 .lt. iadd) then
               ifocus(icnt) = ifocus(icnt)+1
               iadd2 = iadd2 + 1
            endif
  55     continue
         if (iadd2 .lt. iadd) then
            goto 54
         endif
         if(debug) then
            write(iout,
     $        '(/4x,''Sphere stats after focus point addition    '')')
         endif
         nfocus = 0
         do 56 icnt = 1,ncnt
            nfocus = nfocus + ifocus(icnt)
            if(debug) then
               write(iout,'(/4x,''For Center = '',t35,i10)') icnt
               write(iout,'(4x,''Number of focus points = '',t35,i10)')
     &               ifocus(icnt)
               write(iout,
     $               '(4x,''Number of non-buried points = '',t35,i10)')
     &               nark(icnt+1)-nark(icnt)
            endif
  56     continue
         if(debug) then
            write(iout,'(/4x,''Total number of focus points = '',
     &                 t35,i10)') nfocus
         endif
      else
         nfocus = nfocus2
         if(nfocus.gt.nfocmx) then
            write(iout,*) 'nfocus,nfocmx',nfocus,nfocmx
            call lnkerr('must increase number of focus points')
         endif
      endif
c
c     --- Set the identity of the focus points
c         index1(npoints) is an index which determines which of the 
c         xq,yq,zq sets belong to a focus point.  If index1(i) = j then focus 
c         point i is SOBOL point j.
      ido = 0
      do 60 icnt = 1,ncnt
         iadd = 0
         do 65 i=nark(icnt), nark(icnt+1)-1
            iadd=iadd+1
            if (iadd .le. ifocus(icnt)) then
               ido =  ido + 1
               index1(ido)= i
            endif
  65     continue
  60  continue
c
c     --- Set identity of the SOBOL points that belong to each focus 
c         point, i.e., that are in the same plaquette as the SOBOL points.
c         this is determined by index2 (i).  If index2(i) = j and i=j, the focus
c         point is associated SOBOL point i, if index2(i) = j this means
c         that SOBOL point j is associated with focus point i.
c      do 70 i = 1,kount
c         old = hudge
c         do 75 j = 1,nfocus
c            dx = xq(i) - xq(index1(j))
c            dy = yq(i) - yq(index1(j))
c            dz = zq(i) - zq(index1(j))
c             r = sqrt(dx*dx+dy*dy+dz*dz)
c            if ( r .lt. old) then
c               old = r
c               jmin = j
c            endif
c 75     continue
c        index2(i) = jmin
c 70  continue
c
c
      ioff = 1
      do 70 icnt = 1,ncnt
         do 75 i = nark(icnt),nark(icnt+1)-1
            old = hudge
            do 76 j = 1+ioff1,ifocus(icnt)+ioff1
               dx = xq(i) - xq(index1(j))
               dy = yq(i) - yq(index1(j))
               dz = zq(i) - zq(index1(j))
               r = sqrt(dx*dx+dy*dy+dz*dz)
               if ( r .lt. old) then
                  old = r
                  jmin = j
               endif
  76        continue
            index2(i) = jmin
  75     continue
         ioff1 = ioff1 + ifocus(icnt)
  70  continue
c     Example:
c     The scenario below describes the situation for which SOBOL points
c     1 and 3 are FOCUS points.  SOBOL point 2 is associated with FOCUS
c     point 1 and SOBOL point 4 is associated with FOCUS point 2.
c     SOBOL points        Focus points    SOBOL points
c     xq(1)1------------->index1(1)       index2(1) = 1
c     xq(2)                               index2(2) = 1
c     xq(3)3------------->index1(2)       index2(3) = 3
c     xq(4)                               index2(4) = 3
c     .
c     .
c     .
c     xq(kount)          index1(nfocus)   index2(kount)
c
c     --- Now calculate the number of points per plaquetts and store in 
c         index3.
      do 80 i = 1,nfocus
         index3(i) = 0
         do 85 j = 1,kount
            if (index2(j) .eq. i) then
               index3(i)= index3(i) +1
            endif
  85     continue
  80  continue
c
c     --- pack the focus points into an output array
      do 100 i=1,nfocus
         xyzch(1,i)=xq(index1(i))
         xyzch(2,i)=yq(index1(i))
         xyzch(3,i)=zq(index1(i))
  100 continue
c
c
      return
      end
