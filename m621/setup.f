*deck @(#)setup.f	5.1  11/6/94
      subroutine setup(x,y,z,q,qd,qq,xcnt,ycnt,zcnt,rad,xobs,yobs,zobs,
     $                 epsilon,npoints,nsrcs,ncnt,nobs,prnt,
     $                 nfocus,thresh)
c***begin prologue     setup.f
c***date written       yymmdd
c***revision date      11/6/94
c
c   8 january, 1994    rlm at lanl
c      incorporating into mesa
c***keywords           solvent
c***author             tawa,greg(lanl)
c***source             @(#)setup.f	5.1   11/6/94
c***purpose            echoes the input data
c***description
c
c
c
c***references
c
c***routines called
c
c***end prologue       setup.f
      implicit none
c     --- input variables -----
      logical prnt
      integer npoints,nsrcs,ncnt,nobs,nfocus
      real*8 epsilon,thresh
c     --- input arrays (unmodified) ---
      real*8 	x(nsrcs), y(nsrcs), z(nsrcs)
      real*8 	q(nsrcs),qd(3,nsrcs),qq(6,nsrcs)
      real*8 	xcnt(ncnt), ycnt(ncnt), zcnt(ncnt)
      real*8 	rad(ncnt)
      real*8 	xobs(nobs), yobs(nobs), zobs(nobs)
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer i
      integer inp,iout
      logical debug
c
      parameter (debug=.false.)
c
      common/io/inp,iout
c
c
   20 format(/2x,'SET-UP FOR BSJ RUN -- '
     &	//4x,'the  value of the dielectric constant is'
     &	/2x, '                   epsilon = ',f12.4,
     &	//4x,'the number of q-points for this'
     & 	/4x,'calculation is'
     &	/2x, '                   npoints = ',i9,
     &				/)
   30 format(/4x,'there are ',i3,' sources'
     &,	/2x,'  i        x(A)        y(A)        z(A)        q/e')
   40 format(2x,i3,4(6x,f10.5))
   42 format(/4x,'molecular volume is the union of spherical cavities',
     &		/4x,'there are ',i3,' cavity spheres'
     &,	/2x,'  i        xcnt(A)     ycnt(A)     zcnt(A)     radius(A)')
   44 format(2x,i3,4(6x,f10.5))
   50 format(/4x,'there are ',i3,' observation points'
     &,	/2x,'  i        xobs(A)     yobs(A)     zobs(A)')
   60 format(2x,i3,3(6x,f10.5))
  100 format(5x,'epsilon                   ',f15.4)
  110 format(5x,'number of points          ',i15)
  120 format(5x,'number of focus points    ',i15)
  130 format(5x,'reaction field convergence',f15.9)
c
c
      if(debug) then
         write(iout,20) epsilon, npoints
         write(iout, 30) nsrcs
         do 5 i = 1, nsrcs
            write(iout,40) i, x(i), y(i), z(i), q(i)
  5      continue
         write(iout, 42) ncnt
         do 10 i = 1, ncnt
            write(iout,44) i, xcnt(i), ycnt(i), zcnt(i), rad(i)
  10     continue
         write(iout, 50) nobs
         do 15 i = 1, nobs
            write(iout,60) i, xobs(i), yobs(i), zobs(i)
   15    continue
         write(iout,'(/4x,''The code will attempt to use '',i10,
     &               '' focus points'')')
     &                  nfocus
      endif
c
c
      write(iout,100) epsilon
      write(iout,110) npoints
      write(iout,120) nfocus
      write(iout,130) thresh
c
c
      return
      end
