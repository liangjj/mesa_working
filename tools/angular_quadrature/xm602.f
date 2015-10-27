*deck xm602
      program xm602
      implicit none
c
      real*8 p(3,1000),w(1000)
      real*8 q(100)
      logical gauss,lebdev,debug
      integer npts,i,inp,iout,lmax
      integer nleb(23)
      data nleb/2*0,6,5*0,38,0,50,0,74,0,86,0,110,0,146,0,0,0,194/
      data gauss/.false./,lebdev/.true./
      parameter (debug=.false.)
      common/io/inp,iout
c
c     set the order of the quadrature here, will check up to l+l'=9
      lmax=3
      if(gauss) then
         npts=999
         call gauleg(p,w,lmax,npts)
      else if(lebdev) then
         npts=nleb(lmax)
         if(npts.eq.0) write(0,*) 'nth order lebedev not available',lmax
         call lebedev(p,w,lmax,npts)
      endif
      if(debug) then
         do 10 i=1,npts
            write(iout,*) p(1,i),p(2,i),p(3,i),w(i)
   10    continue
      endif
c
c
      call quad(p,w,npts,lmax)
c
c     murray's radial grid
c     call radial(q,w,10,1.0d+00)
c
      stop
      end
