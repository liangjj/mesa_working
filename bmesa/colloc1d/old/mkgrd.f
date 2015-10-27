*deck mkgrd.f
      subroutine mkgrd(r,rmin,rmax,stp,n,offset,prnt)
      implicit integer (a-z)
      real*8 r, rmin, rmax, stp, tol, del
      parameter ( tol=1.d-08 )
      character*80 title
      logical prnt
      dimension r(*)
      common /io/ inp, iout
      r(1)=rmin
      do 10 i=2,n+offset
         r(i)=r(i-1) + stp
   10 continue
      if (prnt) then
          title='radial grid'
          call prntrm(title,r,n+offset,1,n+offset,1,iout)
      endif
      del=abs(r(n)-rmax)
      if(del.gt.tol) then
         write(iout,*) r(n), rmax
         call lnkerr('quit:bad grid generation')
      endif
      return
      end
