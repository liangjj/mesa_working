*deck @(#)mkgrd.f	1.1 9/8/91
      subroutine mkgrd(x,xinv,xs,rmin,rmax,rdel,n,ptbeg,
     1                 last,nospln,toskp)
      implicit integer (a-z)
      real *8 x, xinv, xs, rmin, rmax, rdel, one
      logical nospln
      dimension x(n), xinv(n), xs(*)
      common /io/ inp, iout
      data one/1.d+00/
      x(1)=rmin
      xinv(1)=one/x(1)
      do 10 i=2,n
         x(i)=x(i-1) + rdel
         xinv(i)=one/x(i)
   10 continue
      last=0
      do 20 i=1,n
         if (x(i).gt.rmax) go to 30
   20 continue
      call lnkerr('last point pointer not found')      
   30 last=i
      rmax=x(i)
      if (.not.nospln) then
           nspln=0
           do 40 i=ptbeg,n,toskp
              nspln=nspln+1
              xs(nspln)=x(i)
   40      continue
      endif
      return
      end
