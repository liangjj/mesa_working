*deck @(#)sizes.f	5.2 2/5/95
      subroutine sizes(ngrid,radshls,dengrid,wts,ptrad,ndmat,
     $                 pr,dofr,y2,u,scr,grid,mxgrd,charge,size)
c***begin prologue     sizes.f
c***date written       yymmdd  
c***revision date      2/5/95      
c
c***keywords           
c***author             
c***source             @(#)sizes.f	5.2   2/5/95
c***purpose            estimates atomic radius
c***description
c
c***references
c
c***routines called
c
c***end prologue       sizes.f
      implicit none
c     --- input variables -----
      integer ngrid,radshls,ndmat,mxgrd
      real*8 charge
c     --- input arrays (unmodified) ---
      integer ptrad(0:radshls)
      real*8 grid(mxgrd,3)
      real*8 dengrid(ngrid,ndmat),wts(ngrid)
c     --- input arrays (scratch) ---
      real*8 pr(0:radshls)
      real*8 y2(0:radshls),u(0:radshls)
c     --- output arrays ---
c     --- output variables ---
      real*8 size
c     --- scratch arrays ---
      real*8 dofr(0:radshls),scr(ngrid)
c     --- local variables ---
      integer r,nanglr
      integer inp,iout
      integer nreal,j,nfound,i
      real*8 sdot,big,toler,zero,one,two,six,three
      real*8 t1,t2,t3,t4,h
      real*8 a1,a2,a3,roots(3)
c
      common/io/inp,iout
c
      parameter (zero=0.0d+00,big=1.0d+31,toler=0.05d+00)
      parameter (one=1.0d+00,two=2.0d+00,three=3.0d+00,six=6.0d+00)
c
c     --- determine an estimate of the 'atomic' size
c         this is done by integrating the density over the unit sphere
c         for each radial shell, then spline fitting the resulting
c         radial distribution.
c
c     --- generate the total density as a function of r.
c         the closed shell piece is multiplied by two.
      call vadd(scr,dengrid(1,1),dengrid(1,1),ngrid)
      if(ndmat.eq.2) then
         call vadd(scr,scr,dengrid(1,2),ngrid)
      endif
      dofr(0)=zero
      pr(0)=zero
      do 100 r=1,radshls
c        find radius of each shell.
         pr(r)=sqrt(grid(ptrad(r-1),1)*grid(ptrad(r-1),1)
     $             +grid(ptrad(r-1),2)*grid(ptrad(r-1),2)
     $             +grid(ptrad(r-1),3)*grid(ptrad(r-1),3))
         nanglr=ptrad(r)-ptrad(r-1)
         dofr(r)=dofr(r-1)
     $          +sdot(nanglr,scr(ptrad(r-1)),1,wts(ptrad(r-1)),1)
  100 continue
c
c     --- cubic spline fit to the resulting distribution.
c         note that the first derivative of the interpolating function is
c         set to zero at the origin and last point.
      call spline3(pr,dofr,radshls,zero,zero,y2,u)
c
c     --- march out in r until all but .05e- are accounted for.
c         the variable charge contains the complete integral.
      do 160 r=1,radshls
         if(charge-dofr(r).le.toler) then
            j=r
            go to 165
         endif
  160 continue
  165 continue
c      --- the cutoff value occurs between radial point j and j-1.
c          in order to get a finer estimate, use the cubic spline fit.
c          must first turn the fit into a cubic equation of the form
c          x**3 +a1*x**2 +a2*x +a3 = 0.  see "numerical recipes".
      h=pr(j)-pr(j-1)
c     coefficient in front of B terms (B=(x-x(j))/(x(j+1)-x(j))
c     B**0
      t1=(charge-toler)-dofr(j-1)
c     b**1
      t2=dofr(j)-dofr(j-1)
      t2=t2-(two*h*h*y2(j-1)/six)
      t2=t2-(h*h*y2(j)/six)
c     B**2
      t3=three*h*h*y2(j-1)/six
c     B**3
      t4=(y2(j)-y2(j-1))*h*h/six
c
c     --- now cast these into the a coefficients above.
      a1=t3/t4
      a2=t2/t4
      a3=-t1/t4
      call cubic(a1,a2,a3,roots,nreal)
c
c     --- test the solutions; we want the one between 0.0 and 1.0.
      if(nreal.eq.1) then
         if(roots(1).lt.zero.or.roots(1).gt.one) then
            call lnkerr('problem with cubic solution in sizes')
         endif
         size=pr(j-1)+h*roots(1)
      else if(nreal.eq.3) then
         nfound=0
         do 200 i=1,nreal
            if(roots(i).lt.zero.or.roots(i).gt.one) then
c              outside the range.
            else
               nfound=nfound+1
               roots(nfound)=roots(i)
            endif
  200    continue
         if(nfound.ne.1) then
            call lnkerr('problem with cubic solution in sizes')
         endif
         size=pr(j-1)+h*roots(1)
      else
         call lnkerr('problem with cubic solution in sizes')
      endif
c
c
      return
      end
