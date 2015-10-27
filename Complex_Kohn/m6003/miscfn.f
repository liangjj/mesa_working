*deck @(#)miscfn.f	1.1 9/8/91
      subroutine miscfn(grid,x,cphi,sphi,npnt)
      implicit integer (a-z)
      real *8 grid, r, x, cphi, sphi
      dimension grid(4,npnt), x(npnt), cphi(npnt)
      dimension sphi(npnt)
      do 10 pnt=1,npnt
         r=sqrt(grid(1,pnt)*grid(1,pnt)+grid(2,pnt)*grid(2,pnt)+
     1             grid(3,pnt)*grid(3,pnt))
         x(pnt)=grid(3,pnt)/r
         cphi(pnt)=grid(1,pnt)/(r*sqrt(1.d+00-x(pnt)*x(pnt)))
         sphi(pnt)=grid(2,pnt)/(r*sqrt(1.d+00-x(pnt)*x(pnt)))
   10 continue
      return
      end
