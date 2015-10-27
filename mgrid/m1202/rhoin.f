*deck rhoin 
c***begin prologue     rhoin
c***date written       941208   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m1200, link 1200, multigrid
c***author             schneider, b. i.(nsf)
c***source             m1200
c***purpose            input routine for rhs of poisson equation
c***                   and exact solution for model.
c***description        a sample rhs for the poisson equation is
c***                   input here on a uniform grid in (x,y).
c***references         multigrid methods, siam, s. mccormick ed. 
c
c***routines called    iosys, util and mdutil
c***end prologue       rhoin
      subroutine rhoin(rho,h,grid,n,type,prnt)
c
      implicit integer (a-z)
      real*8 rho, h, x, y, xx, yy, f1, f2, f3, f4
      logical prnt
      character*80 title
      character*2 itoc
      character*(*) type
      dimension rho(n,n)
      common/io/inp, iout
      if (type.eq.'model') then
          x=0.d0      
          do 10 i=1,n
             xx=x*x
             f1=1.d0-6.d0*xx
             f2=xx*(1.d0-xx)
             y=0.d0
             do 20 j=1,n
                yy=y*y
                f3=yy*(1.d0-yy)
                f4=1.d0-6.d0*yy
                rho(i,j)=2.d0*( f1*f3 + f4*f2 )
                y=y+h             
 20          continue
             x=x+h
 10       continue
      elseif (type.eq.'exponential') then
          x=0.d0      
          do 30 i=1,n
             xx=x*x
             y=0.d0
             do 40 j=1,n
                yy=y*y
                rho(i,j)=-exp(-(xx+yy))
                y=y+h             
 40          continue
             x=x+h
 30       continue
      elseif (type.eq.'wave') then
          x=0.d0      
          do 50 i=1,n
             f1=sin(5.d0*x)
             y=0.d0
             do 60 j=1,n
                rho(i,j)=f1*sin(5.d0*y)
                y=y+h             
 60          continue
             x=x+h
 50       continue
      else
          call lnkerr('error in rhs definition')
      endif
      if(prnt) then
         title='input rho on grid = '//itoc(grid)
         call prntrm(title,rho,n,n,n,n,iout)
      endif            
      return
      end





