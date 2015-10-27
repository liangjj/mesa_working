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
      real*8 rho, h, x, y, xx, yy, z, zz, f1, f2, f3, f4, f5, f6
      logical prnt
      character*80 title
      character*2 itoc
      character*(*) type
      dimension rho(n,n,n)
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
                z=0.d0
                do 30 k=1,n
                   zz=z*z
                   f5=zz*(1.d0-zz)
                   f6=1.d0-6.d0*zz
                   rho(i,j,k)=2.d0*( f1*f3*f5 + f4*f2*f5 + f6*f2*f3 )
                   z=z+h
 30             continue                   
                y=y+h             
 20          continue
             x=x+h
 10       continue
      elseif (type.eq.'exponential') then
          x=0.d0      
          do 40 i=1,n
             xx=x*x
             y=0.d0
             do 50 j=1,n
                yy=y*y
                z=0.d0
                do 60 k=1,n
                   zz=z*z
                   rho(i,j,k)=-exp(-(xx+yy+zz))
 60             continue                   
                y=y+h             
 50          continue
             x=x+h
 40       continue
      elseif (type.eq.'wave') then
          x=0.d0      
          do 70 i=1,n
             f1=sin(5.d0*x)
             y=0.d0
             do 80 j=1,n
                f2=sin(5.d0*y)
                z=0.d0
                do 90 k=1,n
                   rho(i,j,k)=f1*f2*sin(5.d0*z)
                   z=z+h             
 90             continue
                y=y+h
 80          continue                
             x=x+h
 70       continue
      else
          call lnkerr('error in rhs definition')
      endif
      if(prnt) then
         title='input rho on grid = '//itoc(grid)
         call prntrm(title,rho,n,n,n,n,iout)
      endif            
      return
      end





