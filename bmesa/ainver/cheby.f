*deck cheby
      subroutine cheby(p,dp,ddp,x,wtx,y,left,right,sx,cx,
     1                   sy,cy,nx,ny,n,prnt)
c***begin prologue     cheby
c***date written       940504   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose             
c***                   .
c***description
c***            
c               
c               
c***references
c
c***routines called
c
c***end prologue       cheby
c
      implicit integer (a-z)
      real*8 p, dp, ddp, x, y, wtx, sx, cx, sy, cy, left, right
      real*8 add, pi, sumwt, prefac, norm, arg 
      character*80 title
      character*3 itoc
      logical prnt
      dimension p(ny,nx), dp(ny,nx), ddp(ny,nx), x(nx), wtx(nx)
      dimension sx(nx,n), cx(nx,n), y(ny), sy(ny,n), cy(ny,n)
      common /io/ inp, iout
      data pi/3.1415926535897932384d0/
      call rzero(p,nx*ny)
      call rzero(dp,nx*ny)
      call rzero(ddp,nx*ny)
      add=(right-left)/(n+1)
      prefac = pi/(right-left)
      norm=sqrt(2.d0/(right-left))
      do 20 i=1,n
         do 30 j=1,ny
            arg = i*prefac*( y(j) - left )
            sy(j,i) = sin(arg )
            cy(j,i) = cos(arg)
 30      continue
 20   continue
      call vscale(sy,sy,norm,ny*n)
      call vscale(cy,cy,norm,ny*n)
      do 40 i=1,nx
         do 50 j=1,ny
            do 60 k=1,n
               p(j,i) = p(j,i) + sy(j,k)*sx(i,k)
               dp(j,i) = dp(j,i) + k*prefac*cy(j,k)*sx(i,k)
               ddp(j,i) = ddp(j,i) - k*k*prefac*prefac*sy(j,k)*sx(i,k)
 60         continue
 50      continue
 40   continue   
      add=sqrt( (right-left) / (n+1) )
      call vscale(p,p,add,nx*ny)
      call vscale(dp,dp,add,nx*ny)
      call vscale(ddp,ddp,add,nx*ny)
      if(prnt) then
         title='DVR box functions'
         call prntrm(title,p,ny,nx,ny,nx,iout)
         title='first derivative of DVR box functions'
         call prntrm(title,dp,ny,nx,ny,nx,iout)
         title='second derivative of DVR box functions'
         call prntrm(title,ddp,ny,nx,ny,nx,iout)
      endif
      return
      end















