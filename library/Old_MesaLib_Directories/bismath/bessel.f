*deck bessel
c***begin prologue     bessel
c***date written       xxxxxx   (yymmdd)
c***revision date      890422   (yymmdd)
c***keywords           m6004, link 6004, bessel, spline
c***author             schneider, barry (lanl)
c***source             m6004
c***purpose            generate bessel functions for integer values
c***description        
c***references       
c
c***routines called
c***end prologue       bessel
      subroutine bessel(x,j,dj,y,dy,fact,np,nmax,prnt)
      implicit integer (a-z)
      parameter (dim=1000)
      common /io/ inp, iout
      real*8 zero, one, two, epsilon, x, j, dj, y, dy, fact
      real*8 half, pi, twonpi, fournpi, xinv
      real*8 norm, large, renorm
      real*8 besj0, besj1, besy0, besy1, j0, j1, y0, y1
      real*8 jnpl1, jn, jnm1
      logical prnt
      character*80 title
      dimension x(np), j(np,0:nmax), dj(np,0:nmax)
      dimension y(np,0:nmax), dy(np,0:nmax)
      dimension fact(0:*)
      data zero, one, two / 0.d+00, 1.d+00, 2.d+00 /
      data half, large, renorm / .5d0, 1.d+250, 1.d-250  /
      data pi / 3.141592653589793238462643d+00 / 
      twonpi=two/pi
      fournpi=two*twonpi
c----------------------------------------------------------------------c
c            upward recursion, downward recursion or series            c
c----------------------------------------------------------------------c 
      do 10 i=1,np
         xinv=1.d0/x(i)
         j0 = besj0(x(i))
         j1 = besj1(x(i))
         y0 = besy0(x(i))
         y1 = besy1(x(i))
         j(i,0) = j0
         y(i,0) = y0
         dj(i,0) = - j1         
         dy(i,0) = - y1
         if(nmax.gt.0) then
            j(i,1) = j1
            y(i,1) = y1
            dj(i,1) = j(i,0) - j(i,1)*xinv
            dy(i,1)= y(i,0) - y(i,1)*xinv        
         endif
         if(nmax.gt.1) then
c
c           test for up or down recursion
c          
            nn=x(i)
            if(nn.gt.nmax-1) then
c
c                  recur up
c
               j(i,0)=j0
               j(i,1)=j1
               nbeg=1
               do 20 n=2,nmax
                  j(i,n)=(nbeg+nbeg)*j(i,n-1)*xinv - j(i,n-2)
                  nbeg=nbeg+1
 20            continue
c
c              series or downward recursion
c
            else
	       if(abs(x(i)).gt.one) then   
c----------------------------------------------------------------------c
c              find the value of upper needed to accurately get        c
c              the needed n values                                     c
c----------------------------------------------------------------------c
                  upper=nmax
                  upper=msta1(x(i),200)
                  if(upper.lt.nmax) then
                     upper=nmax + nmax
                  else
                     upper=msta2(x(i),nmax,15)
                  endif
                  upper=max(upper,nmax+2)
                  jnpl1=zero
                  jn=renorm
                  onelss=upper-1
                  nbeg=onelss
                  do 30 n=onelss-1,0,-1
                     jnm1 = ( nbeg + nbeg )*jn*xinv 
     1                                   - 
     2                                 jnpl1
                     if(abs(jnm1).gt.large) then
                        jnm1=jnm1*renorm
                        jn=jn*renorm
                     endif                     
                     if(n.le.nmax) then
                        j(i,n) = jnm1
                     endif   
                     jnpl1=jn
                     jn=jnm1
                     nbeg = nbeg -1
 30               continue   
c----------------------------------------------------------------------c
c                 normalize the j                                      c
c----------------------------------------------------------------------c
                  if(abs(j0).gt.abs(j1)) then
                     norm=j0/j(i,0)
                  else
                     norm=j1/j(i,1)
                  endif
                  do 40 n=0,nmax
                     j(i,n) = norm*j(i,n)
 40               continue
               else
                 do 50 n=2,nmax
                    call bser(x(i),j(i,n),dj(i,n),fact,n,0)
 50              continue
               endif   
            endif
c----------------------------------------------------------------------c
c           recur upward for y(i,n)                                    c
c----------------------------------------------------------------------c
            nbeg=1
            do 60 n=2,nmax
               y(i,n) = ( nbeg + nbeg )*y(i,n-1)*xinv - y(i,n-2)
               if(abs(y(i,n)).gt.large) then
                  y(i,n)=large
               endif
               nbeg=nbeg+1
 60         continue
c----------------------------------------------------------------------c
c               get derivatives                                        c
c----------------------------------------------------------------------c
            do 70 n=2,nmax
               dj(i,n) = j(i,n-1) - n*j(i,n)*xinv
               dy(i,n) = y(i,n-1) - n*y(i,n)*xinv
 70         continue
         endif
 10   continue
      if (prnt) then
          title='regular bessel functions'
          call prntfm(title,j,np,nmax+1,np,nmax+1,iout)
          title='irrregular bessel functions'
          call prntfm(title,y,np,nmax+1,np,nmax+1,iout)
          title='derivative of regular bessel functions'
          call prntfm(title,dj,np,nmax+1,np,nmax+1,iout)
          title='derivative of irregular bessel functions'
          call prntfm(title,dy,np,nmax+1,np,nmax+1,iout)
      endif
      return
      end



