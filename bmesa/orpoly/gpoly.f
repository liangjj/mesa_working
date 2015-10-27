*deck gpoly.f
c***begin prologue     gpoly
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           orthogonal polynomials
c***author             schneider, barry (nsf)
c***source             math
c***purpose            compute orthogonal polynomials and their first and
c***                   second derivatives on the desired interval.
c***description        the recursion relation used is,
c***
c***                   b p (x) = (x - a ) p   (x) - b   p   (x)
c***                    j j            j   j-1       j-1 j-2
c***                   
c***                   where the p's are normalized and orthogonal.
c***                   the derivatives are simply obtained by 
c***                   differentiating the recursion.
c***references         
c
c***routines called    
c***end prologue       gpoly
      subroutine gpoly(pn,dpn,ddpn,rbeg,rend,x,wts,a,b,norm0,scr,n,
     1                 npts,type,nowts,prnt)
      implicit integer (a-z)
      real*8 pn, dpn, ddpn, rbeg, rend, x, wts, a, b, norm0, scr, dxdy
      real*8 renorm
      character*(*) type
      character*80 title
      logical prnt, nowts
      dimension pn(npts,0:n-1), dpn(npts,0:n-1), ddpn(npts,0:n-1)
      dimension x(npts), wts(n), a(n), b(n), scr(n)
      common/io/inp, iout 

c
      dxdy  = 2.d0/( rend - rbeg )
      renorm=sqrt(dxdy)
c     make first two polynomials, since they are special.
      call vfill(pn(1,0),norm0,npts)
      do 10 i=1,npts
         pn(i,1) = ( x(i) - a(1) )*pn(i,0)/b(1)
 10   continue
c     recur for rest.
      do 20 i=2,n-1
         do 30 j=1,npts
            pn(j,i) = ( ( x(j) - a(i) )*pn(j,i-1)
     1                               - 
     2                           b(i-1)*pn(j,i-2) ) / b(i)
 30      continue   
 20   continue
c     first derivative is zero for p0 for all points
      call rzero(dpn(1,0),npts)
c     make dp1
      do 40 j=1,npts
         dpn(j,1) = dxdy*pn(j,0) / b(1)
 40   continue
c     recur on the rest
      do 50 i=2,n-1
         do 60 j=1,npts
            dpn(j,i) = ( ( x(j) - a(i) )*dpn(j,i-1)
     1                               - 
     2                            b(i-1)*dpn(j,i-2)
     3                               +
     4                              dxdy*pn(j,i-1) ) / b(i)
 60      continue   
 50   continue
c     finish with the second derivatives
c     second derivative is zero for p0 and p1 for all points   
      call rzero(ddpn(1,0),npts)
      call rzero(ddpn(1,1),npts)
c     recur on the rest
      do 70 i=2,n-1
         do 80 j=1,npts
            ddpn(j,i) = ( ( x(j) - a(i) )*ddpn(j,i-1)
     1                                 - 
     2                             b(i-1)*ddpn(j,i-2)
     3                                 +
     4                          2.d0*dxdy*dpn(j,i-1) ) / b(i)
 80      continue   
 70   continue
      call vscale(pn,pn,renorm,n*npts)
      call vscale(dpn,dpn,renorm,n*npts)
      call vscale(ddpn,ddpn,renorm,n*npts)
      title='generating '//type//' polynomials and 1 and 2 derivatives'
      if (type.ne.'hermite'.or.type.ne.'laguerre') then
          write(iout,1) title, rbeg, rend
      else
          if (type.eq.'laguerre') then
              write(iout,2) title
          else
              write(iout,3) title
          endif
      endif
      if(prnt) then
          title='polynomials'
          call prntrm(title,pn(1,0),npts,n,npts,n,iout)
          title='first derivatives'
          call prntrm(title,dpn(1,0),npts,n,npts,n,iout)
          title='second derivatives'
          call prntrm(title,ddpn(1,0),npts,n,npts,n,iout)
      endif
      do 90 i=1,npts
         x(i) = .5d0*( ( rend -rbeg )*x(i) + ( rbeg +rend ) )

 90   continue   
      if(.not.nowts) then
         do 100 i=1,n
            wts(i) = .5d0* ( rend - rbeg )*wts(i)
 100     continue
      endif
      if (prnt) then
          title='modified polynomial points'
          call prntrm(title,x,npts,1,npts,1,iout)
          if(.not.nowts) then
             title='modified polynomial weights'
             call prntrm(title,wts,n,1,n,1,iout)
          endif
      endif
      return
 1    format(/,a80,/,'the interval is left = ',e15.8,1x,'right = ',
     1                                         e15.8)
 2    format(/,a80,/,'the interval is semi-infinite')
 3    format(/,a80,/,'the interval is infinite')
      end       
