*deck diagx.f
c***begin prologue     diagx
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***purpose            transformation to representation where co-ordinate
c***                   is diagonal.
c***                   
c***                                                          
c***references         
c
c***routines called    
c***end prologue       diagx
      subroutine diagx(p,dp,ddp,a,b,pn,dpn,ddpn,d,e,dum,n,
     1                 npts,prn,coord,grid)
      implicit integer (a-z)
      real*8 p, dp, ddp, a, b, pn, dpn, ddpn, d, e, dum
      logical prn
      character*3 itoc
      character*80 title
      dimension p(npts,0:n-1), dp(npts,0:n-1), ddp(npts,0:n-1)
      dimension a(n), b(n) 
      dimension pn(npts,0:n-1), dpn(npts,0:n-1), ddpn(npts,0:n-1)
      dimension d(n), e(n), dum(n,n)
      common/io/inp, iout 
c     diagonalize the tridiagonal matrix in the polynomial basis in order to
c     obtain the eigenvalues and eigenvectors of the co-ordinate operator.
      do 10 i=1,n
         d(i)=a(i)
 10   continue
      do 20 i=1,n-1
         e(i+1)=b(i)
 20   continue
      call rzero(dum,n*n)
      do 30 i=1,n
         dum(i,i)=1.d0
 30   continue
      call imtql2(n,n,d,e,dum,ierr)
      call iosys('write real "coordinate eigenvectors for '// 
     1           'coordinate '//itoc(coord)//' and grid '//
     2            itoc(grid)//'" to lamdat',n*n,dum,0,' ')
      title='eigenvalues of co-ordinate representation'
      call prntrm(title,d,n,1,n,1,iout)
      if (prn) then
          title='eigenvectors of co-ordinate representation'
          call prntrm(title,dum,n,n,n,n,iout)
      endif
c     transform the polynomials and their derivatives to the new
c     representation      
      call ebc(pn(1,0),p(1,0),dum,npts,n,n)
      call ebc(dpn(1,0),dp(1,0),dum,npts,n,n)
      call ebc(ddpn(1,0),ddp(1,0),dum,npts,n,n)
      return
      end       
