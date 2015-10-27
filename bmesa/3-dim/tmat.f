*deck tmat.f
c***begin prologue     tmat
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           kinetic energy
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            form one-dimensional kinetic energy matrix.  
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       tmat
      subroutine tmat(p,dp,ddp,t,wts,der,nply,npts,prnh,prnv)
      implicit integer (a-z)
      real*8 p, dp, ddp, wts, der, t
      character*80 title
      logical prnh, prnv
      dimension p(npts,0:nply-1), dp(npts,0:nply-1), ddp(npts,0:nply-1)
      dimension wts(npts), t(nply,nply), der(2)
      common/io/inp, iout 
c     calculate kinetic energy plus bloch operator portion of hamiltonian 
c     in one dimension.
      call rzero(t,nply*nply)
      do 10 i=1,nply
         do 20 j=1,i
            do 30 k=1,npts
               t(i,j) = t(i,j) - .5d0*wts(k)*p(k,i-1)*ddp(k,j-1)
 30         continue
            t(i,j) = t(i,j) +.5d0 *
     1            ( p(npts,i-1)*( dp(npts,j-1) -der(2)*p(npts,j-1) ) 
     2                          - 
     3              p(1,i-1)*( dp(1,j-1) - der(1)*p(1,j-1) ) )
            t(j,i)=t(i,j)
 20      continue
 10   continue
      if(prnh) then
         title='kinetic energy matrix'   
         call prntrm(title,t,nply,nply,nply,nply,iout)
      endif
      return
      end       

