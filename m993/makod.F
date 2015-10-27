*deck makod.f
c***begin prologue     makod
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            matrix elements between a bridge function in
c***                   region i - 1 and any function in region i.
c***                   if function in region i is not a bridge function
c***                   then there are no surface terms.  otherwise we need
c***                   to add the bloch operator contribution.
c
c***                   note that normalization factors are passed from 0:nfun
c***                   in this routine.  this allows us to use norm(0), which
c***                   is the normalization of the bridge function in region
c***                   i - 1 to make the needed integral.
c***references         
c
c***routines called    
c***end prologue       makod
      subroutine makod(hmat,kemat,norm,npt,nfun,n,prn)
      implicit integer (a-z)
      real*8 hmat, kemat, norm
      logical prn
      character*80 title
      dimension hmat(*), kemat(npt,*), norm(0:nfun)
      common/io/inp, iout
      do 10 i=1,nfun
         hmat(i) = hmat(i) + norm(i)*norm(0) * kemat(i,1) 
 10   continue
      if(prn) then
         title='matrix elements between region (i-1) and region i'
         call prntrm(title,hmat,nfun,1,nfun,1,iout)
      endif
      return
      end       



