*deck makd.f
c***begin prologue     makd
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            matrix elements between functions all in
c***                   region i.
c***references         
c
c***routines called    
c***end prologue       makd
      subroutine makd(hmat,kemat,keadd,norm,bridge,npt,nfun,n,prn)
      implicit integer (a-z)
      real*8 hmat, kemat, norm
      real*8 keadd
      logical bridge, prn
      character*80 title
      dimension hmat(n,*), kemat(npt,npt), norm(nfun)
      common/io/inp, iout
      do 10 i=1,nfun
         do 20 j=1,i
            hmat(i,j) = hmat(i,j)  + norm(i)*norm(j)*kemat(i,j)
 20      continue   
 10   continue
c
c     we need to correct the last element if that element involves
c     a bridge function.
c
      if(bridge) then
         hmat(nfun,nfun) = hmat(nfun,nfun) + 
     #                             norm(nfun)*norm(nfun)*keadd 
      endif
      do 30 i=1,nfun
         do 40 j=1,i
            hmat(j,i) = hmat(i,j)
 40      continue
 30   continue   
      if(prn) then
         title='matrix elements in region i'
         call prntrm(title,hmat,nfun,nfun,n,n,iout)
      endif
      return
      end       



