*deck egnpse
c***begin prologue     egnpse
c***date written       940213   (yymmdd)
c***revision date               (yymmdd)
c***keywords           k-matrix, eigenphase
c***author             schneider, barry (nsf)
c***source             
c***purpose            eigenphases from k-matrix
c***description        
c***references       
c
c***routines called
c***end prologue       egnpse
      subroutine egnpse(kmat,eig,psesum,dum,nc,nopen,energy,prnt)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 kmat, eig, dum, psesum, energy
      logical prnt
      dimension kmat(nc,nc), eig(nc), dum(nc)
      call tred2(nc,nopen,kmat,eig,dum,kmat)
      call tql2(nc,nopen,eig,dum,kmat,ierr)
      psesum=0.d0
      do 10 i=1,nopen
         eig(i)=atan(eig(i))
         psesum=psesum+eig(i)
   10 continue
      write (iout,1) energy
      do 20 i=1,nopen
         write(iout,2) i, eig(i)
         if (prnt) then
             write (iout,3) ( kmat(j,i), j=1,nopen)
         endif
 20   continue   
      write(iout,4) psesum
 1    format (/,'phase shift analysis at energy = ',e15.8)
 2    format(/,'eigenphase ',i3,' = ',e15.8)   
 3    format(/,'eigenphase vector',(/,5e15.8))
 4    format(/,'eigenphase sum = ',e15.8)      
      return
      end
