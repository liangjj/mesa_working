*deck richdn.f
c***begin prologue     richdn
c***date written       950704   (yymmdd)
c***revision date               (yymmdd)
c***keywords           finite difference, band, eigenvalues
c***author             schneider, barry(nsf)
c***source             
c***purpose            richardson extrapolation of eigenvalues/eigenvectors 
c***                   of tridiagonal matrix.
c***
c***references
c
c***routines called
c***end prologue      richdn
      subroutine richdn(eig1,eig2,vec1,vec2,stp,n1,n2,nroots)
      implicit integer (a-z)
      dimension  eig1(n1), eig2(n2), vec1(n1,*), vec2(n2,*), stp(2) 
      real*8 eig1, eig2, vec1, vec2, stp, h1sq, h2sq
      common /io/ inp, iout
      h1sq=stp(1)*stp(1)
      h2sq=stp(2)*stp(2)
      do 10 i=2,nroots+1
         eig1(i)=h2sq*eig1(i)/(h2sq-h1sq) + eig2(i)*h1sq/(h1sq-h2sq)
 10   continue
      write(iout,1) (eig1(i),i=2,nroots+1)   
      return
 1    format(/,'the extrapolated eigenvalues',/,(5e15.8))
      end

