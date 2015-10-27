*deck addend.f
c***begin prologue     addend
c***date written       950704   (yymmdd)
c***revision date               (yymmdd)
c***keywords           finite difference, band, eigenvalues
c***author             schneider, barry(nsf)
c***source             
c***purpose            add the endpoints of the vector based on the
c***                   boundary conditions.
c***routines called
c***end prologue      addend
      subroutine addend(vec,nroots,n,bc)
      implicit integer (a-z)
      dimension vec(n,nroots)
      real*8 vec
      character*(*) bc
      common /io/ inp, iout
      do 10 i=1,nroots
         vec(1,i)=0.d0
         if(bc.eq.'zero-function') then
            vec(n,i)=0.d0
          else
            vec(n,i) = ( 4.d0*vec(n-1,i)-vec(n-2,i) )/3.d0
c             vec(n,i)=vec(n-1,i)
          endif
 10   continue   
      return
      end

