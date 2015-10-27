*deck honv.f
c***begin prologue     honv
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***purpose            multiply the hamiltonian matrix on a vector.
c***                   all special structure is accounted for assuming
c***                   one is using a representaton where r and theta are
c***                   separately diagonal.
c***                   the labelling of the vectors and implictly 
c***                   the matrix is done by fixing the index of the
c***                   basis function labelling the theta coordinate and
c***                   then letting the r index run from 1 to nr.
c***                   stated differently, it is as if we blocked the matrix
c***                   in theta. each of the submatrices would go over the 
c***                   theta index.
c***references         
c
c***routines called    
c***end prologue       honv
      subroutine honv(vecold,vecnew,tr,n,nvc)
      implicit integer (a-z)
      real*8 vecold, vecnew, tr
      character*80 title
      dimension vecold(n,nvc), vecnew(n,nvc), tr(n,n)
      common/io/inp, iout 
      call ebc(vecnew,tr,vecold,n,n,nvc)
c      title='old vectors'
c      call prntrm(title,vecold,n,nvc,n,nvc,iout)
c      title='new vectors'
c      call prntrm(title,vecnew,n,nvc,n,nvc,iout)
      return
      end       
