*deck honc.f
c***begin prologue     honc
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
c***end prologue       honc
      subroutine honc(vecold,vecnew,tr,tth,eigr,eigth,small,nr,nth,prn)
      implicit integer (a-z)
      real*8 vecold, vecnew, tr, tth, eigr, eigth, small
      logical prn
      character*80 title
      dimension vecold(nr*nth), vecnew(nr*nth), tr(nr,nr), tth(nth,nth)
      dimension eigr(nr), eigth(nth), small(nr,nr)
      common/io/inp, iout 
      ntot=nr*nth
      call rzero(vecnew,ntot)
      counti=0
      do 10 i=1,nth
         counti=counti+1
         countj=0
         do 20 j=1,i
c           fill the h(i,j) block with the radial matrix.
            call filham(small,tr,tth(i,j),eigr,i,j,nr,nr)
c           the next business is to avoid double counting 
c                      on diagonal elements
            if(i.eq.j) then
               call vscale(small,small,.5d0,nr*nr)
            endif
            countj=countj+1
c           this is where the double counting occurs when i = j
            call apbc(vecnew(counti),small,vecold(countj),nr,nr,1)
            call apbc(vecnew(countj),small,vecold(counti),nr,nr,1)
            countj=countj+nr
 20      continue                           
         counti=counti+nr 
 10   continue   
      if (prn) then
          title='next vector'
          call prntrm(title,vecnew,ntot,1,ntot,1,iout)
      endif
      return
      end       
