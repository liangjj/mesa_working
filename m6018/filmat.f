*deck @(#)filmat.f	1.1 9/8/91
c***begin prologue     filmat
c***date written       930221   (yymmdd)
c***revision date               (yymmdd)
c***keywords           fill, matrix
c***author             schneider barry (nsf)
c***source             math
c***purpose            fill a matrix using a list
c***routines called
c***end prologue       filmat  
      subroutine filmat(matin,matout,m,ma,n,na,list)
      implicit integer (a-z)
      real*8 matin, matout
      dimension matin(ma,m), matout(na,n), list(n)
      common/io/ inp, iout
      do 10 i=1,n
         do 20 j=1,n
            matout(i,j) = matin(list(i),list(j)) 
   20    continue
   10 continue
      return
      end
