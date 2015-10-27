*deck @(#)mskone.f	1.1 9/8/91
c***begin prologue     m7001
c***date written       920525   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m7001, link 7001, spline
c***author             schneider, b. (nsf)
c***source             m7001
c***purpose            symmetry zeros of one electron integrals
c***                   due to angular momentum considerations.
c***
c***description        a matrix of zeros and ones is constructed
c***                   which is later used in a multiply operation 
c***                   to zero the unwanted integrals.
c***                 
c***references       
c
c***routines called    rzero(mylib)
c***end prologue       m7001
      subroutine mskone(zromat,la,ma,na,lb,mb,nb)
      implicit integer (a-z)
      real *8 zromat
      dimension zromat(na,nb), la(na), ma(na), lb(nb), mb(nb)
      common /io/ inp, iout
      call rzero(zromat,na*nb)
      do 10 i=1,na
         do 20 j=1,nb
            if (la(i).eq.lb(j)) then
                if(ma(i).eq.mb(j)) then
                   zromat(i,j)=1.d0
                endif
            endif
   20    continue         
   10 continue         
      return
      end


