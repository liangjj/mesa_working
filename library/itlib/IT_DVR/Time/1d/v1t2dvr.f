*deck v1t2dvr.f
c***begin prologue     v1t2dvr
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            transform a 1-D vector from the dvr to  
c***                   the H0 representation.
c***                                      
c***references         
c
c***routines called    
c***end prologue       v1t2dvr
      subroutine v1t2dvr(vecin,vecout,u1,n1,nt,nc,nvc)
      implicit integer (a-z)
      real*8 ux, vecout, vecin
      dimension u1(n1,n1), vecin(n1,nt,nc,2,nvc), vecout(n1,nt,nc,2,nvc)
      common/io/inp, iout
      call ebc(vecout,u1,vecin,n1,n1,nt*nc*2*nvc)
      return
      end       
