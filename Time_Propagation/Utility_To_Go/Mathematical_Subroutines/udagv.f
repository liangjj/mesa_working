*deck udagv.f
c***begin prologue     udagv
c***date written       970430   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            perform a unitary transformation on a vector
c                           +
c***description        v = u v 
c***
c                                v = v(n,m) = input and output
c                                u = u(n,n) = unitary transformation matrix
c                                s = s(n,m) = scratch matrix
c***references         
c
c***routines called    
c***end prologue       udagv
      subroutine udagv(v,u,s,n,m)
      implicit integer (a-z)
      real*8 v, u, s
      dimension v(n,m), u(n,n), s(n,m)
      common/io/inp, iout 
      call ebtc(s,u,v,n,n,m)
      call copy(s,v,n*m)
      return
      end       
