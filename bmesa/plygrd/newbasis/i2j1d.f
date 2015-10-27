*deck i2j1d.f
c***begin prologue     i2j1d
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            interpolate a 1d vector from grid i
c***                   to a grid j.  the storage in the vector
c***                   corresponds to a loop structure looking like
c***                                  ........
c***                   
c***references         
c
c***routines called    
c***end prologue       i2j1d
      subroutine i2j1d(p1ji,vi,vj,n1j,n1i,nvc)
      implicit integer (a-z)
      real*8 p1ji, vi, vj
      dimension p1ji(n1j,n1i)
      dimension vi(n1i,nvc), vj(n1j,nvc)
      common/io/inp, iout
      call ebc(vj,p1ji,vi,n1j,n1i,nvc)
      return
      end       

