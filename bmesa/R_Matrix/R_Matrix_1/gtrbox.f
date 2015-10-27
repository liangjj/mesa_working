*deck gtrbox.f
c***begin prologue     gtrbox
c***date written       940213   (yymmdd)
c***revision date               (yymmdd)
c***keywords           external solutions
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***description          
c***references       
c
c***routines called
c***end prologue       gtrbox
      subroutine gtrbox(pgrid,rbox,n)
      implicit integer (a-z)
      integer*8 pgrid
      real*8 rbox, x
      common /io/ inp, iout
      pointer (pgrid,x(1))
      h=1
      vphy=h+n*n
      h0=vphy+n
      srf=h0+n*n
      pt0=srf+2
      last=pt0+n
      rbox=x(last)
      return
      end
