*deck calfun.f
c***begin prologue     calfun
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***                   
c***                                                          
c***references         
c
c***routines called    
c***end prologue       calfun
      subroutine calfun(d,c,over,nf,nc)
      implicit integer (a-z)
      real*8 d, c, over
      dimension d(nc), c(nf), over(nc,nf)
      common/io/inp, iout
      call ebc(d,over,c,nc,nf,1)
      return
      end
