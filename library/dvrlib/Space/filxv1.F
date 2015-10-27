*deck filxv1.f
c***begin prologue     filxv1
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           potential
c***author             schneider, barry (nsf)
c***source             
c***purpose            fill grid and potential matrix
c***                   
c***description        fill the physical grid and one body potential 
c***                   matrix arrays in the dvr representation.
c***                   
c***references         
c
c***routines called    
c***end prologue       filxv1
      subroutine filxv1(pham,px,pv,n,word)
      implicit integer (a-z)
#ifdef DECPOINTER
      integer*8 pham, px, pv
#endif DECPOINTER
#ifdef SGIPOINTER
      integer*4 pham, px, pv
#endif SGIPOINTER
      real*8 ham, x, v
      dimension word(2)
      common/io/inp, iout
      pointer (pham,ham(1))
      pointer (px,x(1))
      pointer (pv,v(1))
c
c     locate the needed arrays
c
      h=1
      vphy=h+n*n
      h0=vphy+n 
      srf=h0+n*n
      q0=srf+2
      q1=q0+1
      need=wpadti(1+n)
      call getmem(need,px,word(1),'x',0)
      need=wpadti(1+n)
      call getmem(need,pv,word(2),'v',0)
      call copy(ham(q1),x,n)
      call copy(ham(vphy),v,n)
      return
      end
