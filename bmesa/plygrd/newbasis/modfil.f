*deck modfil.f
c***begin prologue     modfil
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            fill arrays of lobatto polynomials.
c***                   
c***references         
c
c***routines called    
c***end prologue       modfil
      subroutine modfil(ppoly,ngrid,npt)
      implicit integer (a-z)
      real*8 poly
      dimension npt(ngrid)
      common/io/inp, iout
      pointer (ppoly,poly(1))
      cntply=1
      do 10 grdi=1,ngrid
         do 20 grdj=1,ngrid
            p=cntply
            dp=p+npt(grdi)*npt(grdj)
            ddp=dp+npt(grdi)*npt(grdj)
            pn=ddp+npt(grdi)*npt(grdj)
            dpn=pn+npt(grdi)*npt(grdj) 
            ddpn=dpn+npt(grdi)*npt(grdj)
            cntply=ddpn+npt(grdi)*npt(grdj)
            call copy(poly(p),poly(pn),npt(grdi)*npt(grdj))
            call copy(poly(dp),poly(dpn),npt(grdi)*npt(grdj))
            call copy(poly(ddp),poly(ddpn),npt(grdi)*npt(grdj))
 20      continue   
 10   continue   
      return
      end       
