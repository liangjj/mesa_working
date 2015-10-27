*deck @(#)freefn.f
c***begin prologue     freefn
c***date written       920417   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           freefn, link 6080
c***author             schneider, barry (lanl)
c***source             m6080
c***purpose            free function, first and second derivative
c***                   to outgoing wave
c*** 
c
c***references       
c
c***routines called
c***end prologue       freefn
      subroutine freefn(cfn,ddcfn,pt,energy,npts)
      implicit integer (a-z)
      real *8  energy, k, pt, efac, fac, dfac, ddfac, sqk
      complex *16 cfn, ddcfn, ck, cfac
      dimension cfn(npts,2), ddcfn(npts,2), pt(npts)
      common/io/inp,iout
      k=sqrt(2.d0*energy)
      sqk=1.d0/sqrt(k)
      ck=cmplx(0.d0,k)
      do 10 i=1,npts
         cfn(i,1)=sqk*sin(k*pt(i))
         ddcfn(i,1)=-k*k*cfn(i,1)
         efac=exp(-pt(i))
         cfac=exp(ck*pt(i))
         fac=1.d0-efac
         dfac=efac
         ddfac=-efac 
         cfn(i,2)=sqk*fac*cfac
         ddcfn(i,2)=ddfac+2.d0*ck*dfac-k*k*fac
         ddcfn(i,2)=sqk*cfac*ddcfn(i,2)
   10 continue
      return
      end

















