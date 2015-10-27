*deck @(#)mkfrbn.f	1.1 9/8/91
c***begin prologue     mkfrbn
c***date written       890528   (yymmdd)
c***revision date               (yymmdd)
c***keywords           bound-free
c***author             schneider, barry (lanl)
c***source             m6009
c***purpose            fill bound-free matrices
c*** 
c
c***references         none      
c
c***routines called    none
c***end prologue       mkfrbn
      subroutine mkfrbn(hpb,vfrb,vfib,vbfr,vbfi,matbv,ntchn)
      implicit integer (a-z)
      complex *16 hpb, vfib, vbfi
      real *8 vfrb, vbfr
      dimension hpb(ntchn,matbv), vfrb(ntchn,matbv), vfib(ntchn,matbv)
      dimension vbfr(matbv,ntchn), vbfi(matbv,ntchn)
      do 10 ch1=1,ntchn
         do 20 bfn=1,matbv
            vfrb(ch1,bfn)=imag(hpb(ch1,bfn))
            vfib(ch1,bfn)=hpb(ch1,bfn)
            vbfr(bfn,ch1)=vfrb(ch1,bfn)
            vbfi(bfn,ch1)=vfib(ch1,bfn)
   20    continue
   10 continue
      return
      end 
