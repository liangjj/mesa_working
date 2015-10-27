*deck @(#)mkfree.f	1.1 9/8/91
c***begin prologue     mkfree
c***date written       890528   (yymmdd)
c***revision date               (yymmdd)
c***keywords           free-free
c***author             schneider, barry (lanl)
c***source             m6009
c***purpose            fill free-free matrices
c*** 
c
c***references         none      
c
c***routines called    none
c***end prologue       mkfree
      subroutine mkfree(hpp,hpm,vrr,vii,vri,vir,ntchn)
      implicit integer (a-z)
      complex *16 hpp, hpm, vii, vri, vir
      real *8 vrr
      dimension hpp(ntchn,ntchn), hpm(ntchn,ntchn), vrr(ntchn,ntchn)
      dimension vii(ntchn,ntchn), vri(ntchn,ntchn), vir(ntchn,ntchn)
      do 10 ch1=1,ntchn
         do 20 ch2=1,ch1
            vii(ch1,ch2)=hpp(ch1,ch2)
            vii(ch2,ch1)=vii(ch1,ch2)
            vrr(ch1,ch2)=imag(hpm(ch1,ch2))
            vrr(ch2,ch1)=vrr(ch1,ch2)
            vir(ch1,ch2)=hpm(ch1,ch2)
            vri(ch2,ch1)=vir(ch1,ch2)
            vir(ch2,ch1)=hpm(ch2,ch1)
            vri(ch1,ch2)=vir(ch2,ch1)
   20    continue
   10 continue
      return
      end  
