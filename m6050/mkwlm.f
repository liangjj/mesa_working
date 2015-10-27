*deck @(#)mkwlm.f	
c***begin prologue     mkwlm
c***date written       920402   (yymmdd)
c***revision date      
c***keywords           mkwlm, link 6050
c***authors            Schneider,B (NSF)
c***                   
c***source             m6050
c***purpose            scale the vlamda
c***references       
c***routines called    
c***end prologue       mkwlm
      subroutine mkwlm(vlamda,grid,npt,nlam)
      implicit real *8 (a-h,o-z)
      real *8 grid, pi, pifac, rsq, rfac, pre
      complex *16 vlamda
      dimension vlamda(npt,nlam), grid(4,npt)
      data pi /3.14159265358979323846d+00/
      pifac=1.d0/sqrt(4.d0*pi)
      do 10 i=1,npt
         rsq=grid(1,i)*grid(1,i)+
     1                  grid(2,i)*grid(2,i)+
     2                            grid(3,i)*grid(3,i)
         rfac=1.d0/sqrt(rsq)
         pre=pifac*rfac
         do 20 lam=1,nlam
            vlamda(i,lam)=vlamda(i,lam)*pre
   20    continue
   10 continue     
      return
      end














