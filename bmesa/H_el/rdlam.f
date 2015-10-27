*deck rdlam.f 
c***begin prologue     rdlam
c***date written       020303   (yymmdd)
c***revision date               (yymmdd)
c***keywords           
c***                   
c***author             schneider, b. i.(nsf)
c***source             hcoll
c***purpose            reading file for v_lambda and e_lambda
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       rdlam
      subroutine rdlam(v_lam,e_lam,nopt,n)
c
      implicit integer (a-z)
      real*8 v_lam, e_lam
      dimension v_lam(n,nopt), e_lam(nopt)
      common/io/inp, iout      
      read(inp,*)( e_lam(i),i=1,nopt)
      do 10 i=1,nopt
         read(inp,*) (v_lam(j,i),j=1,n)
 10   continue   
      return
      end


















