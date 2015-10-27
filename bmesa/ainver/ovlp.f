*deck ovlp.f
c***begin prologue     ovlp
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
c***end prologue       ovlp
      subroutine ovlp(over,pfine,pcorse,wt,nf,nc)
      implicit integer (a-z)
      real*8 over, pfine, pcorse, wt
      dimension over(nc,nf), pfine(nf,nf), pcorse(nf,nc), q(nf), wt(nf)
      common/io/inp, iout
      do 10 i=1,nc
         do 20 j=1,nf
            over(i,j) = wt(j)*pfine(j,j)*pcorse(j,i)
 20      continue   
 10   continue   
      return
      end
