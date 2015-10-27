*deck prepfn.f
c***begin prologue     prepfn
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***purpose            prepare orthonormal basis set from polynomials
c***                   
c***                                                          
c***references         
c
c***routines called    
c***end prologue       prepfn
      subroutine prepfn(pppji,ppqi,ni,nj)
      implicit integer (a-z)
      real*8 pti, pji, fac, sumwt
      common/io/inp, iout
      pointer (ppqi,pti(1))
      pointer (pppji,pji(1))
      sumwt=0.d0
      q=1
      wt=q+ni
      p=1
      dp=p+ni*nj
      ddp=dp+ni*nj
      do 10 i=1,ni
         sumwt = sumwt + pti(wt)
         fac=1.d0/sqrt(pti(wt))
         wt=wt+1 
         call smul(pji(p),pji(p),fac,nj)
         call smul(pji(dp),pji(dp),fac,nj)
         call smul(pji(ddp),pji(ddp),fac,nj)
         p=p+nj
         dp=dp+nj
         ddp=ddp+nj
 10   continue
      write(iout,1) sumwt
      return
 1    format(/,5x,'sum of the weights = ',e15.8)      
      end       
