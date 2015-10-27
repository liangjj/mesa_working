*deck vmat.f
c***begin prologue     vmat
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            potential matrix in dvr representation
c***                   
c***description        
c***                     
c***references         
c
c***routines called    
c***end prologue       vmat
      subroutine vmat(v,eigr,n,pottyp)
      implicit integer (a-z)
      real*8 eigr
      complex*16 v, eye
      character*(*) pottyp
      dimension v(n,n), eigr(n) 
      common/io/inp, iout
      data eye /(0.d0,1.d0) / 
      call czero(v,n*n)
      if(pottyp.eq.'harmonic-oscillator') then
         do 10 i=1,n
            v(i,i) = v(i,i) + .5d0*eigr(i)*eigr(i)
 10      continue
      elseif(pottyp.eq.'exponential') then
         do 20 i=1,n
            v(i,i) = v(i,i) - exp(-eigr(i))
 20      continue
      elseif(pottyp.eq.'well') then
         do 30 i=1,n
            v(i,i) = v(i,i) - 1.0d0
 30      continue
      elseif(pottyp.eq.'complex-well') then
         do 40 i=1,n
            v(i,i) = v(i,i) -1.d0 - eye*.2d0
 40      continue            
      endif
      return
      end       
