*deck scaprd.f
c***begin prologue     scaprd
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           scalar product for generalized polynomials
c***author             schneider, barry (nsf)
c***source
c***purpose            
c***                   
c***                   
c***description      
c***                 
c***                                                                       
c***                                                          
c***references         
c
c***routines called    
c***end prologue       scaprd
      function scaprd (va,vb,wt,rwt,n)
      implicit integer (a-z)
      real*8 scaprd, va, vb, wt, rwt
      dimension va(n), vb(n), wt(n), rwt(n)
      common/io/inp, iout 
      scaprd=0.d0
      do 10 i=1,n
         scaprd=scaprd+va(i)*wt(i)*rwt(i)*vb(i)
 10   continue   
      return
      end       
