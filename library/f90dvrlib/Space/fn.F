*deck fn.f
c***begin prologue     fn
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***                   
c***references         
c
c***routines called    
c***end prologue       fn
      subroutine fn(pt,f,n,ftyp,typ)
      implicit integer (a-z)
      real*8 pt, f
      character*(*) ftyp, typ
      dimension pt(n), f(n)
      common/io/inp, iout
      if(ftyp.eq.'negative-exponential') then
         if(typ.eq.'function') then
            do 10 i=1,n
               f(i)=exp(-pt(i))
 10         continue   
         elseif(typ.eq.'first-derivative') then
            do 20 i=1,n
               f(i)=- exp(-pt(i))
 20         continue   
         elseif(typ.eq.'second-derivative') then
            do 30 i=1,n
               f(i)=exp(-pt(i))
 30         continue   
         endif
      elseif(ftyp.eq.'positive-exponential') then
         if(typ.eq.'function') then
            do 40 i=1,n
               f(i)=exp(pt(i))
 40         continue   
         elseif(typ.eq.'first-derivative') then
            do 50 i=1,n
               f(i)= exp(pt(i))
 50         continue   
         elseif(typ.eq.'second-derivative') then
            do 60 i=1,n
               f(i)=exp(pt(i))
 60         continue   
         endif
      elseif(ftyp.eq.'sine') then
         if(typ.eq.'function') then
            do 70 i=1,n
               f(i)=sin(pt(i))
 70         continue   
         elseif(typ.eq.'first-derivative') then
            do 80 i=1,n
               f(i)= cos(pt(i))
 80         continue   
         elseif(typ.eq.'second-derivative') then
            do 90 i=1,n
               f(i)= - sin(pt(i))
 90         continue   
         endif
      elseif(ftyp.eq.'cosine') then
         if(typ.eq.'function') then
            do 100 i=1,n
               f(i)=cos(pt(i))
 100        continue   
         elseif(typ.eq.'first-derivative') then
            do 110 i=1,n
               f(i)= - sin(pt(i))
 110        continue   
         elseif(typ.eq.'second-derivative') then
            do 120 i=1,n
               f(i)= - cos(pt(i))
 120        continue   
         endif
      elseif(ftyp.eq.'x**2') then
         if(typ.eq.'function') then
            do 200 i=1,n
               f(i)=pt(i)*pt(i)
 200        continue   
         elseif(typ.eq.'first-derivative') then
            do 210 i=1,n
               f(i)= 2.d0*pt(i)
 210        continue   
         elseif(typ.eq.'second-derivative') then
            do 220 i=1,n
               f(i)= 2.d0
 220        continue   
         endif
      elseif(ftyp.eq.'x') then
         if(typ.eq.'function') then
            do 300 i=1,n
               f(i)=pt(i)
 300        continue   
         elseif(typ.eq.'first-derivative') then
            do 310 i=1,n
               f(i)= 1.d0
 310        continue   
         elseif(typ.eq.'second-derivative') then
            do 320 i=1,n
               f(i)= 0.d0
 320        continue   
         endif
      endif
      return
      end       



