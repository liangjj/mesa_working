*deck pkmat.f
c***begin prologue     pkmat
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           repack even-odd parity matrices
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***                   
c***references         
c
c***routines called    
c***end prologue       pkmat
      subroutine pkmat(mat,tmp,n,n1,n2,type,dir)
      implicit integer (a-z)
      real*8 mat, tmp
      character*(*) type, dir
      dimension mat(n,*), tmp(*)
      common/io/inp, iout
      count=0
      if(dir.eq.'matrix') then
         if(type.eq.'even') then
            do 10 i=1,n1
               do 20 j=1,n1
                 count=count+1
                 tmp(count)=mat(j,i)
 20            continue
 10         continue
            do 30 i=n1+1,n
               do 40 j=n1+1,n
                  count=count+1
                  tmp(count)=mat(j,i)
 40            continue
 30         continue   
         elseif(type.eq.'odd') then
            do 50 i=n1+1,n
               do 60 j=n1+1,n
                  count=count+1
                  tmp(count)=mat(j,i)
 60            continue
 50         continue
            do 70 i=1,n1
               do 80 j=1,n1
                  count=count+1
                  tmp(count)=mat(j,i)
 80            continue
 70         continue
         else   
            call lnkerr('error in parity type')
         endif
         call copy(tmp,mat,n*n)
      elseif(dir.eq.'vector') then
         if(type.eq.'even') then
            do 100 i=1,n1
               do 110 j=1,n
                 count=count+1
                 tmp(count)=mat(j,i)
 110          continue
 100        continue
            do 120 i=n1+1,n
               do 130 j=1,n
                  count=count+1
                  tmp(count)=mat(j,i)
 130           continue
 120        continue   
         elseif(type.eq.'odd') then
            do 200 i=n1+1,n
               do 210 j=1,n
                  count=count+1
                  tmp(count)=mat(j,i)
 210           continue
 200        continue
            do 220 i=1,n1
               do 230 j=1,n
                  count=count+1
                  tmp(count)=mat(j,i)
 230           continue
 220        continue
         else   
            call lnkerr('error in parity type')
         endif
         call copy(tmp,mat,n*n)
      endif
      return
      end       



