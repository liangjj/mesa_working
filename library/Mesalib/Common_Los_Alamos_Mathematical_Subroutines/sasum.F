*deck @(#)sasum.F	5.1  11/29/95
      function sasum(n,sx,incx)
c***begin prologue  sasum                                               
c***date written   791001   (yymmdd)                                    
c***revision date  820801   (yymmdd)                                    
c***category no.  d1a3a                                                 
c***keywords  add,blas,linear algebra,magnitude,sum,vector              
c***author  lawson, c. l., (jpl)                                        
c           hanson, r. j., (snla)                                       
c           kincaid, d. r., (u. of texas)                               
c           krogh, f. t., (jpl)                                         
c***purpose  sum of magnitudes of s.p vector components                 
c***description                                                         
c                                                                       
c                b l a s  subprogram                                    
c    description of parameters                                          
c                                                                       
c     --input--                                                         
c        n  number of elements in input vector(s)                       
c       sx  single precision vector with n elements                     
c     incx  storage spacing between elements of sx                      
c                                                                       
c     --output--                                                        
c    sasum  single precision result (zero if n .le. 0)                  
c                                                                       
c     returns sum of magnitudes of single precision sx.                 
c     sasum = sum from 0 to n-1 of  abs(sx(1+i*incx))                   
c***references  lawson c.l., hanson r.j., kincaid d.r., krogh f.t.,     
c                 *basic linear algebra subprograms for fortran usage*, 
c                 algorithm no. 539, transactions on mathematical       
c                 software, volume 5, number 3, september 1979, 308-323 
c***routines called  (none)                                             
c***end prologue  sasum                                                 
c                                                                       
      real*8 sx(*),sasum
c***first executable statement  sasum                                   
#ifdef LOCALBLAS
      real*8 temp
      sasum = 0.0d0
      temp = 0.0d0
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        temp = temp + abs(sx(i))
   10 continue
      sasum = temp
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,6)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        temp = temp + abs(sx(i))
   30 continue
      if( n .lt. 6 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,6
        temp = temp + abs(sx(i)) + abs(sx(i + 1)) + abs(sx(i + 2))
     *  + abs(sx(i + 3)) + abs(sx(i + 4)) + abs(sx(i + 5))
   50 continue
   60 sasum = temp
      return
      end
#endif LOCALBLAS
#ifdef MACHINEBLAS
      real*8 dasum
      sasum=dasum(n,sx,incx)
      return                                                            
      end                                                               
#endif MACHINEBLAS
