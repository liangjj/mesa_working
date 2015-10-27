*deck %W%  %G%
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
      real*8 sx(*),sasum,dasum
      sasum=dasum(n,sx,incx)
      return                                                            
      end                                                               
