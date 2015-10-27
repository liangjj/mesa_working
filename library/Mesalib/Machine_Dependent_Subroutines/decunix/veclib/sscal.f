*deck @(#)sscal.f	1.1  11/30/90
      subroutine sscal(n,sa,sx,incx)                                    
c***begin prologue  sscal                                               
c***date written   791001   (yymmdd)                                    
c***revision date  820801   (yymmdd)                                    
c***category no.  d1a6                                                  
c***keywords  blas,linear algebra,scale,vector                          
c***author  lawson, c. l., (jpl)                                        
c           hanson, r. j., (snla)                                       
c           kincaid, d. r., (u. of texas)                               
c           krogh, f. t., (jpl)                                         
c***purpose  s.p. vector scale x = a*x                                  
c***description                                                         
c                                                                       
c                b l a s  subprogram                                    
c    description of parameters                                          
c                                                                       
c     --input--                                                         
c        n  number of elements in input vector(s)                       
c       sa  single precision scale factor                               
c       sx  single precision vector with n elements                     
c     incx  storage spacing between elements of sx                      
c                                                                       
c     --output--                                                        
c       sx  single precision result (unchanged if n .le. 0)             
c                                                                       
c     replace single precision sx by single precision sa*sx.            
c     for i = 0 to n-1, replace sx(1+i*incx) with  sa * sx(1+i*incx)    
c***references  lawson c.l., hanson r.j., kincaid d.r., krogh f.t.,     
c                 *basic linear algebra subprograms for fortran usage*, 
c                 algorithm no. 539, transactions on mathematical       
c                 software, volume 5, number 3, september 1979, 308-323 
c***routines called  (none)                                             
c***end prologue  sscal                                                 
c                                                                       
      implicit real*8 (a-h,o-z)
      real*8 sa,sx(1)                                                   
c***first executable statement  sscal                                   
      call dscal(n,sa,sx,incx)                                    
      return                                                            
      end                                                               
