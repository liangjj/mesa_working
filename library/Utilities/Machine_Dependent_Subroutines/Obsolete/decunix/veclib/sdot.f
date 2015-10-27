*deck @(#)sdot.f	1.1  11/30/90
      function sdot(n,sx,incx,sy,incy)                           
c***begin prologue  sdot                                                
c***date written   791001   (yymmdd)                                    
c***revision date  820801   (yymmdd)                                    
c***category no.  d1a4                                                  
c***keywords  blas,inner product,linear algebra,vector                  
c***author  lawson, c. l., (jpl)                                        
c           hanson, r. j., (snla)                                       
c           kincaid, d. r., (u. of texas)                               
c           krogh, f. t., (jpl)                                         
c***purpose  s.p. inner product of s.p. vectors                         
c***description                                                         
c                                                                       
c                b l a s  subprogram                                    
c    description of parameters                                          
c                                                                       
c     --input--                                                         
c        n  number of elements in input vector(s)                       
c       sx  single precision vector with n elements                     
c     incx  storage spacing between elements of sx                      
c       sy  single precision vector with n elements                     
c     incy  storage spacing between elements of sy                      
c                                                                       
c     --output--                                                        
c     sdot  single precision dot product (zero if n .le. 0)             
c                                                                       
c     returns the dot product of single precision sx and sy.            
c     sdot = sum for i = 0 to n-1 of  sx(lx+i*incx) * sy(ly+i*incy),    
c     where lx = 1 if incx .ge. 0, else lx = (-incx)*n, and ly is       
c     defined in a similar way using incy.                              
c***references  lawson c.l., hanson r.j., kincaid d.r., krogh f.t.,     
c                 *basic linear algebra subprograms for fortran usage*, 
c                 algorithm no. 539, transactions on mathematical       
c                 software, volume 5, number 3, september 1979, 308-323 
c***routines called  (none)                                             
c***end prologue  sdot                                                  
c                                                                       
      implicit real*8 (a-h,o-z)
      real*8 sdot, dsdot
      real*8 sx(1),sy(1)                                                
c***first executable statement  sdot                             
      sdot= dsdot(n,sx,incx,sy,incy)                                  
      return                                                            
      end                                                               
