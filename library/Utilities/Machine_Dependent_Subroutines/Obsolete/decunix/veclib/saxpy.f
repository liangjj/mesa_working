*deck @(#)saxpy.f	1.1  11/30/90
      subroutine saxpy(n,sa,sx,incx,sy,incy)                            
c***begin prologue     saxpy                                               
c***date written       791001   (yymmdd)                                    
c***revision date      820801   (yymmdd)                                    
c***category no.       d1a7                                                  
c***keywords           blas,linear algebra,triad,vector
c***author             lawson, c. l., (jpl) 
c                      hanson, r. j., (snla)
c                      kincaid, d. r., (u. of texas)
c                      krogh, f. t., (jpl)
c***purpose            s.p. computation y = a*x + y
c***description                                                         
c                                                                       
c                b l a s  subprogram                                    
c    description of parameters                                          
c                                                                       
c     --input--                                                         
c        n  number of elements in input vector(s)                       
c       sa  single precision scalar multiplier                          
c       sx  single precision vector with n elements                     
c     incx  storage spacing between elements of sx                      
c       sy  single precision vector with n elements                     
c     incy  storage spacing between elements of sy                      
c                                                                       
c     --output--                                                        
c       sy  single precision result (unchanged if n .le. 0)             
c                                                                       
c     overwrite single precision sy with single precision sa*sx +sy.    
c     for i = 0 to n-1, replace  sy(ly+i*incy) with sa*sx(lx+i*incx) +  
c       sy(ly+i*incy), where lx = 1 if incx .ge. 0, else lx = (-incx)*n 
c       and ly is defined in a similar way using incy.                  
c***references         lawson c.l., hanson r.j., kincaid d.r., krogh f.t.,
c                 *basic linear algebra subprograms for fortran usage*, 
c                 algorithm no. 539, transactions on mathematical       
c                 software, volume 5, number 3, september 1979, 308-323 
c***routines called    (none)                                             
c***end prologue       saxpy                                                 
c                                                                       
      implicit real*8 (a-h,o-z)
      real*8 sx(1),sy(1),sa                                             
      call daxpy(n,sa,sx,incx,sy,incy)                            
      return                                                            
      end                                                               
