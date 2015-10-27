*deck %W%  %G%
      subroutine sswap(n,sx,incx,sy,incy)                               
c***begin prologue  sswap                                               
c***date written   791001   (yymmdd)                                    
c***revision date  820801   (yymmdd)                                    
c***category no.  d1a5                                                  
c***keywords  blas,interchange,linear algebra,vector                    
c***author  lawson, c. l., (jpl)                                        
c           hanson, r. j., (snla)                                       
c           kincaid, d. r., (u. of texas)                               
c           krogh, f. t., (jpl)                                         
c***purpose  interchange s.p vectors                                    
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
c       sx  input vector sy (unchanged if n .le. 0)                     
c       sy  input vector sx (unchanged if n .le. 0)                     
c                                                                       
c     interchange single precision sx and single precision sy.          
c     for i = 0 to n-1, interchange  sx(lx+i*incx) and sy(ly+i*incy),   
c     where lx = 1 if incx .ge. 0, else lx = (-incx)*n, and ly is       
c     defined in a similar way using incy.                              
c***references  lawson c.l., hanson r.j., kincaid d.r., krogh f.t.,     
c                 *basic linear algebra subprograms for fortran usage*, 
c                 algorithm no. 539, transactions on mathematical       
c                 software, volume 5, number 3, september 1979, 308-323 
c***routines called  (none)                                             
c***end prologue  sswap                                                 
c                                                                       
      implicit real*8 (a-h,o-z)
      real*8 sx(*),sy(*),stemp1,stemp2,stemp3                           42
      call dswap(n,sx,incx,sy,incy)                               
      return                                                            
      end                                                               
