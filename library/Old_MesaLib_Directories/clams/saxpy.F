*deck @(#)saxpy.F	5.1  11/29/95
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
c***first executable statement  saxpy                                   
#ifdef LOCALBLAS
      if(n.le.0)return
      if (sa .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        sy(iy) = sy(iy) + sa*sx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sy(i) = sy(i) + sa*sx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        sy(i) = sy(i) + sa*sx(i)
        sy(i + 1) = sy(i + 1) + sa*sx(i + 1)
        sy(i + 2) = sy(i + 2) + sa*sx(i + 2)
        sy(i + 3) = sy(i + 3) + sa*sx(i + 3)
   50 continue
      return
      end
#endif LOCALBLAS
#ifdef MACHINEBLAS
      call daxpy(n,sa,sx,incx,sy,incy)
      return                                                            
      end                                                               
#endif MACHINEBLAS
