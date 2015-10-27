*deck @(#)isamax.F	5.1  11/29/95
      function isamax(n,sx,incx)                                
c***begin prologue       isamax                                              
c***date written         791001   (yymmdd)                                    
c***revision date        860401   (yymmdd)                                    
c***category no.         d1a2                                                  
c***keywords             blas,linear algebra,maximum component,vector
c***author               lawson, c. l., (jpl)
c                        hanson, r. j., (snla) 
c                        kincaid, d. r., (u. of texas)
c
c***purpose              find the smallest index of an element of maximum
c                        magnitude of a vector.
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
c   isamax  smallest index (zero if n .le. 0)                           
c                                                                       
c     find smallest index of maximum magnitude of single precision sx.  
c     isamax =  first i, i = 1 to n, to minimize  abs(sx(1-incx+i*incx) 
c                                                                       
c***references           lawson c.l., hanson r.j., kincaid d.r., krogh f.t.,
c                 *basic linear algebra subprograms for fortran usage*, 
c                 algorithm no. 539, transactions on mathematical       
c                 software, volume 5, number 3, september 1979, 308-323 
c***routines called      (none)                                             
c***end prologue         isamax                                                
c                                                                       
      implicit real*8 (a-h,o-z)
      integer isamax
      real*8 sx(*)
c***first executable statement  isamax                                  


c
c if there's a local blas library (presumably optimized), use idamax, which
c should be the equivalent of this routine
#ifdef LOCALBLAS
      real*8 dmax
      isamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      isamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      dmax = dabs(sx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(sx(ix)).le.dmax) go to 5
         isamax = i
         dmax = dabs(sx(ix))
    5    ix = ix + incx
   10 continue
      return
c        code for increment equal to 1
c
   20 dmax = dabs(sx(1))
      do 30 i = 2,n
         if(dabs(sx(i)).le.dmax) go to 30
         isamax = i
         dmax = dabs(sx(i))
   30 continue
      return
      end
#endif LOCALBLAS
#ifdef MACHINEBLAS
      integer idamax
      isamax=idamax(n,sx,incx)
      return                                                            
      end                                                               
#endif MACHINEBLAS
