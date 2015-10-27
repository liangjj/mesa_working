*deck %W%  %G%
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
      real*8 sx(*),smax,xmag 
c***first executable statement  isamax                                  
      isamax = 0                                                        
      if(n.le.0) return                                                 
      isamax = 1                                                        
      if(n.le.1)return                                                  
      if(incx.eq.1)goto 20                                              
c                                                                       
c        code for increments not equal to 1.                            
c                                                                       
      smax = abs(sx(1))                                                 
      ns = n*incx                                                       
      ii = 1                                                            
          do 10 i=1,ns,incx                                             
          xmag = abs(sx(i))                                             
          if(xmag.le.smax) go to 5                                      
          isamax = ii                                                   
          smax = xmag                                                   
    5     ii = ii + 1                                                   
   10     continue                                                      
      return                                                            
c                                                                       
c        code for increments equal to 1.                                
c                                                                       
   20 smax = abs(sx(1))                                                 
      do 30 i = 2,n                                                     
         xmag = abs(sx(i))                                              
         if(xmag.le.smax) go to 30                                      
         isamax = i                                                     
         smax = xmag                                                    
   30 continue                                                          
      return                                                            
      end                                                               
