      SUBROUTINE  MONDER (N, X, F, D, H, SLOPE)                         
      implicit real*8 (a-h,o-z)
      integer  n                                                        
      real*8  x(n), f(n), d(n), slope(n), h(n)                            
c                                                                       
c-----------------------------------------------------------------------
c                                                                       
c     monder uses the fritsch-carlson forumlas to set derivative        
c     values for a piecewise cubic interpolant to the data (x, f)       
c     so that the interpolant is monotone on any subinterval on         
c     which the data are monotone.                                      
c                                                                       
c     this version uses:                                                
c       1. three-point difference formulas to initialize derivatives    
c          (including endpoints).                                       
c       2. region s(2).                                                
c       3. algorithm a for moving a point into region.                  
c       4. any negative alpha or beta (indicating a change in mono-     
c          tonicity of the data) is set to zero to insure the strict    
c          piecewise monotonicity of the interpolant.                   
c                                                                       
c     subroutine pwcfev may be used to evaluate the resulting        
c     piecewise cubic function.                                        
c                                                                       
c     reference:  f. n. fritsch and r. e. carlson, monotone piecewise  
c       cubic interpolation, lawrence livermore laboratory report     
c       ucrl-82453 (january 1979).                                      
c                                                                       
c-----------------------------------------------------------------------
c                                                                    
c     on input:                                                         
c        n      is the number of data points.                           
c               restriction:  n.ge.4  (not checked).                    
c        x      is the array of independent variable values.            
c               restriction:  x must be strictly increasing, that is    
c                  x(i) .lt. x(i+1), i=1(1)n-1  (not checked).          
c        f      is the array of dependent variable values.              
c                                                                       
c     on output:                                                        
c        d      will be set to the desired derivative values.           
c        h      will be the array of interval lengths,                  
c                  h(i) = x(i+1) - x(i), i=1(1)n-1.                     
c        slope  will be the array of slopes of chords,                  
c                  slope(i) = (f(i+1) - f(i))/h(i), i=1(1)n-1.          
c                                                                       
c     note:  arrays h and slope are no longer needed after the call to  
c            monder.                                                    
c                                                                       
c     fortran intrinsics used:  abs.                                    
c                                                                       
c-----------------------------------------------------------------------
c                                                                     
c     algorithm by:  f. n. fritsch, lawrence livermore laboratory, and  
c                    r. e. carlson, grove city college, pa.             
c     programmed by: f. n. fritsch.                                     
c     date last changed:  18 july 1979  (fnf)                          
c                                                                       
c     change record:                                                    
c        78-12-07   minor cosmetic changes to get ready for library.    
c           1. removed argument icount.                         
c                   2. changed argument y to f (to be consistent with   
c                      pwcfev).                                         
c           1. changed treatment of interval adjacent to change 
c                      in monotonicity of data (see item 4, above).     
c                   2. minor additions to comment section.              
c           changed to use region s(2) instead of s(3).         
c                                                                       
c-----------------------------------------------------------------------
c                                                                       
c        local declarations.                                            
c                                                                       
      integer  i, nless1                                                
      real*8  alpha, beta, delta, fuzz, tau                               
      data  fuzz /1.0e-14/                                              
c                                                                       
c        initialize.                                                    
c                                                                       
      nless1 = n - 1                                                    
c                                                                       
c        compute interval lengths and slopes.                           
c                                                                       
      do 10  i = 1, nless1                                              
         h(i) = x(i+1) - x(i)                                           
         slope(i) = (f(i+1) - f(i))/h(i)                                
   10 continue                                                          
c                                                                       
c        initialize d(1) via non-centered three-point formula.          
c                                                                       
      d(1) = ((h(1)+h(1)+h(2))*slope(1) - h(1)*slope(2))/(h(1)+h(2))    
      if (d(1)*slope(1) .lt. 0.)  d(1) = 0.                             
c                                                                       
c        cycle through all intervals.                                   
c                                                                       
      do 50  i = 1, nless1                                              
         if (i .lt. nless1)  go to 20                                   
c                                                                       
c        special case of right endpoint.                                
            d(n) = ((h(n-1)+h(n-1)+h(n-2))*slope(n-1)                   
     *                           - h(n-1) *slope(n-2))/(h(n-2)+h(n-1))  
            if (d(n)*slope(n-1) .lt. 0.)  d(n) = 0.                     
            go to 25                                                    
   20    continue                                                       
c                                                                       
c        use three-point formula to initialize right-hand               
c        derivative for interval (x(i), x(i+1)) .                       
c                                                                       
         d(i+1) = (h(i+1)*slope(i) + h(i)*slope(i+1))/(h(i)+h(i+1))     
c                                                                       
   25    continue                                                       
c                                                                       
c        adjust d(i) and/or d(i+1), if necessary to insure monotonicity 
c        on interval (x(i), x(i+1)) .                                   
c                                                                       
c           take care of flat data.                                     
c                                                                       
         if (abs(slope(i)) .gt. fuzz)  go to 30                         
            alpha = 0.                                                  
            beta = 0.                                                   
            go to 45                                                    
   30    continue                                                       
c                                                                       
c           compute scaled derivatives.                                 
c                                                                       
         alpha = d(i) / slope(i)                                        
         beta = d(i+1) / slope(i)                                       
c                                                                       
c           take care of nonmonotone data.                              
c                                                                       
c        assertion:  if either of the following tests is satisfied,     
c                    (alpha,beta) is not in first quadrant, which       
c                    means that slope changes sign at one or both ends  
c                    of interval.                                       
         if (alpha .lt. 0.)  alpha = 0.                                 
         if ( beta .lt. 0.)   beta = 0.                                 
c                                                                       
c        assertion:  alpha and beta are now both nonnegative.           
c                                                                       
c           make    sqrt( alpha**2 + beta**2 ) .le. 3 .                 
c                                                                       
         delta = sqrt(alpha**2 + beta**2)                               
         if (delta .le. 3.)  go to 45                                   
c        assertion:  point is outside the circle.  need to adjust.      
            tau = 3./delta                                              
            alpha = tau*alpha                                           
            beta  = tau*beta                                            
c                                                                       
c           recompute derivative values.                                
c                                                                       
   45    continue                                                       
         d(i) = alpha*slope(i)                                          
         d(i+1) = beta*slope(i)                                         
c                                                                       
   50 continue                                                          
c                                                                       
c        end of derivative assignment.                                  
c                                                                       
      return                                                            
      end  
