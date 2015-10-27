************************************************************************
      SUBROUTINE  CUBFEV (T0, C0, C1, C2, C3, K, T, IDERIV, CVAL)       
      implicit real*8(a-h,o-z)
      integer  k, ideriv                                                
      real*8  t0, c0, c1, c2, c3, t(k), cval(ideriv,k)                    
c                                                                       
c-----------------------------------------------------------------------
c                                                                       
c     evaluate the cubic function                                       
c                                                                       
c        c(t) = c0 + c1*(t-t0) + c2*(t-t0)**2 + c3*(t-t0)**3            
c                                                                       
c     and its first (ideriv-1) derivatives at t(i), i=1(1)k.            
c                                                                       
c-----------------------------------------------------------------------
c                                                                       
c     on input:                                                         
c        t0      is the expansion point (usually the left endpoint of   
c                an interval).                                          
c        c0,...,c3  are the cubic coefficients.                         
c        k       is the number of evaluation points.                    
c                it is assumed that  k.ge.1 .  (not checked)            
c        t       is the array of evaluation points.                     
c        ideriv  is the number of functions to be evaluated.            
c                it is assumed that ideriv = 1, 2, or 3.  (not checked)
c                                                                       
c     on output:                                                        
c        cval    will contain the desired function values, as follows. 
c                                                                      
c                 cval(1,i) = c(t(i))                       )         
c                 cval(2,i) = (d/dt) c(t(i))  (ideriv.ge.2) )  i=1(1)k.
c                 cval(3,i) = (d/dt)**2 c(t(i))  (ideriv=3) )          
c                                                                      
c--------------------------------------------------------------------
c                                                                    
c     algorithm by:  f. n. fritsch, lawrence livermore laboratory,  a
c                    r. p. dickinson, jr., chabot college.           
c     programmed by:  f. n. fritsch.                                 
c     date last changed:  14 march 1979  (fnf)                      
c                                                                    
c     change record:                                                   
c          initial implementation.                          
c                                                                    
c---------------------------------------------------------------------
c                                                                    
c        local declarations.                                          
c                                                                    
      integer  i                                                      
      real*8  d2, d3, delta                                            
c                                                                    
c        evaluation loop.                                           
c                                                                    
      do 500  i = 1, k                                             
         delta = t(i) - t0                                          
         cval(1,i) = c0 + delta*(c1 + delta*(c2 + delta*c3))       
c                                                                    
c        note..  if derivatives will never be needed, may             
c                delete from here to (but not including) statement 50
c                may also delete argument ideriv and the first index  
c                of cval in this case.                                
         if (ideriv .eq. 1)  go to 500                             
c                                                                     
c           evaluate first derivatives.                               
c                                                                     
         d2 = 2.*c2                                                
         d3 = 3.*c3                                                 
         cval(2,i) = c1 + delta*(d2 + delta*d3)                     
c                                                                    
c        note..  if second derivatives will never be needed, may        
c                delete from here to (but not including) statement 500 .
         if (ideriv .eq. 2)  go to 500                                  
c                                                                       
c           compute second derivatives.                                 
c                                                                       
         cval(3,i) = d2 + 2.*delta*d3                                   
  500 continue                                                          
c                                                                       
      return                                                            
c                                                                       
      end                                                                
