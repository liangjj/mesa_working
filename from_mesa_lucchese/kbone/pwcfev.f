************************************************************************
      SUBROUTINE  PWCFEV (IDERIV, N, T, F, D, NE, TE, FE)               
      implicit real*8(a-h,o-z)
      integer  ideriv, n, ne                                          
      real*8  t(n), f(n), d(n), te(ne), fe(ideriv,ne)                     
c                                                                   
c----------------------------------------------------------------------
c                                                                       
c        evaluate the first (ideriv-1) derivatives of the piecewise     
c        cubic function defined by  n, t, f, d  at the points  te(i),   
c        i = 1(1)ne.                                                   
c                                                                       
c     this is version (d).  it:                                         
c     (1) converts on the fly to polynomial form.                       
c     (2) collects all points in segment il before calling cubfev.      
c                                                                      
c-----------------------------------------------------------------------
c                                                                      
c     on input:                                                        
c        ideriv  indicates how many derivatives are desired.          
c                restriction:  1 .le. ideriv .le. 3 (not checked).     
c                note:  ideriv=1 implies only function values requested
c        n       is the number of data points.                        
c                restriction:  n .ge. 2  (not checked).                
c        t       is the array of independent variable values.           
c                the search procedure assumes that t is strictly        
c                increasing.  (not checked)                          
c        f       is the corresponding array of function values.         
c        d       is the corresponding array of derivative values.      
c        ne      is the number of points at which evaluation is desired.
c        te      is the array of evaluation points.                     
c                the search procedure assumes that te is monotone       
c                increasing.  (not checked)                            
c                                                                       
c     on return.                                                        
c        fe      contains the function values, as follows:         
c                  fe(id,i) is the value of the (id-1)-st derivative of 
c                           the piecewise cubic at te(i), i=1(1)ne,     
c                           id=1(1)ideriv.                            
c               ---------                                               
c                caution:  the order of the indices of fe has been      
c               ---------  changed since the 78-12-22 version of pwcfev.
c                                                                      
c     other routines used:  cubfev.                                     
c                                                                 
c----------------------------------------------------------------------
c                                                                      
c     algorithm by:  f. n. fritsch, lawrence livermore laboratory.      
c     programmed by: f. n. fritsch.                                  
c     date last changed:  15 march 1979 (fnf)                           
c                                                                       
c     change record:                                                   
c        78-12-07   minor cosmetic changes to get ready for library.   
c        78-12-21   minor modifications to search routine.              
c        78-12-22   changed to unnormalized hermite basis functions.  
c                   also changed argument name ye to fe.              
c        a  changed order of indices of fe, to allow variable   
c                   number of evaluation points in calling program      
c                   without requiring another argument.                 
c                   also changed argument names x to t and xe to te.    
c        c  modified to convert on the fly to polynomial form.  
c        d  changed search procedure to only compute polynomial 
c                   coefficients once per interval.                     
c        79-03-15   changed index il to ir-1.                         
c                                                                      
c-----------------------------------------------------------------------
c                                                                      
c        local declarations.                                            
c                                                                     
      integer  ir, j, jfirst, ki                                        
      real*8  c2, c3, delta, r, rh, s                                     
c                                                                     
c        loop over intervals.  (  interval index is  il = ir-1  .)     
c                              (interval is t(il).le.t.lt.t(ir) .)    
      jfirst = 1                                                        
      do 50  ir = 2, n                                                  
c                                                                       
c           skip out of loop if have processed all evaluation points.   
c                                                                       
         if (jfirst .gt. ne)  go to 51                                  
c                                                                       
c           locate all points in interval.                              
c                                                                       
         do 20  j = jfirst, ne                                          
            if (te(j) .ge. t(ir))  go to 30                             
   20    continue                                                       
         j = ne + 1                                                     
         go to 40                                                       
c                                                                       
c           have located first point beyond interval.                   
c                                                                       
   30    continue                                                       
         if (ir .eq. n)  j = ne + 1                                     
c                                                                       
   40    continue                                                       
         ki = j - jfirst                                                
c                                                                       
c           skip evaluation if no points in interval.                   
c                                                                       
         if (ki .eq. 0)  go to 50                                       
c                                                                       
c           compute polynomial coefficients.                            
c                                                                       
         rh = 1./(t(ir) - t(ir-1))                                    
         delta = (f(ir) - f(ir-1))*rh                                  
         s = d(ir-1) - delta                                           
         r = s + (d(ir) - delta)                                     
         c2 = -(s + r)*rh                                               
         c3 = r*rh*rh                                                   
c                                                                       
c           evaluate cubic and appropriate derivatives at te(i),      
c           i = jfirst (1) j-1 .                                       
c                                                                       
         call cubfev (t(ir-1), f(ir-1), d(ir-1), c2, c3,              
     *                ki, te(jfirst), ideriv, fe(1,jfirst))             
c                                                                     
         jfirst = j                                                   
   50 continue                                                          
c                                                                    
   51 continue                                                          
      return                                                          
      end    

