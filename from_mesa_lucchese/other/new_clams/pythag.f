*deck %W%  %G%
      function pythag(a,b)                                       
c***begin prologue     pythag
c***date written       yymmdd  
c***revision date      ymmdd
c
c***keywords           
c***author             unknown 
c***source             %W%   %G%
c***purpose            
c                                                                       
c     finds sqrt(a**2+b**2) without overflow or destructive underflow   
c***description
c     
c    
c
c***references
c
c***routines called (none)
c
c***end prologue       pythag
      implicit real*8 (a-h,o-z)
      real*8 a,b,pythag
c                                                                       
      real*8 p,q,r,s,t                                                  
c***first executable statement  pythag                                  
      p = max(abs(a),abs(b))
      q = min(abs(a),abs(b))
      if (q .eq. 0.0d0) go to 20                                        
   10 continue                                                          
         r = (q/p)**2                                                   
         t = 4.0d0 + r                                                  
         if (t .eq. 4.0d0) go to 20                                     
         s = r/t                                                        
         p = p + 2.0d0*p*s                                              
         q = q*s                                                        
      go to 10                                                          
   20 pythag = p                                                        
      return                                                            
      end                                                               
