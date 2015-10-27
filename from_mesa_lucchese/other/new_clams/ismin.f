*deck %W% %G%
      function ismin(n,sx,incx)
c***begin prologue  ismin     
c***date written   790614   (yymmdd)
c***revision date  860401   (yymmdd)
c***category no.  d1a2             
c***keywords  vector,minimum,index
c***author  kahaner, d. k., los alamos national laboratory  
c***purpose  find the smallest index of a minimum element of a vector. 
c***description                                                       
c                                                                    
c   this function finds the smallest index of a minimum element of a
c   real array sx whose n elements are stored sequentially with      
c   spacing incx >= 1.  if n <= 0, the value zero is returned.      
c   thus, if i = ismin(n,sx,1), then sx(i) is an element of array sx 
c   of minimum value.                                               
c                                                                  
c   description of parameters                                     
c                                                                
c    --input--                                                  
c        n  number of elements in input vector                 
c       sx  single precision vector with n elements           
c     incx  storage spacing between elements of sx   
c                                                   
c    --output--                                    
c    ismin  smallest index (zero if n .le. 0)     
c                                                
c***references  (none)                          
c***routines called  (none)                    
c***end prologue  ismin                       
      integer ismin
      real*8 sx(*),smin  
      integer i,incx,ix,n 
c***first executable statement  ismin   
      ismin = 0                        
      if( n .lt. 1 ) return           
      ismin = 1                      
      if(n.eq.1)return              
      if(incx.eq.1)go to 20        
c                                 
c        code for increment not equal to 1
c                                        
      ix = 1                            
      smin = (sx(1))                   
      ix = ix + incx                  
      do 10 i = 2,n                  
         if((sx(ix)).ge.smin) go to 5                     
         ismin = i                                       
         smin = (sx(ix))                                
    5    ix = ix + incx                                
   10 continue                                        
      return                                         
c                                                   
c        code for increment equal to 1             
c                                                 
   20 smin = (sx(1))                             
      do 30 i = 2,n                             
         if((sx(i)).ge.smin) go to 30          
         ismin = i                            
         smin = (sx(i))                      
   30 continue                              
      return                               
      end                                 
