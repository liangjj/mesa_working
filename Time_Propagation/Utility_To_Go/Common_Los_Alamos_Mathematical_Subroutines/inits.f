      function inits(os,nos,eta)                                        
c***begin prologue  inits                                               
c***date written   770401   (yymmdd)                                    
c***revision date  820801   (yymmdd)                                    
c***category no.  c3a2                                                  
c***keywords  initialize,orthogonal series,special function             
c***author  fullerton, w., (lanl)                                       
c***purpose  initializes an orthogonal series so that it defines the    
c            number of terms to carry in the series to meet a specified 
c            error.                                                     
c***description                                                         
c                                                                       
c initialize the orthogonal series so that inits is the number of terms 
c needed to insure the error is no larger than eta.  ordinarily, eta    
c will be chosen to be one-tenth machine precision.                     
c                                                                       
c             input arguments --                                        
c os     array of nos coefficients in an orthogonal series.             
c nos    number of coefficients in os.                                  
c eta    requested accuracy of series.                                  
c***references  (none)                                                  
c***routines called  xerror                                             
c***end prologue  inits                                                 
c
      implicit real*8 (a-h,o-z)
c
      dimension os(nos)                                                 
c***first executable statement  inits                                   
      if (nos.lt.1) call lnkerr ( 'inits   number of coefficients lt 1')
c                                                                       
      err = 0.d0
      do 10 ii=1,nos                                                    
        i = nos + 1 - ii                                                
        err = err + abs(os(i))                                          
        if (err.gt.eta) go to 20                                        
 10   continue                                                          
c                                                                       
 20   if (i.eq.nos) call lnkerr ( 'inits   eta may be too small')
      inits = i                                                         
c                                                                       
      return                                                            
      end                                                               
