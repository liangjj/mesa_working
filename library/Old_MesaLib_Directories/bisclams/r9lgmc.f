      function r9lgmc(x)                                                
c***begin prologue  r9lgmc                                              
c***date written   770801   (yymmdd)                                    
c***revision date  820801   (yymmdd)                                    
c***category no.  c7e                                                   
c***keywords  correction factor,log gamma,special function              
c***author  fullerton, w., (lanl)                                       
c***purpose  computes the log gamma correction factor so that           
c            log(gamma(x)) = log(sqrt(2*pi)) + (x-.5)*log(x) - x     
c            + r9lgmc(x)                                                
c***description                                                         
c                                                                       
c compute the log gamma correction factor for x .ge. 10.0 so that       
c  log (gamma(x)) = log(sqrt(2*pi)) + (x-.5)*log(x) - x + r9lgmc(x)  
c                                                                       
c series for algm       on the interval  0.          to  1.00000d-02    
c                                        with weighted error   3.40d-16 
c                                         log weighted error  15.47     
c                               significant figures required  14.39     
c                                    decimal places required  15.86     
c***references  (none)                                                  
c***routines called  csevl,inits,r1mach,lnkerr                          
c***end prologue  r9lgmc                                                
c
      implicit real*8 (a-h,o-z)
c
      dimension algmcs(6)                                               
      data algmcs( 1) /    .1666389480 45186d0 /                        
      data algmcs( 2) /   -.0000138494 817606d0 /                       
      data algmcs( 3) /    .0000000098 108256d0 /                       
      data algmcs( 4) /   -.0000000000 180912d0 /                       
      data algmcs( 5) /    .0000000000 000622d0 /                       
      data algmcs( 6) /   -.0000000000 000003d0 /                       
      data nalgm, xbig, xmax / 0, 2*0.0d0 /                               
c***first executable statement  r9lgmc                                  
      if (nalgm.ne.0) go to 10                                          
      nalgm = inits (algmcs, 6, r1mach(3))                              
      xbig = 1.0d0/sqrt(r1mach(3))                                        
      xmax = exp (min(log(r1mach(2)/12.0d0), -log(12.0d0*r1mach(1))) )  
c                                                                       
 10   if (x.lt.10.0d0) call lnkerr ( 'r9lgmc  x must be ge 10')
      if (x.ge.xmax) go to 20                                           
c                                                                       
      r9lgmc = 1.0d0/(12.0d0*x)                                             
      if (x.lt.xbig) r9lgmc = csevl (2.0d0*(10.d0/x)**2-1.d0,
     1                               algmcs, nalgm)/x
      return                                                            
c                                                                       
 20   r9lgmc = 0.0d0
      call lnkerr ( 'r9lgmc  x so big r9lgmc underflows')
      return                                                            
c                                                                       
      end                                                               
