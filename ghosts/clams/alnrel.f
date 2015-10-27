*deck %W%  %G%
      function alnrel(x)                                                
c***begin prologue  alnrel                                              
c***date written   770401   (yymmdd)                                    
c***revision date  820801   (yymmdd)                                    
c***category no.  c4b                                                   
c***keywords  elementary function,logarithm,relative                    
c***author  fullerton, w., (lanl)                                       
c***purpose  evaluates ln(1+x) accurate in the sense of relative error. 
c***description                                                         
c                                                                       
c alnrel(x) evaluates ln(1+x) accurately in the sense of relative       
c error when x is very small.  this routine must be used to             
c maintain relative error accuracy whenever x is small and              
c accurately known.                                                     
c                                                                       
c series for alnr       on the interval -3.75000d-01 to  3.75000d-01    
c                                        with weighted error   1.93d-17 
c                                         log weighted error  16.72     
c                               significant figures required  16.44     
c                                    decimal places required  17.40     
c***references  (none)                                                  
c***routines called  csevl,inits,r1mach,xerror                          
c***end prologue  alnrel                                                
c
      implicit real*8 (a-h,o-z)
c
      dimension alnrcs(23)                                              
      data alnrcs( 1) /   1.0378693562 743770d0 /                       
      data alnrcs( 2) /   -.1336430150 4908918d0 /                      
      data alnrcs( 3) /    .0194082491 35520563d0 /                     
      data alnrcs( 4) /   -.0030107551 12753577d0 /                     
      data alnrcs( 5) /    .0004869461 47971548d0 /                     
      data alnrcs( 6) /   -.0000810548 81893175d0 /                     
      data alnrcs( 7) /    .0000137788 47799559d0 /                     
      data alnrcs( 8) /   -.0000023802 21089435d0 /                     
      data alnrcs( 9) /    .0000004164 04162138d0 /                     
      data alnrcs(10) /   -.0000000735 95828378d0 /                     
      data alnrcs(11) /    .0000000131 17611876d0 /                     
      data alnrcs(12) /   -.0000000023 54670931d0 /                     
      data alnrcs(13) /    .0000000004 25227732d0 /                     
      data alnrcs(14) /   -.0000000000 77190894d0 /                     
      data alnrcs(15) /    .0000000000 14075746d0 /                     
      data alnrcs(16) /   -.0000000000 02576907d0 /                     
      data alnrcs(17) /    .0000000000 00473424d0 /                     
      data alnrcs(18) /   -.0000000000 00087249d0 /                     
      data alnrcs(19) /    .0000000000 00016124d0 /                     
      data alnrcs(20) /   -.0000000000 00002987d0 /                     
      data alnrcs(21) /    .0000000000 00000554d0 /                     
      data alnrcs(22) /   -.0000000000 00000103d0 /                     
      data alnrcs(23) /    .0000000000 00000019d0 /                     
      data nlnrel, xmin /0, 0./                                         
c***first executable statement  alnrel                                  
      if (nlnrel.ne.0) go to 10                                         
      nlnrel = inits (alnrcs, 23, 0.1*r1mach(3))                        
      xmin = -1.0 + sqrt(r1mach(4))                                     
c                                                                       
 10   if (x.le.(-1.0)) call xerror ( 'alnrel  x is le -1', 18, 2, 2)    
      if (x.lt.xmin) call xerror ( 'alnrel  answer lt half precision bec
     1ause x too near -1', 54,    1, 1)                                 
c                                                                       
      if (abs(x).le.0.375) alnrel = x*(1. -                             
     1  x*csevl (x/.375, alnrcs, nlnrel))                               
      if (abs(x).gt.0.375) alnrel = log (1.0+x)                        
c                                                                       
      return                                                            
      end                                                               
