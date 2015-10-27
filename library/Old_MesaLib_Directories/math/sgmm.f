*deck  @(#)sgmm.f	5.1 11/6/94 
      subroutine sgmm(m,n,l,a,ia,b,ib,c,ic,jtrpos,job)                 
c***begin prologue          sgmm                                               
c***date written                     (yymmdd)
c***revision date           910802   (yymmdd)
c   2 august  1991  rlm at lanl
c     changing calls to sgmv to the new blas sgemv
c***category no.            d1b6 
c***keywords                matrix multiplication 
c***author                  jordan, tom, los alamos national laboratory 
c***purpose                 performs matrix multiplication with all forms
c                           of transposition of the input and output matrices
c                           and with four forms of accumulation.
c***description                                                         
c                                                                       
c     given storage arrays a,b,c, this routine may                      
c                                  t       t                            
c     be used to compute ab=c, ab=c ,----,a b =c                        
c     (see the jtrpos parameter).  in addition, for each form           
c     of transposition, the following four forms of accumulation        
c     are possible:  ab=c, c+ab=c, -ab=c and c-ab=c (see                
c     the job parameter).  much explicit transposition can be           
c     avoided by using this routine and the transposition is            
c     achieved at no cost.                                              
c     this is a general routine for multiplying two arbitrary matrices  
c     a and b with four forms of accumulation (see the job parameter)   
c     and with any of the combination of transposes (see the jtrpos     
c     parameter) including the result matrix c.                         
c                                                                       
c     -arguments-                                                       
c                                                                       
c      on input                                                         
c       m = the number of rows (columns) of a (a(transpose)) and        
c           the number of rows (columns) of c (c(transpose)).           
c       n = the number of colunms (rows) of a (a(transpose)) and the    
c           the number of rows (columns) of b (b(transpose)).           
c       l = the number of columns (rows) of b (b(transpose)) and the    
c           number of columns (rows) of c (c(transpose)).               
c       a, b, c = a 2-d array for a (a(transpose), etc.                 
c       ia, ib, ic = the first dimension of a, b, c in the calling      
c                    program.                                           
c       jtrpos = a transposition indicator defined as follows:          
c                                                                       
c                  value                transposition                   
c                    0                  a   *b   =c                     
c                    1                  a   *b   =c(t)                  
c                    2                  a   *b(t)=c                     
c                    3                  a   *b(t)=c(t)                  
c                    4                  a(t)*b   =c                     
c                    5                  a(t)*b   =c(t)                  
c                    6                  a(t)*b(t)=c                     
c                    7                  a(t)*b(t)=c(t)                  
c                                                                       
c         job = a task designator                                       
c              -1   c =  -ab                                            
c              -2   c = c-ab                                            
c              +1   c =  +ab                                            
c              +2   c = c+ab                                            
c                                                                       
c      on output                                                        
c      c = a 2-d array containing the specified product                 
c           c (c(transpose)).                                           
c                                                                       
c***references              (none)
c***routines called         sgemv                                              
c***special implementation  cftmath                                     
c***end prologue            sgmm                                                 
      implicit real*8 (a-h,o-z)
      dimension a(ia,n),b(ib,l),c(ic,l)                                 
      real*8 zero,one
c
      parameter (zero=0.0d+0,one=1.0d+0)
c
c***first executable statement  sgmm                                   
      if (l .le. 0)                              return                 
      jb = isign(iabs(job)+2*(jtrpos/4),job)                            
c
c     set up appropriate multipliers for sgemv
      facab=one
      if(jb.lt.0) facab=-one
      facc=zero
      if(mod(jb,2).eq.0) facc=one
c
      jump = jtrpos+1                                                   
      go to (10,20,30,40,50,60,70,80),jump                              
   10 continue                                                          
c                                                                       
c** this routine forms the matrix product c=ab, the routine calls a cal 
c   kernel to multiply a matrix times a vector -                        
cinput                                                                  
c   m = the number of rows of a and c.                                  
c   n = the number of columns of a and rows of b                        
c   l = the number of columns of b and c.                               
c   a,b,c = 2-d arrays containing the matrices a,b and c                
c   ia,ib,ic= the first dimensions of arrays a,b and c                  
coutput                                                                 
c   c = the product matrix c=ab                                         
c                                                                       
      do  15    j=1,l                                                   
c        call sgmv (m,n,a,ia,b(1,j),1,c(1,j),1,jb)                    
         call sgemv('n',m,n,facab,a,ia,b(1,j),1,facc,c(1,j),1)
   15 continue
      return                                                            
   50 continue                                                          
c                                                                       
c** this routine forms the matrix product a(transpose)b=c               
cinput                                                                  
c  m = the number of columns of a and rows of c.                        
c  n = the number of rows of a and b.                                   
c  l = the number of columns of b and c                                 
c  a,b,c = 2-d arrays containing matrices a,b and c.                    
c  ia,ib,ic = the first dimensions of arrays a,b,c.                     
coutput                                                                 
c  c = the product matrix c=a(transpose)b                               
c                                                                       
      do 55    j=1,l                                                    
c        call sgmv (m,n,a,ia,b(1,j),1,c(1,j),1,jb)                    
         call sgemv ('t',n,m,facab,a,ia,b(1,j),1,facc,c(1,j),1)
   55 continue
      return                                                            
   30 continue                                                          
c                                                                       
ccc this routine forms the matrix product c=ab(transpose)               
c                                                                       
cinput                                                                  
c  m = the number of rows of a and c.                                   
c  n = the number of columns of a and c.                                
c  l = the number of rows of b and columns of c.                        
c  a,b,c= 2-d arrays containing the matrices a,b and c.                 
c  ia,ib,ic= the first dimension of arrays a,b,c.                       
coutput                                                                 
c  c = the product matrix c = ab(transpose)                             
c                                                                       
      do 35    j=1,l                                                    
c        call sgmv (m,n,a,ia,b(j,1),ib,c(1,j),1,jb)                   
         call sgemv('n',m,n,facab,a,ia,b(j,1),ib,facc,c(1,j),1)
   35 continue
      return                                                            
   20 continue                                                          
c                                                                       
ccc this routine forms the matrix product ab and stores its transpose   
c   in c.                                                               
c                                                                       
cinput                                                                  
c  m = the number of rows of a and columns of c.                        
c  n = the number of columns of a and rows of b.                        
c  l = the number of columns of b and rows of c.                        
c  a,b,c = 2-d arrays containing matrices a,b and c.                    
c  ia,ib,ic = the first dimension of arrays a,b,c.                      
coutput                                                                 
c  c  contains the transpose of ab.                                     
c                                                                       
      do 25    j=1,l                                                    
c        call sgmv (m,n,a,ia,b(1,j),1,c(j,1),ic,jb)                   
         call sgemv ('n',m,n,facab,a,ia,b(1,j),1,facc,c(j,1),ic)
   25 continue
      return                                                            
   70 continue                                                          
c                                                                       
ccc  this routine forms the matrix product c=a(transpose)b(transpose)   
c                                                                       
cinput                                                                  
c  m = the number of columns of a and rows of c.                        
c  n = the number of rows of a and columns of b.                        
c  l = the number of rows of b and columns of c                         
c  a,b,c = 2-d arrays containing matrices a,b and c.                    
c  ia,ib,ic= the first dimensions of array a,b,c.                       
coutput                                                                 
c  c = the product matrix c=a(transpose)b(transpose)                    
c                                                                       
      do 75    j=1,l                                                    
c        call sgmv (m,n,a,ia,b(j,1),ib,c(1,j),1,jb)                   
         call sgemv ('t',n,m,facab,a,ia,b(j,1),ib,facc,c(1,j),1)
   75 continue
      return                                                            
   60 continue                                                          
c                                                                       
c** this routine forms the matrix product a(transpose)b and stores      
c   its transpose in c.                                                 
cinput                                                                  
c  m = the number of columns of a and c.                                
c  n = the number of rows of a and b.                                   
c  l = the number of columns of b and rows of c.                        
c  a,b,c = 2-d arrays containing matrices a,b and c.                    
c  ia,ib,ic = the first dimensions of arrays a,b,c.                     
coutput                                                                 
c  c contains the transpose of the product matrix a(transpose)b         
c                                                                       
      do 65    j=1,l                                                    
c        call sgmv (m,n,a,ia,b(1,j),1,c(j,1),ic,jb)                   
         call sgemv ('t',n,m,facab,a,ia,b(1,j),1,facc,c(j,1),ic) 
   65 continue
      return                                                            
   40 continue                                                          
c                                                                       
c** this routine forms the matrix product ab(transpose) and stores      
c   its transpose in c                                              
cinput                                                                  
c  m = the number of rows of a and columns of c.                        
c  n = the number of columns of a and b.                                
c  l = the number of rows of b and c                                    
c  a,b,c = 2-d arrays containing matrices a,b and c.                    
c  ia,ib,ic = the first dimensions of arrays a,b,c.                     
coutput                                                                 
c  c contains the transpose of the product ab(transpose)                
c                                                                       
      do 45    j=1,l                                                    
c        call sgmv (m,n,a,ia,b(j,1),ib,c(j,1),ic,jb)                  
         call sgemv ('n',m,n,facab,a,ia,b(j,1),ib,facc,c(j,1),ic)
   45 continue
      return                                                            
   80 continue                                                          
c                                                                       
c** this routine forms the matrix product a(transpose)b(transpose) and  
c   stores its transpose in c.                                          
cinput                                                                  
c  m = the number of columns of a and c.                                
c  n = the number of rows of a and columns of b.                        
c  l = the number of rows of b and c.                                   
c  a,b,c = 2-d arrays containing matrices a,b and c.                    
c  ia,ib,ic = the first dimensions of arrays a,b,c.                     
coutput                                                                 
c  c contains the transpose of the product a(transpose)b(transpose)     
      do 85    j=1,l                                                    
c        call sgmv (m,n,a,ia,b(j,1),ib,c(j,1),ic,jb)                  
         call sgemv ('t',n,m,facab,a,ia,b(j,1),ib,facc,c(j,1),ic)
   85 continue
      return                                                            
      end                                                               
