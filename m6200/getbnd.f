*deck getbnd.f
c***begin prologue     getbnd
c***date written       930802   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           size, link m6200
c***author             schneider, barry (nsf)
c***source             m6200
c***purpose            data for b nodes in lebedev quadrature
c***references         lebedev paper
c
c***routines called
c***end prologue       size
      subroutine getbnd (pt,wt,l,m,bigb,nleb)
      implicit integer (a-z)
      real*8 pt, wt, l, m
      real*8 bigb
      dimension pt(3,*), wt(*), l(4), m(4), bigb(4)
      common /io/ inp, iout
c          eight (l,l,m)
      ii=0
      call xyzw(pt(1,ii+1),wt(ii+1),l(1),l(1),m(1),bigb(1))          
      call xyzw(pt(1,ii+2),wt(ii+2),l(1),l(1),-m(1),bigb(1)) 
      call xyzw(pt(1,ii+3),wt(ii+3),l(1),-l(1),m(1),bigb(1))          
      call xyzw(pt(1,ii+4),wt(ii+4),l(1),-l(1),-m(1),bigb(1)) 
      call xyzw(pt(1,ii+5),wt(ii+5),-l(1),l(1),m(1),bigb(1)) 
      call xyzw(pt(1,ii+6),wt(ii+6),-l(1),l(1),-m(1),bigb(1))
      call xyzw(pt(1,ii+7),wt(ii+7),-l(1),-l(1),m(1),bigb(1))      
      call xyzw(pt(1,ii+8),wt(ii+8),-l(1),-l(1),-m(1),bigb(1))
c         eight (l,m,l)
      ii=8
      call xyzw(pt(1,ii+1),wt(ii+1),l(1),m(1),l(1),bigb(1))          
      call xyzw(pt(1,ii+2),wt(ii+2),l(1),m(1),-l(1),bigb(1)) 
      call xyzw(pt(1,ii+3),wt(ii+3),l(1),-m(1),l(1),bigb(1))          
      call xyzw(pt(1,ii+4),wt(ii+4),l(1),-m(1),-l(1),bigb(1)) 
      call xyzw(pt(1,ii+5),wt(ii+5),-l(1),m(1),l(1),bigb(1)) 
      call xyzw(pt(1,ii+6),wt(ii+6),-l(1),m(1),-l(1),bigb(1))
      call xyzw(pt(1,ii+7),wt(ii+7),-l(1),-m(1),l(1),bigb(1))      
      call xyzw(pt(1,ii+8),wt(ii+8),-l(1),-m(1),-l(1),bigb(1))  
c         eight (m,l,l)
      ii=16
      call xyzw(pt(1,ii+1),wt(ii+1),m(1),l(1),l(1),bigb(1))          
      call xyzw(pt(1,ii+2),wt(ii+2),m(1),l(1),-l(1),bigb(1)) 
      call xyzw(pt(1,ii+3),wt(ii+3),m(1),-l(1),l(1),bigb(1))          
      call xyzw(pt(1,ii+4),wt(ii+4),m(1),-l(1),-l(1),bigb(1)) 
      call xyzw(pt(1,ii+5),wt(ii+5),-m(1),l(1),l(1),bigb(1)) 
      call xyzw(pt(1,ii+6),wt(ii+6),-m(1),l(1),-l(1),bigb(1))
      call xyzw(pt(1,ii+7),wt(ii+7),-m(1),-l(1),l(1),bigb(1))      
      call xyzw(pt(1,ii+8),wt(ii+8),-m(1),-l(1),-l(1),bigb(1))
      if (nleb.eq.13) then
          return
      endif                         
      if(nleb.gt.11) then
c         eight (l,l,m)
          ii=24
          call xyzw(pt(1,ii+1),wt(ii+1),l(2),l(2),m(2),bigb(2))          
          call xyzw(pt(1,ii+2),wt(ii+2),l(2),l(2),-m(2),bigb(2)) 
          call xyzw(pt(1,ii+3),wt(ii+3),l(2),-l(2),m(2),bigb(2))          
          call xyzw(pt(1,ii+4),wt(ii+4),l(2),-l(2),-m(2),bigb(2)) 
          call xyzw(pt(1,ii+5),wt(ii+5),-l(2),l(2),m(2),bigb(2)) 
          call xyzw(pt(1,ii+6),wt(ii+6),-l(2),l(2),-m(2),bigb(2))
          call xyzw(pt(1,ii+7),wt(ii+7),-l(2),-l(2),m(2),bigb(2))      
          call xyzw(pt(1,ii+8),wt(ii+8),-l(2),-l(2),-m(2),bigb(2))          
c         eight (l,m,l)
          ii=32
          call xyzw(pt(1,ii+1),wt(ii+1),l(2),m(2),l(2),bigb(2))          
          call xyzw(pt(1,ii+2),wt(ii+2),l(2),m(2),-l(2),bigb(2)) 
          call xyzw(pt(1,ii+3),wt(ii+3),l(2),-m(2),l(2),bigb(2))          
          call xyzw(pt(1,ii+4),wt(ii+4),l(2),-m(2),-l(2),bigb(2)) 
          call xyzw(pt(1,ii+5),wt(ii+5),-l(2),m(2),l(2),bigb(2)) 
          call xyzw(pt(1,ii+6),wt(ii+6),-l(2),m(2),-l(2),bigb(2))
          call xyzw(pt(1,ii+7),wt(ii+7),-l(2),-m(2),l(2),bigb(2))      
          call xyzw(pt(1,ii+8),wt(ii+8),-l(2),-m(2),-l(2),bigb(2))    
c         eight (m,l,l)
          ii=40
          call xyzw(pt(1,ii+1),wt(ii+1),m(2),l(2),l(2),bigb(2))          
          call xyzw(pt(1,ii+2),wt(ii+2),m(2),l(2),-l(2),bigb(2)) 
          call xyzw(pt(1,ii+3),wt(ii+3),m(2),-l(2),l(2),bigb(2))          
          call xyzw(pt(1,ii+4),wt(ii+4),m(2),-l(2),-l(2),bigb(2)) 
          call xyzw(pt(1,ii+5),wt(ii+5),-m(2),l(2),l(2),bigb(2)) 
          call xyzw(pt(1,ii+6),wt(ii+6),-m(2),l(2),-l(2),bigb(2))
          call xyzw(pt(1,ii+7),wt(ii+7),-m(2),-l(2),l(2),bigb(2))      
          call xyzw(pt(1,ii+8),wt(ii+8),-m(2),-l(2),-l(2),bigb(2))              
c         eight (l,l,m)
          ii=48
          call xyzw(pt(1,ii+1),wt(ii+1),l(3),l(3),m(3),bigb(3))          
          call xyzw(pt(1,ii+2),wt(ii+2),l(3),l(3),-m(3),bigb(3)) 
          call xyzw(pt(1,ii+3),wt(ii+3),l(3),-l(3),m(3),bigb(3))          
          call xyzw(pt(1,ii+4),wt(ii+4),l(3),-l(3),-m(3),bigb(3)) 
          call xyzw(pt(1,ii+5),wt(ii+5),-l(3),l(3),m(3),bigb(3)) 
          call xyzw(pt(1,ii+6),wt(ii+6),-l(3),l(3),-m(3),bigb(3))
          call xyzw(pt(1,ii+7),wt(ii+7),-l(3),-l(3),m(3),bigb(3))      
          call xyzw(pt(1,ii+8),wt(ii+8),-l(3),-l(3),-m(3),bigb(3))          
c         eight (l,m,l)
          ii=56
          call xyzw(pt(1,ii+1),wt(ii+1),l(3),m(3),l(3),bigb(3))          
          call xyzw(pt(1,ii+2),wt(ii+2),l(3),m(3),-l(3),bigb(3)) 
          call xyzw(pt(1,ii+3),wt(ii+3),l(3),-m(3),l(3),bigb(3))          
          call xyzw(pt(1,ii+4),wt(ii+4),l(3),-m(3),-l(3),bigb(3)) 
          call xyzw(pt(1,ii+5),wt(ii+5),-l(3),m(3),l(3),bigb(3)) 
          call xyzw(pt(1,ii+6),wt(ii+6),-l(3),m(3),-l(3),bigb(3))
          call xyzw(pt(1,ii+7),wt(ii+7),-l(3),-m(3),l(3),bigb(3))      
          call xyzw(pt(1,ii+8),wt(ii+8),-l(3),-m(3),-l(3),bigb(3))    
c         eight (m,l,l)
          ii=64
          call xyzw(pt(1,ii+1),wt(ii+1),m(3),l(3),l(3),bigb(3))          
          call xyzw(pt(1,ii+2),wt(ii+2),m(3),l(3),-l(3),bigb(3)) 
          call xyzw(pt(1,ii+3),wt(ii+3),m(3),-l(3),l(3),bigb(3))          
          call xyzw(pt(1,ii+4),wt(ii+4),m(3),-l(3),-l(3),bigb(3)) 
          call xyzw(pt(1,ii+5),wt(ii+5),-m(3),l(3),l(3),bigb(3)) 
          call xyzw(pt(1,ii+6),wt(ii+6),-m(3),l(3),-l(3),bigb(3))
          call xyzw(pt(1,ii+7),wt(ii+7),-m(3),-l(3),l(3),bigb(3))      
          call xyzw(pt(1,ii+8),wt(ii+8),-m(3),-l(3),-l(3),bigb(3))                
          if (nleb.gt.17) then
c         eight (l,l,m)
              ii=72
              call xyzw(pt(1,ii+1),wt(ii+1),l(4),l(4),m(4),bigb(4))          
              call xyzw(pt(1,ii+2),wt(ii+2),l(4),l(4),-m(4),bigb(4)) 
              call xyzw(pt(1,ii+3),wt(ii+3),l(4),-l(4),m(4),bigb(4))          
              call xyzw(pt(1,ii+4),wt(ii+4),l(4),-l(4),-m(4),bigb(4)) 
              call xyzw(pt(1,ii+5),wt(ii+5),-l(4),l(4),m(4),bigb(4)) 
              call xyzw(pt(1,ii+6),wt(ii+6),-l(4),l(4),-m(4),bigb(4))
              call xyzw(pt(1,ii+7),wt(ii+7),-l(4),-l(4),m(4),bigb(4))      
              call xyzw(pt(1,ii+8),wt(ii+8),-l(4),-l(4),-m(4),bigb(4))          
c         eight (l,m,l)
              ii=80
              call xyzw(pt(1,ii+1),wt(ii+1),l(4),m(4),l(4),bigb(4))          
              call xyzw(pt(1,ii+2),wt(ii+2),l(4),m(4),-l(4),bigb(4)) 
              call xyzw(pt(1,ii+3),wt(ii+3),l(4),-m(4),l(4),bigb(4))          
              call xyzw(pt(1,ii+4),wt(ii+4),l(4),-m(4),-l(4),bigb(4)) 
              call xyzw(pt(1,ii+5),wt(ii+5),-l(4),m(4),l(4),bigb(4)) 
              call xyzw(pt(1,ii+6),wt(ii+6),-l(4),m(4),-l(4),bigb(4))
              call xyzw(pt(1,ii+7),wt(ii+7),-l(4),-m(4),l(4),bigb(4))      
              call xyzw(pt(1,ii+8),wt(ii+8),-l(4),-m(4),-l(4),bigb(4))   
c         eight (m,l,l)
              ii=88
              call xyzw(pt(1,ii+1),wt(ii+1),m(4),l(4),l(4),bigb(4))          
              call xyzw(pt(1,ii+2),wt(ii+2),m(4),l(4),-l(4),bigb(4)) 
              call xyzw(pt(1,ii+3),wt(ii+3),m(4),-l(4),l(4),bigb(4))          
              call xyzw(pt(1,ii+4),wt(ii+4),m(4),-l(4),-l(4),bigb(4)) 
              call xyzw(pt(1,ii+5),wt(ii+5),-m(4),l(4),l(4),bigb(4)) 
              call xyzw(pt(1,ii+6),wt(ii+6),-m(4),l(4),-l(4),bigb(4))
              call xyzw(pt(1,ii+7),wt(ii+7),-m(4),-l(4),l(4),bigb(4))      
              call xyzw(pt(1,ii+8),wt(ii+8),-m(4),-l(4),-l(4),bigb(4))                 
          endif    
      endif                                           
      return
      end

