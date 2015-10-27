*deck %W%  %G%
c***begin prologue     getdnd
c***date written       930802   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           size, link m6200
c***author             schneider, barry (nsf)
c***source             m6200
c***purpose            data for d nodes in lebedev quadrature
c***references         lebedev paper
c
c***routines called
c***end prologue       size
      subroutine getdnd (pt,wt,r,u,w,bigd)
      implicit integer (a-z)
      real*8 pt, wt, r, u, w, bigd
      dimension pt(3,*), wt(*), bigd(1), r(1), u(1), w(1)
      common /io/ inp, iout
      call xyzw(pt(1,1),wt(1),r(1),u(1),w(1),bigd(1))
      call xyzw(pt(1,2),wt(2),r(1),-u(1),w(1),bigd(1))      
      call xyzw(pt(1,3),wt(3),-r(1),u(1),w(1),bigd(1))
      call xyzw(pt(1,4),wt(4),-r(1),-u(1),w(1),bigd(1))            
      call xyzw(pt(1,5),wt(5),r(1),u(1),-w(1),bigd(1))     
      call xyzw(pt(1,6),wt(6),r(1),-u(1),-w(1),bigd(1))     
      call xyzw(pt(1,7),wt(7),-r(1),u(1),-w(1),bigd(1))     
      call xyzw(pt(1,8),wt(8),-r(1),-u(1),-w(1),bigd(1))     
      call xyzw(pt(1,9),wt(9),r(1),w(1),u(1),bigd(1))     
      call xyzw(pt(1,10),wt(10),r(1),-w(1),u(1),bigd(1))  
      call xyzw(pt(1,11),wt(11),-r(1),w(1),u(1),bigd(1))        
      call xyzw(pt(1,12),wt(12),-r(1),-w(1),u(1),bigd(1))         
      call xyzw(pt(1,13),wt(13),r(1),w(1),-u(1),bigd(1))        
      call xyzw(pt(1,14),wt(14),r(1),-w(1),-u(1),bigd(1))  
      call xyzw(pt(1,15),wt(15),-r(1),w(1),-u(1),bigd(1))  
      call xyzw(pt(1,16),wt(16),-r(1),-w(1),-u(1),bigd(1))
      call xyzw(pt(1,17),wt(17),u(1),r(1),w(1),bigd(1))     
      call xyzw(pt(1,18),wt(18),u(1),-r(1),w(1),bigd(1))  
      call xyzw(pt(1,19),wt(19),-u(1),r(1),w(1),bigd(1))        
      call xyzw(pt(1,20),wt(20),-u(1),-r(1),w(1),bigd(1))         
      call xyzw(pt(1,21),wt(21),u(1),r(1),-w(1),bigd(1))        
      call xyzw(pt(1,22),wt(22),u(1),-r(1),-w(1),bigd(1))  
      call xyzw(pt(1,23),wt(23),-u(1),r(1),-w(1),bigd(1))  
      call xyzw(pt(1,24),wt(24),-u(1),-r(1),-w(1),bigd(1))
      call xyzw(pt(1,25),wt(25),u(1),w(1),r(1),bigd(1))     
      call xyzw(pt(1,26),wt(26),u(1),-w(1),r(1),bigd(1))  
      call xyzw(pt(1,27),wt(27),-u(1),w(1),r(1),bigd(1))        
      call xyzw(pt(1,28),wt(28),-u(1),-w(1),r(1),bigd(1))         
      call xyzw(pt(1,29),wt(29),u(1),w(1),-r(1),bigd(1))        
      call xyzw(pt(1,30),wt(30),u(1),-w(1),-r(1),bigd(1))  
      call xyzw(pt(1,31),wt(31),-u(1),w(1),-r(1),bigd(1))  
      call xyzw(pt(1,32),wt(32),-u(1),-w(1),-r(1),bigd(1))            
      call xyzw(pt(1,33),wt(33),w(1),u(1),r(1),bigd(1))     
      call xyzw(pt(1,34),wt(34),w(1),-u(1),r(1),bigd(1))  
      call xyzw(pt(1,35),wt(35),-w(1),u(1),r(1),bigd(1))        
      call xyzw(pt(1,36),wt(36),-w(1),-u(1),r(1),bigd(1))         
      call xyzw(pt(1,37),wt(37),w(1),u(1),-r(1),bigd(1))        
      call xyzw(pt(1,38),wt(38),w(1),-u(1),-r(1),bigd(1))  
      call xyzw(pt(1,39),wt(39),-w(1),u(1),-r(1),bigd(1))  
      call xyzw(pt(1,40),wt(40),-w(1),-u(1),-r(1),bigd(1))
      call xyzw(pt(1,41),wt(41),w(1),r(1),u(1),bigd(1))     
      call xyzw(pt(1,42),wt(42),w(1),-r(1),u(1),bigd(1))  
      call xyzw(pt(1,43),wt(43),-w(1),r(1),u(1),bigd(1))        
      call xyzw(pt(1,44),wt(44),-w(1),-r(1),u(1),bigd(1))         
      call xyzw(pt(1,45),wt(45),w(1),r(1),-u(1),bigd(1))        
      call xyzw(pt(1,46),wt(46),w(1),-r(1),-u(1),bigd(1))  
      call xyzw(pt(1,47),wt(47),-w(1),r(1),-u(1),bigd(1))  
      call xyzw(pt(1,48),wt(48),-w(1),-r(1),-u(1),bigd(1))               
      return
      end

