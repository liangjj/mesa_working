*deck %W%  %G%
c***begin prologue     getcnd
c***date written       930802   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           size, link m6200
c***author             schneider, barry (nsf)
c***source             m6200
c***purpose            data for c nodes in lebedev quadrature
c***references         lebedev paper
c
c***routines called
c***end prologue       size
      subroutine getcnd (pt,wt,p,q,bigc)
      implicit integer (a-z)
      real*8 pt, wt, p, q, bigc
      dimension pt(3,*), wt(*), bigc(1), p(1), q(1)
      common /io/ inp, iout
      call xyzw(pt(1,1),wt(1),p(1),q(1),0.d0,bigc(1))
      call xyzw(pt(1,2),wt(2),p(1),-q(1),0.d0,bigc(1))
      call xyzw(pt(1,3),wt(3),-p(1),q(1),0.d0,bigc(1))
      call xyzw(pt(1,4),wt(4),-p(1),-q(1),0.d0,bigc(1))
      call xyzw(pt(1,5),wt(5),p(1),0.d0,q(1),bigc(1))
      call xyzw(pt(1,6),wt(6),p(1),0.d0,-q(1),bigc(1))
      call xyzw(pt(1,7),wt(7),-p(1),0.d0,q(1),bigc(1))
      call xyzw(pt(1,8),wt(8),-p(1),0.d0,-q(1),bigc(1))            
      call xyzw(pt(1,9),wt(9),0.d0,p(1),q(1))
      call xyzw(pt(1,10),wt(10),0.d0,p(1),-q(1))
      call xyzw(pt(1,11),wt(11),0.d0,-p(1),q(1))
      call xyzw(pt(1,12),wt(12),0.d0,-p(1),-q(1))
      call xyzw(pt(1,13),wt(13),q(1),p(1),0.d0,bigc(1))
      call xyzw(pt(1,14),wt(14),q(1),-p(1),0.d0,bigc(1))
      call xyzw(pt(1,15),wt(15),-q(1),p(1),0.d0,bigc(1))
      call xyzw(pt(1,16),wt(16),-q(1),-p(1),0.d0,bigc(1))
      call xyzw(pt(1,17),wt(17),q(1),0.d0,p(1),bigc(1))
      call xyzw(pt(1,18),wt(18),q(1),0.d0,-p(1),bigc(1))
      call xyzw(pt(1,19),wt(19),-q(1),0.d0,p(1),bigc(1))
      call xyzw(pt(1,20),wt(20),-q(1),0.d0,-p(1),bigc(1))            
      call xyzw(pt(1,21),wt(21),0.d0,q(1),p(1),bigc(1))
      call xyzw(pt(1,22),wt(22),0.d0,q(1),-p(1),bigc(1))
      call xyzw(pt(1,23),wt(23),0.d0,-q(1),p(1),bigc(1))
      call xyzw(pt(1,24),wt(24),0.d0,-q(1),-p(1),bigc(1)) 
      return
      end

