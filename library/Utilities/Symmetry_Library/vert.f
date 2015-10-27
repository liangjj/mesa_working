*deck @(#)vert.f	5.1  11/6/94
      subroutine vert(maxap3,n,natoms,nop,maxop,trans,nperm,
     $     a,b,d,lperm,iprmut)
      implicit real*8(a-h,o-z)
c
c     generate the n operations of a set of verticle planes.  one of
c     the planes is assumed to be coincident with the yz plane.
c
      dimension a(maxap3,3), b(maxap3,3), d(maxap3,3), q(3,3), r(3,3),
     $     s(3,3), t(3,3), trans(3,3,maxop), nperm(maxap3,maxop),
     $     lperm(1), iprmut(natoms,1)
      data half,one,eight/0.5d+00,1.0d+00,8.0d+00/
      save half,one,eight
c
c     call rtrace(6hvert  ,1)
      phi = eight * atan(one) / float(n)
      if (mod(n,2) .eq. 0) phi = half * phi
      theta = - phi
      do 100 iop=1,n
         theta = theta + phi
         call rotate(maxap3,a,b,natoms,q,3,theta)
         call reflct(maxap3,b,d,natoms,r,1)
         chi = - theta
         call rotate(maxap3,d,b,natoms,s,3,chi)
         call ebc(t,r,q,3,3,3)
         call ebc(q,s,t,3,3,3)
         call fill(maxap3,natoms,nop,maxop,q,trans,nperm,a,b,
     $        lperm,iprmut)
 100  continue
      return
      end
