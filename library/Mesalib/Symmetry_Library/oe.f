*deck @(#)oe.f	5.1  11/6/94
      subroutine oe(out,op)
c***begin prologue     oe.f
c***date written       930715  (yymmdd)
c***revision date      11/6/94
c***keywords           O point group representation matrix E
c***author             RUSSO, thomas (lanl) (although I'm not sure I want to 
c                        admit it.
c***source             @(#)oe.f	5.1  11/6/94
c***purpose            generate E representation matrix 
c***description        given a matrix which is a representation of an operator
c                      in the basis (x,y,z), generate a matrix in a basis
c                      (3z**2-r**2,x**2-y**2).
c                      For the point group O this means generate
c                      a representation matrix of the operator in the
c                      irreducible representation E
c 
c          real*8 op   matrix of operator in basis (x,y,z)
c          real*8 out  output matrix of op in basis (3z**2-r**2,x**2-y**2)
c
c***references
c     mostly snagged out of trmat.  This is a real kludge, believe you me.
c     the first part of this is snagged out of ot2
c
c***routines called
c***end prologue       oe.f
c
      implicit none
      real*8 out(2,2),op(3,3),temp(3,3),tmp2(3,2)
      real*8 a(2,3),b(3,2)
      integer p1,p2,i,j,d,dum,food,ntot,ltot,mtot
      integer nxp(3),nyp(3),nzp(3),nxd(3),nyd(3),nzd(3),ds(3,2)
      integer inp,iout
      common /io/inp,iout
c
c the n*p arrays are just the n,l,m values of p functions from which the products
c will be made
c
      data nxp/1,0,0/
      data nyp/0,1,0/
      data nzp/0,0,1/
c
c the n*d arrays are are the nlm values of xz,yz,xy
c
      data nxd/2,0,0/
      data nyd/0,2,0/
      data nzd/0,0,2/
c
c ds tells us that function 1 is x**2, function 2 is y**...
c
      data ds/1,2,3, 1,2,3/
c
c A is part of the (x**2,y**2,z**2)->(3z**2-r**2,x**2-y**2,r**2) matrix
c
      a(1,1)=-1.0d0
      a(1,2)=-1.0d0
      a(1,3)=2.0d0
      a(2,1)=1.0d0
      a(2,2)=-1.0d0
      a(2,3)=0.0d0
c
c B is part of the inverse of the above
c
      b(1,1)=-1.0d0/6.0d0
      b(2,1)=-1.0d0/6.0d0
      b(3,1)=1.0d0/3.0d0
      b(1,2)=0.5d0
      b(2,2)=-0.5d0
      b(3,2)=0.0d0

      call rzero(temp,9)
c
c this one is a two stepper.  We generate the representation in basis
c x**2,y**2,z**2, which is reducible into E and A1.  Then we reduce it with 
c a change of basis.
c
      do 10 d=1,3
         p1=ds(d,1)
         p2=ds(d,2)
         do 20 i=1,3
            do 30 j=1,3
               ntot=nxp(i)+nxp(j)
               ltot=nyp(i)+nyp(j)
               mtot=nzp(i)+nzp(j)
               food=4
               do 40 dum=1,3
                  if (nxd(dum).eq.ntot .and.
     $                nyd(dum).eq.ltot .and.
     $                nzd(dum).eq.mtot) then
                     food=dum
                     goto 41
                  endif
 40            continue 
 41            if (food.eq.4) goto 30
               temp(d,food)=temp(d,food)+op(p1,i)*op(p2,j)
 30         continue 
 20      continue 
 10   continue 
c
c now we have the 3x3 matrix in basis x**2,y**2,z**2.  
c to reduce this representation we need to transform to a basis 
c (3z**2-r**2,x**2-y**2,r**2).  The matrix which takes the latter into the
c former is
c
c   -1  -1  2                        -1/6  .5  1/3
c A= 1  -1  0   and its inverse is B=-1/6 -.5  1/3
c    1   1  1                         1/3   0  1/3
c but since we only want the E part of the transformed matrix we can 
c just throw away the last row of the former and the last column of the latter
c Doing AopB gives us the desired rep.
      call ebc(tmp2,temp,b,3,3,2)
      call ebc(out,a,tmp2,2,3,2)
      return
      end
                  
