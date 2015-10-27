*deck @(#)ot2.f	5.1  11/6/94
      subroutine ot2(out,op)
c***begin prologue     ot2.f
c***date written       930715  (yymmdd)
c***revision date      11/6/94
c***keywords           O point group representation matrix T2 
c***author             RUSSO, thomas (lanl) (although I'm not sure I want to 
c                        admit it.
c***source             @(#)ot2.f	5.1  11/6/94
c***purpose            generate T2 representation matrix 
c***description        given a matrix which is a representation of an operator
c                      in the basis (x,y,z), generate a matrix in a basis
c                      (xz,yz,xy).  For the point group O this means generate
c                      a representation matrix of the operator in the irreducible
c                      representation T2.
c 
c          real*8 op   matrix of operator in basis (x,y,z)
c          real*8 out  output matrix of op in basis (xz,yz,xy)
c
c***references
c     mostly snagged out of trmat.  This is a real kludge, believe you me.
c
c***routines called
c***end prologue       ot2.f
      implicit none
      real*8 out(3,3),op(3,3)
      integer p1,p2,i,j,d,dum,food,ntot,ltot,mtot
      integer nxp(3),nyp(3),nzp(3),nxd(3),nyd(3),nzd(3),ds(3,2)
      integer inp,iout
      common /io/inp,iout
c
c the n*p arrays are just the n,l,m values of p functions from which the xz,yz,xy will be made
c
      data nxp/1,0,0/
      data nyp/0,1,0/
      data nzp/0,0,1/
c
c the n*d arrays are are the nlm values of xz,yz,xy
c
      data nxd/1,0,1/
      data nyd/0,1,1/
      data nzd/1,1,0/
c
c ds tells us that function 1 is xz, function 2 is yz...
c
      data ds/1,2,1, 3,3,2/

      call rzero(out,9)
c
c here's what we're gonna do.  For each product d=xz,yz,xy,  see how the
c pieces transform by looking at the elements of the matrix op
c This is snagged right out of trmat, but we don't bother with any x**2 elements
c or the like, as is done in trmat.
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
               out(d,food)=out(d,food)+op(p1,i)*op(p2,j)
 30         continue 
 20      continue 
 10   continue 
      return
      end
                  
