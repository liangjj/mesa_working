*deck @(#)mcvcc.f	1.1  11/30/90
      subroutine mcvcc(aii,bii,aij,bij,bijt,deg,nsym,linear)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)mcvcc.f	1.1   11/30/90
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue
c
      implicit real*8(a-h,o-z)
      common / number / zero,pt5,one,two,four,eight
      dimension aii(10,2),bii(10,2),aij(10,10,2),bij(10,10,2)
      dimension bijt(10,10,2),deg(10)
c
      if(linear.ne.0) go to 300
c
      do 200 i=1,nsym
         deg(i)=one
         aii(i,1)=-one
         bii(i,1)=one+two
         do 210 j=1,nsym
            aij(i,j,1)=-two
            bij(i,j,1)=eight
            bijt(i,j,1)=-two
 210     continue
 200  continue
c
      return
c
 300  continue
c
      deg(1)=one
      aii(1,1)=-one
      bii(1,1)=one+two
      aij(1,1,1)=-two
      bij(1,1,1)=eight
      bijt(1,1,1)=-two
c
      do 410 j=2,nsym
         aij(j,1,1)=-four
         aij(j,1,2)=zero
         bij(j,1,1)=eight+eight
         bij(j,1,2)=zero
         bijt(j,1,1)=-four
         bijt(j,1,2)=zero
 410  continue
c
      if(nsym.eq.1) return
c
      do 500 i=2,nsym
         deg(i)=two
         aii(i,1)=-two
         aii(i,2)=-two
         bii(i,1)=eight+four+two
         bii(i,2)=-two
         do 510 j=2,nsym
            aij(i,j,1)=-four
            aij(i,j,2)=-four
            bij(i,j,1)=eight+eight+eight+eight
            bij(i,j,2)=zero
            bijt(i,j,1)=-four
            bijt(i,j,2)=-four
 510     continue
 500  continue
c
      return
      end
