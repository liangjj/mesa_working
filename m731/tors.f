*deck @(#)tors.f	5.1  11/6/94
      subroutine tors(noint,i,j,k,l,b,ib,c)
c***begin prologue     tors.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             binkley, et al. (g82)
c***source             @(#)tors.f	5.1   11/6/94
c***purpose            
c***description
c
c        adapted from the normal coordinate analysis program of
c        schachtschneider, shell development .
c
c***references
c
c***routines called
c
c***end prologue       tors.f
      implicit none
c     --- input variables -----
      integer noint,i,j,k,l
c     --- input arrays (unmodified) ---
      real*8 c(*)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      integer ib(4,2)
      real*8 b(3,4,2)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer iaind,jaind,laind,kaind,m
      integer inp,iout
      real*8 rij(3),rjk(3),rkl(3),eij(3),ejk(3),ekl(3),cr1(3),cr2(3)
      real*8 zero,one,dijsq,djksq,dklsq,dij,djk,dkl
      real*8 dotpj,dotpk,sinpj,sinpk,smi,smj,sml,f1,f2
      parameter (zero=0.0d+00,one=1.0d+00)
c
      common/io/inp,iout
c
c
      iaind=3*(i-1)
      jaind=3*(j-1)
      kaind=3*(k-1)
      laind=3*(l-1)
      ib(1,noint)=i
      ib(2,noint)=j
      ib(3,noint)=k
      ib(4,noint)=l
      dijsq=zero
      djksq=zero
      dklsq=zero
      do 10 m=1,3
         rij(m)=c(m+jaind)-c(m+iaind)
         dijsq=dijsq+rij(m)**2
         rjk(m)=c(m+kaind)-c(m+jaind)
         djksq=djksq+rjk(m)**2
         rkl(m)=c(m+laind)-c(m+kaind)
         dklsq=dklsq+rkl(m)**2
   10 continue
      dij=sqrt(dijsq)
      djk=sqrt(djksq)
      dkl=sqrt(dklsq)
c
      do 20 m=1,3
         eij(m)=rij(m)/dij
         ejk(m)=rjk(m)/djk
         ekl(m)=rkl(m)/dkl
   20 continue
c
      cr1(1)=eij(2)*ejk(3)-eij(3)*ejk(2)
      cr1(2)=eij(3)*ejk(1)-eij(1)*ejk(3)
      cr1(3)=eij(1)*ejk(2)-eij(2)*ejk(1)
      cr2(1)=ejk(2)*ekl(3)-ejk(3)*ekl(2)
      cr2(2)=ejk(3)*ekl(1)-ejk(1)*ekl(3)
      cr2(3)=ejk(1)*ekl(2)-ejk(2)*ekl(1)
      dotpj=-(eij(1)*ejk(1)+eij(2)*ejk(2)+eij(3)*ejk(3))
      dotpk=-(ejk(1)*ekl(1)+ejk(2)*ekl(2)+ejk(3)*ekl(3))
      sinpj=sqrt(one-dotpj**2)
      sinpk=sqrt(one-dotpk**2)
c
      do 30 m=1,3
         smi=-cr1(m)/(dij*sinpj*sinpj)
         b(m,1,noint)=smi
         f1=(cr1(m)*(djk-dij*dotpj))/(djk*dij*sinpj*sinpj)
         f2=(dotpk*cr2(m))/(djk*sinpk*sinpk)
         smj=f1-f2
         b(m,2,noint)=smj
         sml= cr2(m)/(dkl*sinpk*sinpk)
         b(m,4,noint)=sml
         b(m,3,noint)=(-smi-smj-sml)
   30 continue
c
c
      return
      end
