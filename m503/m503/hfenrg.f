*deck  @(#)hfenrg.f	5.1 11/6/94
      subroutine hfenrg(jmat,kmat,dmat,nnp,nbf,ncoul,nexch,ndmat,h,
     $                  enuc,energy,c,gamma,calc,eonel,etwoel,
     $                  nshell,fcoef,alpha,beta)
c
c***begin prologue     hfenrg
c***date written       870521   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           scf energy
c***author             saxe, paul (lanl)
c***source             @(#)hfenrg.f	5.1 11/6/94 
c
c***purpose            to calculate the energy.
c
c***description
c
c***references
c
c***routines called    sdot
c
c***end prologue       hfenrg
c
      implicit integer (a-z)
c
      real*8 jmat(nnp,ncoul),kmat(nnp,nexch),dmat(nnp,ndmat)
      real*8 h(nnp),c(2),gamma,eonel,etwoel
      real*8 enuc,energy,two,half,four,one,three,l,w
      real*8 haa,hbb,hab
      real*8 fcoef(nshell),alpha(nshell,nshell),beta(nshell,nshell)
      character*(*) calc
c
      real*8 sdot
c
      parameter (two=2.0d+00, half=1.0d+00/two,four=4.0d+00)
      parameter (one=1.0d+00 ,three=3.0d+00)
c
c     ----- double the off-diagonal part of the density matrices
c           to account for symmetry
c
      do 3 i=1,ndmat
         do 1 j=1,nnp
            dmat(j,i)=dmat(j,i)*two
    1    continue
         do 2 j=1,nbf
            jj=j*(j+1)/2
            dmat(jj,i)=dmat(jj,i)*half
    2    continue
    3 continue
c
      if (calc.eq.'closed') then
c
c        ----- closed shell scf -----
c            e= 2*h(ii) + 2*j(ii,jj) - k(ii,jj)
c
         eonel=two*sdot(nnp,h,1,dmat,1)
         etwoel=two*sdot(nnp,jmat,1,dmat,1)-
     $              sdot(nnp,kmat,1,dmat,1)
c
      else if (calc.eq.'open') then
c
c        ----- high-spin open-shell scf -----
c
c        e = 2 * h(ii) +  2  * j(ii,jj) -       k(ii,jj)
c              + h(aa) + 1/2 * j(aa,bb) - 1/2 * k(aa,bb)
c              +          2  * j(ii,aa) -       k(ii,aa)
c
c          i,j = d.o. ; a,b = h.o
c
         eonel=two*sdot(nnp,h,1,dmat(1,1),1)+
     $             sdot(nnp,h,1,dmat(1,2),1)
         etwoel= two*sdot(nnp,jmat(1,1),1,dmat(1,1),1)-
     $               sdot(nnp,kmat(1,1),1,dmat(1,1),1)+
     $          half*sdot(nnp,jmat(1,2),1,dmat(1,2),1)-
     $          half*sdot(nnp,kmat(1,2),1,dmat(1,2),1)+
     $           two*sdot(nnp,jmat(1,1),1,dmat(1,2),1)-
     $               sdot(nnp,kmat(1,1),1,dmat(1,2),1)
c
      else if (calc.eq.'gvb') then
c
c        ----- tcscf -----
c
c        haa = 2 * h(ii) + 2 * j(ii,jj) -     k(ii,jj)
c            + 2 * h(aa) + 2 * j(aa,aa) -     k(aa,aa)
c            +             4 * j(ii,aa) - 2 * k(ii,aa)
c
c        hbb = 2 * h(ii) + 2 * j(ii,jj) -     k(ii,jj)
c            + 2 * h(bb) + 2 * j(bb,bb) -     k(bb,bb)
c            +             4 * j(ii,bb) - 2 * k(ii,bb)
c
c        hab = k(aa,bb)
c
c        l=hbb-haa
c
c        w=1/sqrt(l**2+hab**2)
c
c        c1**2 = 1/2(1+lw)
c
c        c2**2 = 1/2(1-lw)   (c2 is < 0)
c
c        gamma = (3-2*c1*c2) * k(aa,bb)
c              - (1+2*c1*c2) * (j(aa,bb) - j(aa,aa) + k(aa,aa))
c
         haa=two*sdot(nnp,h,1,dmat(1,1),1)+
     $       two*sdot(nnp,jmat(1,1),1,dmat(1,1),1)-
     $           sdot(nnp,kmat(1,1),1,dmat(1,1),1)+
     $       two*sdot(nnp,h,1,dmat(1,2),1)+
     $       two*sdot(nnp,jmat(1,2),1,dmat(1,2),1)-
     $           sdot(nnp,kmat(1,2),1,dmat(1,2),1)+
     $      four*sdot(nnp,jmat(1,1),1,dmat(1,2),1)-
     $       two*sdot(nnp,kmat(1,1),1,dmat(1,2),1)
c
         hbb=two*sdot(nnp,h,1,dmat(1,1),1)+
     $       two*sdot(nnp,jmat(1,1),1,dmat(1,1),1)-
     $           sdot(nnp,kmat(1,1),1,dmat(1,1),1)+
     $       two*sdot(nnp,h,1,dmat(1,3),1)+
     $       two*sdot(nnp,jmat(1,3),1,dmat(1,3),1)-
     $           sdot(nnp,kmat(1,3),1,dmat(1,3),1)+
     $      four*sdot(nnp,jmat(1,1),1,dmat(1,3),1)-
     $       two*sdot(nnp,kmat(1,1),1,dmat(1,3),1)
c
         hab=sdot(nnp,kmat(1,2),1,dmat(1,3),1)
c
         l=hbb-haa
         w=one/sqrt(l**2+four*hab**2)
c
         c(1)= sqrt(half*(1+l*w))
         c(2)=-sqrt(half*(1-l*w))
c
         gamma=(three-two*c(1)*c(2))*hab-
     $         (one  +two*c(1)*c(2))*(sdot(nnp,jmat(1,2),1,dmat(1,3),1)-
     $                                sdot(nnp,jmat(1,2),1,dmat(1,2),1)+
     $                                sdot(nnp,kmat(1,2),1,dmat(1,2),1))
c
c        e = c1**2 * haa  +  c2**2 * hbb  +  2 * c1 * c2 * hab
c
         eonel=c(1)**2*(two*sdot(nnp,h,1,dmat(1,1),1)+
     $                  two*sdot(nnp,h,1,dmat(1,2),1))+
     $         c(2)**2*(two*sdot(nnp,h,1,dmat(1,1),1)+
     $                  two*sdot(nnp,h,1,dmat(1,3),1))
         etwoel=c(1)**2*haa+c(2)**2*hbb+two*c(1)*c(2)*hab-eonel
c
      else if(calc.eq.'general') then
c
         eonel=0.d0
         etwoel=0.d0
         do 999 i=1,nshell-1
            eonel=eonel+two*fcoef(i)*sdot(nnp,h,1,dmat(1,i),1)
            etwoel=etwoel+alpha(i,i)*sdot(nnp,jmat(1,i),1,dmat(1,i),1)+
     $                   beta(i,i)*sdot(nnp,kmat(1,i),1,dmat(1,i),1)
 999     continue
c
         if(nshell.gt.2) then
            do 998 i=2,nshell-1
               do 997 j=1,i-1
                  etwoel= etwoel
     $                 +alpha(i,j)*sdot(nnp,jmat(1,i),1,dmat(1,j),1)+
     $                  alpha(j,i)*sdot(nnp,jmat(1,j),1,dmat(1,i),1)+
     $                  beta(i,j)*sdot(nnp,kmat(1,i),1,dmat(1,j),1)+
     $                  beta(j,i)*sdot(nnp,kmat(1,j),1,dmat(1,i),1)
 997           continue
 998        continue
         end if
      end if
c
c
      energy=enuc+eonel+etwoel
c
c     ----- halve the off-diagonal part of the density matrices
c           to retrun them to normal
c
      do 6 i=1,ndmat
         do 4 j=1,nnp
            dmat(j,i)=dmat(j,i)*half
    4    continue
         do 5 j=1,nbf
            jj=j*(j+1)/2
            dmat(jj,i)=dmat(jj,i)*two
    5    continue
    6 continue
c
c
      return
      end
