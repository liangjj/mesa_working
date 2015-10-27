*deck @(#)dprim.f	5.1  11/6/94
      subroutine dprim(nroots,nv,rhotsq,urho,rho,dc00,dcp00,alpha,a,b,
     #                 d1exp,c,camcb,nat,dercen,ndcen,npass,
     #                 nmax,mmax,d2exp,nd2,nderiv)
c
c***begin prologue     dprim
c***date written       850101   (yymmdd)
c***revision date      871011   (yymmdd)
c
c   11 november 1987   pws at lanl
c      adding 'd2exp' for the second derivatives.
c
c***keywords
c***author             saxe, paul (lanl)
c***source             @(#)dprim.f	5.1   11/6/94
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       dprim
c
      implicit integer (a-z)
c
      real*8 rhotsq(nroots,nv),urho(nroots,nv),rho(nv)
      real*8 dc00(nroots,nv,ndcen),dcp00(nroots,nv,ndcen),alpha(nv,4)
      real*8 a(nv),b(nv),d1exp(nroots,nv,3,ndcen),c(3,nat)
      real*8 camcb(nv,3),cimcj,ckmcl
      real*8 d2exp(nroots,nv,nd2)
      integer dercen(4)
c
      do 2 root=1,nroots
         do 1 i=1,nv
            rhotsq(root,i)=urho(root,i)*rho(i)/(urho(root,i)+rho(i))
    1    continue
    2 continue
c
      if (nmax.gt.0) then
         do 4 root=1,nroots
            do 3 i=1,nv
               dc00(root,i,1)=alpha(i,1)/a(i)*
     #            (1.0d+00-rhotsq(root,i)/a(i))-1.0d+00
    3       continue
    4    continue
c
         if (npass.le.2) then
            do 6 root=1,nroots
               do 5 i=1,nv
                  dc00(root,i,2)=alpha(i,2)/a(i)*
     #               (1.0d+00-rhotsq(root,i)/a(i))
    5          continue
    6       continue
         else if (npass.eq.3) then
            do 8 root=1,nroots
               do 7 i=1,nv
                  dc00(root,i,2)=alpha(i,3)/b(i)*rhotsq(root,i)/a(i)
    7          continue
    8       continue
         end if
c
         if (npass.eq.1) then
            do 10 root=1,nroots
               do 9 i=1,nv
                  dc00(root,i,3)=alpha(i,3)/b(i)*rhotsq(root,i)/a(i)
    9          continue
   10       continue
         end if
      end if
c
c
      if (mmax.gt.0) then
         do 14 root=1,nroots
            do 13 i=1,nv
               dcp00(root,i,1)=alpha(i,1)/a(i)*rhotsq(root,i)/b(i)
   13       continue
   14    continue
c
         if (npass.le.2) then
            do 16 root=1,nroots
               do 15 i=1,nv
                  dcp00(root,i,2)=alpha(i,2)/a(i)*rhotsq(root,i)/b(i)
   15          continue
   16       continue
         else if (npass.eq.3) then
            do 18 root=1,nroots
               do 17 i=1,nv
                  dcp00(root,i,2)=alpha(i,3)/b(i)*
     #               (1.0d+00-rhotsq(root,i)/b(i))-1.0d+00
   17          continue
   18       continue
         end if
c
         if (npass.eq.1) then
            do 20 root=1,nroots
               do 19 i=1,nv
                  dcp00(root,i,3)=alpha(i,3)/b(i)*
     #               (1.0d+00-rhotsq(root,i)/b(i))-1.0d+00
   19          continue
   20       continue
         end if
      end if
c
c     ----- form the derivatives of the exponential term -----
c
      do 30 coord=1,3
         cimcj=c(coord,dercen(1))-c(coord,dercen(2))
         ckmcl=c(coord,dercen(3))-c(coord,dercen(4))
         do 22 root=1,nroots
            do 21 i=1,nv
               d1exp(root,i,coord,1)=(-cimcj*alpha(i,2)-
     #             camcb(i,coord)*rhotsq(root,i))*alpha(i,1)/a(i)*
     #                2.0d+00
   21       continue
   22    continue
c
         if (npass.le.2) then
            do 24 root=1,nroots
               do 23 i=1,nv
                  d1exp(root,i,coord,2)=(cimcj*alpha(i,1)-
     #                camcb(i,coord)*rhotsq(root,i))*alpha(i,2)/a(i)*
     #                   2.0d+00
   23          continue
   24       continue
         else if (npass.eq.3) then
            do 26 root=1,nroots
               do 25 i=1,nv
                  d1exp(root,i,coord,2)=(-ckmcl*alpha(i,4)+
     #                camcb(i,coord)*rhotsq(root,i))*alpha(i,3)/b(i)*
     #                   2.0d+00
   25          continue
   26       continue
         end if
         if (npass.eq.1) then
            do 28 root=1,nroots
               do 27 i=1,nv
                  d1exp(root,i,coord,3)=(-ckmcl*alpha(i,4)+
     #                camcb(i,coord)*rhotsq(root,i))*alpha(i,3)/b(i)*
     #                   2.0d+00
   27          continue
   28       continue
         end if
   30 continue
c
c     ----- and the second-derivatives of the exponential if needed
c
      if (nderiv.ge.2) then
         do 32 root=1,nroots
            do 31 i=1,nv
               d2exp(root,i,1)=(-alpha(i,1)*alpha(i,2)/a(i)-
     $              alpha(i,1)**2*rhotsq(root,i)/a(i)**2)*2.0d+00
 31         continue
 32      continue
c
         if (npass.le.2) then
            do 34 root=1,nroots
               do 33 i=1,nv
                  d2exp(root,i,2)=(alpha(i,1)*alpha(i,2)/a(i)-
     $                 alpha(i,1)*alpha(i,2)*rhotsq(root,i)/a(i)**2)*
     $                 2.0d+00
                  d2exp(root,i,3)=(-alpha(i,1)*alpha(i,2)/a(i)-
     $                 alpha(i,2)**2*rhotsq(root,i)/a(i)**2)*2.0d+00
 33            continue
 34         continue
         else if (npass.eq.3) then
            do 36 root=1,nroots
               do 35 i=1,nv
                  d2exp(root,i,2)=2.0d+00*alpha(i,1)*alpha(i,3)*
     $                 rhotsq(root,i)/(a(i)*b(i))
                  d2exp(root,i,3)=(-alpha(i,3)*alpha(i,4)/b(i)-
     $                 alpha(i,3)**2*rhotsq(root,i)/b(i)**2)*2.0d+00
 35            continue
 36         continue
         end if
         if (npass.eq.1) then
            do 38 root=1,nroots
               do 37 i=1,nv
                  d2exp(root,i,4)=2.0d+00*alpha(i,1)*alpha(i,3)*
     $                 rhotsq(root,i)/(a(i)*b(i))
                  d2exp(root,i,5)=2.0d+00*alpha(i,2)*alpha(i,3)*
     $                 rhotsq(root,i)/(a(i)*b(i))
                  d2exp(root,i,6)=(-alpha(i,3)*alpha(i,4)/b(i)-
     $                 alpha(i,3)**2*rhotsq(root,i)/b(i)**2)*2.0d+00
 37            continue
 38         continue
         end if
      end if
c
c
      return
      end
