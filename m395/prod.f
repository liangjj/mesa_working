*deck  @(#)prod.f	2.1 10/10/91
      subroutine prod(alpha,center,alf12,cen12,pre12,nbf,ntri)
c***begin prologue     prod
c***date written       920324   (yymmdd)
c***revision date               (yymmdd)
c
c
c***keywords
c***author             schneider, barry  (nsf)
c***source             %W% %G% 
c***purpose
c                      product of two gaussians 
c***description
c***references
c
c***routines called
c
c***end prologue       prod
c
      implicit integer (a-z)
      real *8 alpha, center, alf12, cen12, alfrat, pre12
      dimension alpha(nbf), center(3,nbf), alf12(ntri), cen12(3,ntri)
      dimension pre12(ntri)
c
      icnt=0
      do 10 i=1,nbf
         do 20 j=1,i
            icnt=icnt+1
            alf12(icnt)=alpha(i)+alpha(j)
            alfrat=alpha(i)*alpha(j)/alf12(icnt)
            pre12(icnt)=0.d0 
            do 30 k=1,3
               cen12(k,icnt)=(alpha(i)*center(k,i)+
     1                        alpha(j)*center(k,j))/alf12(icnt)
               pre12(icnt)=pre12(icnt)+(center(k,i)-center(k,j))
     1                                *(center(k,i)-center(k,j))
   30       continue
            pre12(icnt)=exp(-pre12(icnt)*alfrat)
   20    continue
   10 continue     
      return
      end

