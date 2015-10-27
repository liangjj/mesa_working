*deck  @(#)gtrans	2.1 10/10/91
      subroutine gtrans(alpha,center,pre,kx,ky,kz,ftran,ntri,n)
c***begin prologue     gtrans
c***date written       920324   (yymmdd)
c***revision date               (yymmdd)
c
c
c***keywords
c***author             schneider, barry  (nsf)
c***source             %W% %G% 
c***purpose
c                      fourier transform of gaussian 
c***description
c***references
c
c***routines called
c
c***end prologue       gtrans
c
      implicit integer (a-z)
      real *8 alpha, center, kx, ky, kz, ksq, pre, pre1, k12, k22
      complex *16 ftran, eye, comfac
      dimension alpha(ntri), center(3,ntri), pre(ntri), ftran(ntri,*)
      dimension kx(n), ky(n), kz(n)
      data eye /(0.d0,1.d0)/
c
      kcnt=0
      do 10 k1=1,n
         k12=kx(k1)*kx(k1)
         do 20 k2=1,n
            k22=ky(k2)*ky(k2)
            do 30 k3=1,n
               kcnt=kcnt+1
               ksq=k12+k22+kz(k3)*kz(k3)
               icnt=0
               do 40 i=1,nbf
                  do 50 j=1,i
                     icnt=icnt+1 
                     pre1=sqrt(1.d0/(2.d0*alpha(icnt)))
     1                                    /(2.d0*alpha(icnt))
                     comfac=eye*(kx(k1)*center(1,icnt)+
     1                           ky(k2)*center(2,icnt)+
     2                           kz(k3)*center(3,icnt))
                     comfac=exp(comfac)
                     ftran(icnt,kcnt)=pre(icnt)*pre1*comfac*
     1                                exp(-ksq/(4.d0*alpha(icnt)))
   50             continue
   40          continue
   30       continue      
   20    continue
   10 continue     
      return
      end
