*deck @(#)nuctwo.f	5.1  11/6/94
      subroutine nuctwo(zan,c,nat,d2e,nd2e)
c
c***begin prologue     nuctwo
c***date written       861215  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords
c***author             saxe, paul (lanl)
c***source             @(#)nuctwo.f	5.1   11/6/94
c***purpose            calculation of nuclear repulsion second derivatives
c***description
c
c***references
c***routines called
c***end prologue       nuctwo
c
      implicit integer (a-z)
c
      real*8 zan(nat),c(3,nat),r,d2e(nd2e)
c
      common /io/     inp,iout
c
      iadd(i,j)=i*(i-1)/2+j
c
c     ----- calculate the repulsion energy -----
c
      call rzero(d2e,nd2e)
c
      do 10 i=1,nat
         do 9 j=1,i-1
c
            r=sqrt((c(1,i)-c(1,j))**2+(c(2,i)-c(2,j))**2+
     #             (c(3,i)-c(3,j))**2)
            if (r.lt.1.0d-06) then
               write (iout,1) i,j
    1          format (//,' ##### nuctwo: atoms',i3,' and',i3,' are',
     #                 ' too close',//)
               call lnkerr('atoms on top of each other')
            end if
c
            do 8 ic=1,3
               do 7 jc=1,3
                  ij=iadd(ic+3*(i-1),jc+3*(j-1))
                  ii=iadd(ic+3*(i-1),jc+3*(i-1))
                  jj=iadd(ic+3*(j-1),jc+3*(j-1))
                  if (ic.eq.jc) then
                     d2e(ii)=d2e(ii)-
     #                             zan(i)*zan(j)/r**3 +
     #                             3.0d+00*zan(i)*zan(j)*
     #                             (c(ic,i)-c(ic,j))**2/r**5
                     d2e(jj)=d2e(jj)-
     #                             zan(i)*zan(j)/r**3 +
     #                             3.0d+00*zan(i)*zan(j)*
     #                             (c(ic,i)-c(ic,j))**2/r**5
                     d2e(ij)=d2e(ij)+
     #                             zan(i)*zan(j)/r**3 -
     #                             3.0d+00*zan(i)*zan(j)*
     #                             (c(ic,i)-c(ic,j))**2/r**5
                  else
                     if (ic.ge.jc) then
                        d2e(ii)=d2e(ii)+
     #                             3.0d+00*zan(i)*zan(j)*
     #                             (c(ic,i)-c(ic,j))*
     #                             (c(jc,i)-c(jc,j))/r**5
                        d2e(jj)=d2e(jj)+
     #                             3.0d+00*zan(i)*zan(j)*
     #                             (c(ic,i)-c(ic,j))*
     #                             (c(jc,i)-c(jc,j))/r**5
                     end if
                     d2e(ij)=d2e(ij)-
     #                             3.0d+00*zan(i)*zan(j)*
     #                             (c(ic,i)-c(ic,j))*
     #                             (c(jc,i)-c(jc,j))/r**5
                  end if
    7          continue
    8       continue
    9    continue
   10 continue
c
c
      return
      end
