c $Header $
*deck legint.f
c***begin prologue     legint
c***date written       921221   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           legint, link m6201, legendre integrals
c***author             schneider, barry (nsf)
c***source             m6201
c***purpose            calculate overlap integrals of legendre functions
c***references         none
c
c***routines called
c***end prologue       legint
      subroutine legint (plm,phifn,scr,scr1,nth,nph,num,lmax,m)
      implicit integer (a-z)
      real*8 plm, phifn, scr, scr1
      character*80 title
      character*2 itoc
      dimension plm(nth,m:lmax), phifn(nph,num)
      dimension scr(m:lmax,m:lmax), scr1(num,num)
      common /io/ inp, iout
      call rzero(scr,(lmax+1)*(lmax+1))
      do 10 l1=m,lmax
         do 20 l2=l1,lmax
            do 30 n=1,nth
               scr(l1,l2)=scr(l1,l2) + plm(n,l1)*plm(n,l2)
   30       continue
            scr(l2,l1)=scr(l1,l2)          
   20    continue
   10 continue 
      title='legendre overlap integrals for m = '//itoc(m)
      dim=lmax-m+1
      call prntfm(title,scr,dim,dim,dim,dim,iout)
      do 40 i=1,num
         do 50 j=1,i
            scr1(i,j)=0.d0
            do 60 n=1,nph
               scr1(i,j)=scr1(i,j) + phifn(n,i)*phifn(n,j)
   60       continue
            scr1(j,i)=scr1(i,j)
   50    continue
   40 continue
      title='phi overlap integrals for m = '//itoc(m)
      call prntfm(title,scr1,num,num,num,num,iout)
      return
      end

