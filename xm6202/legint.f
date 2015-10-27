*deck legint.f
c***begin prologue     legint
c***date written       930922   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           legint, link m6201, legendre integrals
c***author             schneider, barry (nsf)
c***source             m6201
c***purpose            calculate overlap integrals of legendre functions
c***references         none
c
c***routines called
c***end prologue       legint
      subroutine legint (plm,phifn,scr,nth,nph,lmax,m,nm)
      implicit integer (a-z)
      real*8 plm, phifn, scr, phint
      character*80 title
      character*2 itoc
      dimension plm(nth,m:lmax), phifn(nph,nm)
      dimension scr(m:lmax,m:lmax)
      common /io/ inp, iout
      dim=lmax-m+1
      call rzero(scr,dim*dim)
      do 10 l1=m,lmax
         do 20 l2=l1,lmax
            do 30 n=1,nth
               scr(l1,l2)=scr(l1,l2) + plm(n,l1)*plm(n,l2)
   30       continue
            scr(l2,l1)=scr(l1,l2)          
   20    continue
   10 continue 
      title='legendre overlap integrals for m = '//itoc(m)
      call prntfm(title,scr,dim,dim,dim,dim,iout)
      do 40 number=1,nm
         phint=0.d0
         do 50 n=1,nph
            phint = phint + phifn(n,number)*phifn(n,number)
   50    continue
         write (iout,1) m, phint
   40 continue  
    1 format(/,5x,'phi overlap integral for m = ',i4,1x,e15.8)    
      return
      end

