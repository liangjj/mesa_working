*deck vecprd.f
c***begin prologue     vecprd
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           copy
c***author             schneider, barry (nsf)
c***source             
c***purpose            vector outer product
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       vecprd
      subroutine vecprd(vout,va,vb,ni,nj,n,m,nc)
      implicit integer (a-z)
      real*8 vout, va, vb
      dimension vout(n,m,nc,nc), va(m), vb(n)
      common/io/inp, iout
      do 10 i=1,m
         do 20 j=1,n
            vout(j,i,ni,nj) = vout(j,i,ni,nj) + va(i)*vb(j)
 20      continue
 10   continue   
      return
      end       


