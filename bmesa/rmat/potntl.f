*deck potntl
c***begin prologue     potntl
c***date written       940213   (yymmdd)
c***revision date               (yymmdd)
c***keywords           potential, matrix
c***author             schneider, barry (nsf)
c***source             
c***purpose            potential matrix elements
c***description        
c***references       
c
c***routines called
c***end prologue       potntl
      subroutine potntl(fi,fj,vcij,pot,vij,pt,range,ni,nj,npts,
     1                  type,chni,chnj,ntri,prnt)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 fi, fj, vcij, vij, pot, pt, range
      character*(*) type
      dimension fi(npts,ni), fj(npts,nj)
      dimension vcij(ntri), vij(ni,nj), pot(npts), pt(npts)
      logical prnt
      index=chni*(chni-1)/2+chnj
      call filpot(vcij(index),pot,pt,range,npts,chni,chnj,type)
      do 10 i=1,ni
         do 20 j=1,nj
            vij(i,j)=0.d0
            do 30 k=1,npts
               vij(i,j)=vij(i,j)+fi(k,i)*pot(k)*fj(k,j)
   30       continue
   20    continue
   10 continue
      return
      end



