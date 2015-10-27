*deck @(#)linear.f	5.1  11/6/94
      subroutine linear(b,g,s,mdim)
c
c***begin prologue     linear
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)linear.f	5.1   11/6/94
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       linear
c
      implicit real*8(a-h,o-z)
c
      dimension b(2),g(2),s(2)
c
      common /io/ inp,iout
c
      data zero/1.d-12/
      save zero
c
      if(mdim.eq.1)go to 200
c
c
c       this routine solves a set of simultaneous
c       linear equations via symmetric gaussian elimination
c       see  roothaan & bagus  methods of comp. physics  vol. ii
c
c
      imq=mdim+1
      mq=mdim-1
c
c         begin  elimination
c
      itot=mdim*(mdim+1)/2
      do 70 k=1,mq
 9971    format(/,3(2x,e15.7))
         imq=imq-1
         lmq=imq*(imq-1)/2
         adiag=1.d0/b(lmq+imq)
         gmq=g(imq)
         lend=imq-1
c
         do 60 l=1,lend
            lx=lend-l+1
            lnq=lx*(lx-1)/2
            diag=b(lmq+lx)*adiag
c
            do 65 m=1,lx
               b(lnq+m)=b(lnq+m)-diag*b(lmq+m)
 65         continue
            g(lx)=g(lx)-diag*gmq
 60      continue
 70   continue
c
c     end  elimination
c
c
c     begin   substitution
c
      ix=1
      s(1)=g(1)/b(1)
      do 110 i=2,mdim
         i1=i-1
         xx=0.d0
         do 100 k=1,i1
            xx=xx+s(k)*b(ix+k)
 100     continue
         ix=ix+i
         if(abs(b(ix)).lt.zero)go to 500
         s(i)=(g(i)-xx)/b(ix)
 110  continue
c
c        finished
c
      return
c
 200  continue
      s(1)=g(1)/b(1)
      return
c
 500  continue
      s(i)=0.d0
 600  format(//,'   singular matrix in linear ',d15.7)
      if(i.eq.mdim) return
      call lnkerr('singular matrix in linear ')
c
c
      end
