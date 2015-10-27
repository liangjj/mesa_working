      subroutine renorm(vl,vr,n,m)
      implicit integer(a-z)
      complex*16 vl, vr, cdotc, valc
      real*8 tst
      dimension vl(n,m), vr(n,m)
      common/io/inp, iout
      character*80 title
      do 10 i=1,m
         tst=1.d-15
         do 20 j=1,m
            tst=max(tst,dreal(vl(j,i)))
 20      continue
         if(tst.lt.1.d-15) then
            do 30 j=1,m
               tst=dimag(vl(j,i))
               vl(j,i)=tst
 30         continue   
         endif
 10   continue   
      do 40 i=1,m
         tst=1.d-15
         do 50 j=1,m
            tst=max(tst,dreal(vr(j,i)))
 50      continue
         if(tst.lt.1.d-15) then
            do 60 j=1,m
               tst=dimag(vr(j,i))
               vr(j,i)=tst
 60         continue   
         endif
 40   continue   
      do 70 i=1,m
         valc=cdotc(m,vl(1,i),1,vr(1,i),1)
c         write(iout,*) valc
         valc=1.d0/sqrt(valc)
c         write(iout,*) valc
         call cscal(m,valc,vr(1,i),1)
         call cscal(m,conjg(valc),vl(1,i),1)
 70   continue
      return
      end


