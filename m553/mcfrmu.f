*deck @(#)mcfrmu.f	5.1  11/6/94
      subroutine mcfrmu(hold,u,del,lok,len,mix,nmo,nocc,l)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)mcfrmu.f	5.1   11/6/94
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue
c
      implicit real*8(a-h,o-z)
      dimension del(2),lok(2),len(2),mix(2)
      dimension hold(nmo,nmo),u(nmo,nmo)
c
      common /io/ inp,iout
c
c
      do 10 i=1,nmo
         do 10 j=1,nmo
            u(i,j)=0.0d0
            hold(i,j)=0.0d0
 10      continue
c
         l = 0
         do 20 i=1,nocc
            m1 = lok(i) + 1
            m2 = m1 + len(i) - 1
            if (m2 .lt. m1) go to 20
            do 22 m = m1, m2
               if (i .ge. mix(m)) go to 22
c      write (iout,2200) i, m, mix(m), l, del(l)
c 2200 format('0*mcfrmu i m mix l del ',4(1x,i6),1x,f16.8)
               l = l + 1
               u(mix(m),i) = del(l)
               u(i,mix(m)) = - del(l)
 22         continue
 20      continue
c
         call mcconv(u,u,nmo,nmo,hold)
         do 40 i=1,nmo
            do 50 j=1,nmo
               u(j,i)=u(j,i)+hold(j,i)*0.5d0
 50         continue
            u(i,i)=u(i,i)+1.0d0
 40      continue
c
c     write (iout,9000) u
c9000 format(' *mcfrmu u '//4(1x,f16.8))
         return
         end
