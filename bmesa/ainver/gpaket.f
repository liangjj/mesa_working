*deck gpaket.f
c***begin prologue     gpaket
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            calculate zero time gaussian wavepacket 
c***                   
c***references         
c
c***routines called    
c***end prologue       gpaket
      subroutine gpaket(u01,u02,u03,p1,p2,p3,eigc1,eigc2,eigc3,
     1                  wt1,wt2,wt3,psi0,n,nd,dim,coord)
      implicit integer (a-z)
      real*8 u01, u02, u03, p1, p2, p3, eigc1, eigc2, eigc3
      real*8 wt1, wt2, wt3, psi0
      real*8 z
      character*(*) coord
      character*80 title
      dimension nd(3)
      dimension p1(nd(1),*), p2(nd(2),*), p3(nd(3),*)
      dimension u01(nd(1),*), u02(nd(2),*), u03(nd(3),*)
      dimension eigc1(nd(1)), eigc2(nd(2)), eigc3(nd(3))
      dimension wt1(nd(1)), wt2(nd(2)), wt3(nd(3))
      dimension psi0(n)
      common/io/inp, iout
      pointer (p,z(1))
      if(dim.eq.1) then
         i1=1
         i2=i1+nd(1)
         need=wpadti(i2+nd(2))
      elseif(dim.eq.2) then
         i1=1
         i2=i1+nd(1)
         i3=i2+nd(2)
         i4=i3+nd(1)
         need=wpadti( i4 + nd(2) )
      elseif(dim.eq.3) then
         i1=1
         i2=i1+nd(1)
         i3=i2+nd(2)
         i4=i3+nd(3)
         i5=i4+nd(1)
         i6=i5+nd(2) 
         need=wpadti( i6 + nd(3) )
      endif  
      call memory(need,p,ngot,'gpaket',0)
      write(iout,1)
c     
c     first do projection of time independent packet on basis states.
c     then calculate projection onto unperturbed states by matrix
c     multiplication.
c
      if(coord.eq.'cartesian') then
         write(iout,2)
         if(dim.eq.1) then
            icnt=i1            
            do 10 i=1,n
               z(icnt) = p1(i,i)*wt1(i)*
     1                   exp(-eigc1(i)*eigc1(i))
               icnt=icnt+1
 10         continue
            call ebtc(z(i2),u01,z(i1),n,n,1)
            call copy(z(i2),psi0,n)
         elseif(dim.eq.2) then
            icnt=i1
            do 20 i=1,nd(1)
               z(icnt) = p1(i,i)*wt1(i)*exp(-eigc1(i)*eigc1(i))
               icnt=icnt+1
 20         continue
            icnt=i2   
            do 30 i=1,nd(2)
               z(icnt) = p2(i,i)*wt2(i)*exp(-eigc2(i)*eigc2(i))
               icnt=icnt+1
 30         continue   
            call ebtc(z(i3),u01,z(i1),nd(1),nd(1),1)
            call ebtc(z(i4),u02,z(i2),nd(2),nd(2),1)
            count=0
            icnt=i3
            do 40 i=1,nd(1)
               jcnt=i4
               do 50 j=1,nd(2)
                  count=count+1
                  psi0(count) = z(icnt)*z(jcnt)
                  jcnt=jcnt+1 
 50            continue
               icnt=icnt+1
 40         continue   
         elseif(dim.eq.3) then
            icnt=i1
            do 60 i=1,nd(1)
               z(icnt) = p1(i,i)*wt1(i)*exp(-eigc1(i)*eigc1(i))
               icnt=icnt+1
 60         continue
            icnt=i2   
            do 70 i=1,nd(2)
               z(icnt) = p2(i,i)*wt2(i)*exp(-eigc2(i)*eigc2(i))
               icnt=icnt+1
 70         continue
            icnt=i3               
            do 80 i=1,nd(3)
               z(icnt) = p3(i,i)*wt3(i)*exp(-eigc3(i)*eigc3(i))
               icnt=icnt+1
 80         continue               
            call ebtc(z(i4),u01,z(i1),nd(1),nd(1),1)
            call ebtc(z(i5),u02,z(i2),nd(2),nd(2),1)
            call ebtc(z(i6),u03,z(i3),nd(3),nd(3),1)
            count=0
            icnt=i4
            do 100 i=1,nd(1)
               jcnt=i5
               do 110 j=1,nd(2)
                  kcnt=i6
                  do 120 k=1,nd(3)
                     count=count+1
                     psi0(count)=z(icnt)*z(jcnt)*z(kcnt)
                     kcnt=kcnt+1
 120              continue   
                  jcnt=jcnt+1
 110           continue
               icnt=icnt+1
 100        continue   
         endif
      else
         write(iout,3)
         if(dim.eq.1) then
            icnt=i1
            do 200 i=1,nd(1)
               z(icnt) = eigc1(i)*p1(i,i)*wt1(i)*exp(-eigc1(i)*eigc1(i))
               icnt=icnt+1
 200        continue
            call ebtc(z(i2),u01,z(i1),n,n,1)        
            call copy(z(i2),psi0,n)
         elseif(dim.eq.2) then
            icnt=i1
            do 210 i=1,nd(1)
               z(icnt) = eigc1(i)*p1(i,i)*wt1(i)*exp(-eigc1(i)*eigc1(i))
               icnt=icnt+1
 210        continue
            jcnt=i2
            do 220 i=1,nd(2)
               z(jcnt) = p2(i,i)*wt2(i)
               jcnt=jcnt+1
 220        continue
            call ebtc(z(i3),u01,z(i1),nd(1),nd(1),1)
            call ebtc(z(i4),u02,z(i2),nd(2),nd(2),1)
            count=0
            icnt=i3
            do 230 i=1,nd(1)
               jcnt=i4
               do 240 j=1,nd(2)
                  count=count+1
                  psi0(count) = z(icnt)*z(jcnt)
                  jcnt=jcnt+1 
 240           continue
               icnt=icnt+1
 230        continue   
         elseif(dim.eq.3) then
            icnt=i1
            do 300 i=1,nd(1)
               z(icnt) = eigc1(i)*p1(i,i)*wt1(i)*exp(-eigc1(i)*eigc1(i))
               icnt=icnt+1
 300        continue
            icnt=i2
            do 310 i=1,nd(2)
               z(icnt) = p2(i,i)*wt2(i)
               icnt=icnt+1
 310        continue
            icnt=i3
            do 320 i=1,nd(3)
               z(icnt) = p3(i,i)*wt3(i)
               icnt=icnt+1              
 320        continue
            call ebtc(z(i4),u01,z(i1),nd(1),nd(1),1)
            call ebtc(z(i5),u02,z(i2),nd(2),nd(2),1)
            call ebtc(z(i6),u03,z(i3),nd(3),nd(3),1)
            count=0
            icnt=i4
            do 400 i=1,nd(1)
               jcnt=i5
               do 410 j=1,nd(2)
                  kcnt=i6
                  do 420 k=1,nd(3)
                     count=count+1
                     psi0(count)=z(icnt)*z(jcnt)*z(kcnt)
                     kcnt=kcnt+1
 420              continue   
                  jcnt=jcnt+1
 410           continue
               icnt=icnt+1
 400        continue   
         endif
      endif
      call memory(-ngot,p,idum,'gpaket',idum) 
      return
 1    format(/,5x,'initial wavepacket at t=0')
 2    format(/,5x,'the form of the radial packet is:',///,15x,
     1            'psi = exp(-r*r)')
 3    format(/,5x,'the form of the radial packet is:',///,15x,
     1            'psi = r * exp(-r*r)')
      end       
