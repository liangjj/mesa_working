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
      subroutine gpaket(p1,p2,p3,eigc1,eigc2,eigc3,
     1                  wt1,wt2,wt3,psi0,v1,v2,v3,t,
     2                  n,n1,n2,n3,dim,coord,phase)
      implicit integer (a-z)
      real*8 p1, p2, p3, eigc1, eigc2, eigc3
      real*8 wt1, wt2, wt3, psi0
      complex*16 v1, v2, v3, t, cfac
      character*(*) coord
      character*80 str
      logical phase
      dimension p1(n1,*), p2(n2,*), p3(n3,*)
      dimension eigc1(n1), eigc2(n2), eigc3(n3)
      dimension wt1(n1), wt2(n2), wt3(n3)
      dimension psi0(n,2), v1(*), v2(*), v3(*)
      dimension t(*)
      common/io/inp, iout
      if(coord.eq.'cartesian') then
         str='the initial wavepacket is cartesian'
         write(iout,1) str
         if(phase) then
            write(iout,2)
         else
            write(iout,3)
         endif
         if(dim.eq.1) then
            call cfactr(p1,eigc1,wt1,v1,t,n,phase)
            do 10 i=1,n
               psi0(i,1) = real(v1(i))
               psi0(i,2) = imag(v1(i))
 10         continue                      
         elseif(dim.eq.2) then
            call cfactr(p1,eigc1,wt1,v1,t,n1,phase)
            call cfactr(p2,eigc2,wt2,v2,t,n2,phase)
            count=0
            do 20 i=1,n1
               do 30 j=1,n2
                  count=count+1
                  cfac = v1(i)*v2(j)
                  psi0(count,1) = real(cfac)
                  psi0(count,2) = imag(cfac)
 30            continue
 20         continue   
         elseif(dim.eq.3) then
            call cfactr(p1,eigc1,wt1,v1,t,n1,phase)
            call cfactr(p2,eigc2,wt2,v2,t,n2,phase)
            call cfactr(p3,eigc3,wt3,v3,t,n3,phase)
            count=0
            do 40 i=1,n1
               do 50 j=1,n2
                  do 60 k=1,n3
                     cfac = v1(i)*v2(j)*v3(k)
                     count=count+1
                     psi0(count,1) = real(cfac)
                     psi0(count,2) = imag(cfac)
 60               continue
 50            continue
 40         continue   
         endif
      else
         str='the initial wavepacket is a radial packet'
         write(iout,1) str
         if(phase) then
            write(iout,4)
         else
            write(iout,5)
         endif
         if(dim.eq.1) then
            call rfactr(p1,eigc1,wt1,v1,t,n1,phase)
            do 100 i=1,n
               psi0(i,1) = real(v1(i))
               psi0(i,2) = imag(v1(i))
 100        continue                        
         elseif(dim.eq.2) then
            call rfactr(p1,eigc1,wt1,v1,t,n1,phase)
            call rfactr(p2,eigc2,wt2,v2,t,n2,phase)
            count=0
            do 200 i=1,n1
               do 210 j=1,n2
                  count=count+1
                  cfac = v1(i)*v2(j)
                  psi0(count,1) = real(cfac)
                  psi0(count,2) = imag(cfac)                  
 210           continue
 200        continue   
         elseif(dim.eq.3) then
            call rfactr(p1,eigc1,wt1,v1,t,n1,phase)
            call rfactr(p2,eigc2,wt2,v2,t,n2,phase)
            call rfactr(p3,eigc3,wt3,v3,t,n3,phase)
            count=0
            do 300 i=1,n1
               do 310 j=1,n2
                  do 320 k=1,n3
                     count=count+1
                     cfac = v1(i)*v2(j)*v3(k)
                     psi0(count,1) = real(cfac)
                     psi0(count,2) = imag(cfac)
 320              continue
 310           continue   
 300        continue       
         endif
      endif
      return
 1    format(/,5x,a80)
 2    format(/,5x,'the form of the packet is:',///,15x,
     1            'psi = exp(-icos(r)-r*r)')
 3    format(/,5x,'the form of the packet is:',///,15x,
     1            'psi = exp(-r*r)')
 4    format(/,5x,'the form of the packet is:',///,15x,
     1            'psi = r * exp(-icos(r)-r*r)')
 5    format(/,5x,'the form of the packet is:',///,15x,
     1            'psi = r * exp(-r*r)')
       end       





