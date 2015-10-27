*deck vmat.f
c***begin prologue     vmat
c***date written       970797   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           3-d schroedinger equation
c***author             schneider, barry (nsf)
c***source             trap3d
c***purpose            interaction potential in dvr representation
c***                   
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       vmat
      subroutine vmat(eigc1,eigc2,eigc3,v,work,hbar,mass,scale,n,
     1                n1,n2,n3,numdim,prnt,ops)
      implicit integer (a-z)
      real*8 eigc1, eigc2, eigc3, v, work
      real*8 hbar, mass, scale, zij
      real*8 pi, zero, half, one, two, three, ten
      real*8 timau
      real*8 aij, sij, omegij, vij
      real*8 facij, tmp, fpkey
      character*32 pot, chrkey
      character*80 cpass
      character*320 card
      character*(*) ops
      character*1 itoc
      logical prnt, dollar, diag, odiag, toau, useau, logkey
      dimension eigc1(n1), eigc2(n2), eigc3(n3), v(n)
      dimension work(n,3)
      dimension aij(3,3), sij(3,3), omegij(3,3), vij(3,3), facij(3,3)
      dimension zij(3,3), tmp(3)
      common/io/inp, iout
      data pi/3.14159265358979323844d0/
      data zero, half, one, two, three, ten / 0.d0, .5d0, 1.d0, 2.d0, 
     1                                        3.d0, 10.d0 /
      data timau / 2.418884d-17 /
c
c
c
c     lets assume pairwise interactions of some sort
c
c        if the model potential is a harmonic oscillator the form is;
c
c                    m*omega*omega*r*r
c              v =   _________________
c                           2
c        if the model potential is an exponential the form is:
c
c        here the potential is
c             v = - a*exp(-b*r)
c
c        if the model potential is a coulomb potential the form is:
c        
c             v = z/r
c
c
      if ( dollar('$vint',card,cpass,inp) ) then
           pot=chrkey(card,'potential-type','none',' ')
           call rzero(vij,9)
           call rzero(omegij,9)
           call rzero(zij,9)
           call rzero(aij,9)
           call rzero(sij,9)
           scale=one 
           diag=logkey(card,'diagonal-potential',.false.,' ')
           odiag=logkey(card,'no-off-diagonal-potential',.false.,' ')
           if(diag) then
              do 10 i=1,numdim
                 vij(i,i)=fpkey(card,'v'//itoc(i)//itoc(i),zero,' ')
                 omegij(i,i)=fpkey(card,'omega'//itoc(i)//itoc(i),
     1                             zero,' ')
                 zij(i,i)=fpkey(card,'z'//itoc(i)//itoc(i),-1.d0,' ')
                 aij(i,i)=fpkey(card,'a'//itoc(i)//itoc(i),zero,' ')
                 sij(i,i)=fpkey(card,'s'//itoc(i)//itoc(i),zero,' ')
 10           continue
           endif
           if(odiag) then
              do 20 i=1,numdim                  
                 do 30 j=i+1,numdim
                    vij(i,j)=fpkey(card,'v'//itoc(i)//itoc(j),one,' ')
                    vij(j,i)=vij(i,j)
                    omegij(i,j)=fpkey(card,'omega'//itoc(i)//itoc(j),
     1                                ten,' ')
                    omegij(i,j)=omegij(i,j)*two*pi
                    omegij(j,i)=omegij(i,j)
                    zij(i,j)=fpkey(card,'z'//itoc(i)//itoc(j),zero,' ')
                    zij(j,i)=zij(i,j)
                    aij(i,j)=fpkey(card,'a'//itoc(i)//itoc(j),one,' ')
                    aij(j,i)=aij(i,j)
                    sij(i,j)=fpkey(card,'s'//itoc(i)//itoc(j),one,' ')
                    sij(i,j)=sij(j,i)
 30              continue
 20           continue
           endif                     
           if(pot.eq.'harmonic-oscillator') then
              if(numdim.eq.1) then
                 scale=one/(omegij(1,1)*hbar)
              endif
              if(numdim.eq.2) then
                 scale=one/(sqrt(omegij(1,1)*omegij(2,2))*hbar)
              endif
              if(numdim.eq.3) then
                 scale=omegij(1,1)*omegij(2,2)*omegij(3,3)
                 scale=scale**(one/three)
                 scale=one/(scale*hbar)
              endif
           endif             
      endif
      toau=logkey(ops,'to-atomic-units',.false.,' ')
      useau=logkey(ops,'use-atomic-units',.false.,' ')
      if(toau) then
         do 40 i=1,numdim
            do 50 j=1,i
               omegij(i,j)=omegij(i,j)*timau
               omegij(j,i)=omegij(i,j)
 50         continue
 40      continue   
      endif
      if(useau) then
         do 60 i=1,numdim
            do 70 j=1,i
               omegij(i,j)=one
               omegij(j,i)=omegij(i,j)
 70         continue
 60      continue   
      endif
      call rzero(v,n)
      if(pot.eq.'harmonic-oscillator') then
         do 80 i=1,numdim
            do 90 j=1,i
               facij(i,j)=mass*omegij(i,j)*omegij(i,j)*half
               facij(j,i)=facij(i,j)
 90         continue   
 80      continue   
         do 100 i=1,numdim
            if(diag) then
               write(iout,1) i, omegij(i,i)
            endif
            if(.not.odiag) then
               do 110 j=i+1,numdim
                  write(iout,2) i, j, omegij(i,j)
 110           continue   
            endif
 100     continue   
         if(diag) then    
            do 120 i=1,n1
               work(i,1)=facij(1,1)*eigc1(i)*eigc1(i)
 120        continue
            if(numdim.le.1) then
               call copy(work(1,1),v,n1)
               return
            endif
            do 130 i=1,n2
               work(i,2)=facij(2,2)*eigc2(i)*eigc2(i)      
 130        continue
            if(numdim.le.2) then
                call v2d(v,work,n,n1,n2)
            elseif(numdim.eq.3) then
                do 140 i=1,n3
                   work(i,3)=facij(3,3)*eigc3(i)*eigc3(i)      
 140            continue
                call v3d(v,work,n,n1,n2,n3)
            else
                call lnkerr('error in number of dimensions')
            endif
         endif
         if(.not.odiag) then    
            do 150 i=1,n1
               ii=n2*(i-1)
               do 160 j=1,n2
                  ij=ii+j
                  work(ij,1) = ( eigc1(i) - eigc2(j) ) *
     1                         ( eigc1(i) - eigc2(j) ) * facij(1,2)
 160           continue
 150        continue
            if(numdim.le.2) then
               call v2od(v,work,n,n1,n2)
            elseif(numdim.eq.3) then
               do 170 i=1,n1
                  ii=n3*(i-1)
                  do 180 j=1,n3
                     ij=ii+j
                     work(ij,2) = (eigc1(i) - eigc3(j)) * 
     1                            (eigc1(i) - eigc3(j)) * facij(1,3) 
 180              continue   
 170           continue
               do 190 i=1,n2
                  ii=n3*(i-1)
                  do 200 j=1,n3
                     ij=ii+j
                     work(ij,3) = (eigc2(i) - eigc3(j)) *
     1                            (eigc2(i) - eigc3(j)) *facij(2,3)
 200              continue
 190           continue
               call v3od(v,work,n,n1,n2,n3)
            else
                call lnkerr('error in number of dimensions')
            endif               
         endif   
      elseif(pot.eq.'exponential') then
         do 300 i=1,numdim
            if(diag) then
               write(iout,3) i, aij(i,i), sij(i,i)
            endif
            if(.not.odiag) then
               do 310 j=i+1,numdim
                  write(iout,4) i, j, aij(i,j), sij(i,j)
 310           continue   
            endif
 300     continue   
         if(diag) then
            do 320 i=1,n1
               work(i,1)=aij(1,1)*exp(-sij(1,1)*abs(eigc1(i))) 
 320        continue   
            if(numdim.le.1) then
               call copy(work(1,1),v,n1)
               return
            endif
            do 330 i=1,n2
               work(i,2)=aij(2,2)*exp(-sij(2,2)*abs(eigc2(i))) 
 330        continue
            if(numdim.le.2) then
                call v2d(v,work,n,n1,n2)
            elseif(numdim.eq.3) then
                do 340 i=1,n3
                   work(i,3)=aij(3,3)*exp(-sij(3,3)*abs(eigc3(i))) 
 340            continue 
                call v3d(v,work,n,n1,n2,n3)              
            else
                call lnkerr('error in number of dimensions')
            endif               
         endif
         if(.not.odiag) then
            do 350 i=1,n1
               ii=n2*(i-1)
               tmp(1)=abs(eigc1(i))
               do 360 j=1,n2    
                  ij=ii+j
                  tmp(2)=abs ( tmp(1) - abs(eigc2(j)) )
                  work(ij,1)=aij(1,2)*exp(-sij(1,2)*tmp(2))
 360           continue
 350        continue
            if(numdim.le.2) then
               call v2od(v,work,n,n1,n2)            
            elseif(numdim.eq.3) then
               do 370 i=1,n1
                  ii=n3*(i-1)
                  tmp(1)=abs(eigc1(i))
                  do 380 j=1,n3    
                     ij=ii+j
                     tmp(2)=abs ( tmp(1) - abs(eigc3(j)) )
                     work(ij,2)=aij(1,3)*exp(-sij(1,3)*tmp(2))
 380              continue
 370           continue
               do 390 i=1,n2
                  ii=n3*(i-1)
                  tmp(1)=abs(eigc2(i))
                  do 400 j=1,n3    
                     ij=ii+j
                     tmp(2)=abs ( tmp(1) - abs(eigc3(j)) )
                     work(ij,3)=aij(2,3)*exp(-sij(2,3)*tmp(2))
 400              continue
 390           continue
               call v3od(v,work,n,n1,n2,n3)    
            else
               call lnkerr('error in number of dimensions')    
            endif            
         endif            
      elseif(pot.eq.'well') then
         do 500 i=1,numdim
            if(diag) then
               write(iout,5) i, vij(i,i)
            endif
            if(.not.odiag) then
               do 510 j=i+1,numdim
                  write(iout,6) i, j, vij(i,j)
 510           continue   
            endif
 500     continue            
         if(diag) then
            do 520 i=1,n1
               work(i,1) = vij(1,1)
 520        continue                  
            if(numdim.le.1) then
               call copy(work(1,1),v,n1)
               return
            endif
            do 530 i=1,n2
               work(i,2) = vij(2,2)
 530        continue
            if(numdim.le.2) then
               call v2d(v,work,n,n1,n2)
            elseif(numdim.eq.3) then
               do 540 i=1,n3
                  work(i,3)= vij(3,3)
 540           continue   
               call v3d(v,work,n,n1,n2,n3)
            else
               call lnkerr('error in number of dimensions')    
            endif
         endif   
         if(.not.odiag) then
            do 600 i=1,n1
               ii=n2*(i-1)
               do 610 j=1,n2
                  ij=ii+j
                  work(ij,1) = vij(1,2)
 610           continue   
 600        continue
            if(numdim.le.2) then
               call v2od(v,work,n,n1,n2)
            elseif(numdim.eq.3) then
               do 620 i=1,n1
                  ii=n3*(i-1)
                  do 630 j=1,n3    
                     ij=ii+j
                     work(ij,2)=vij(1,3)
 630              continue
 620           continue
               do 640 i=1,n2
                  ii=n3*(i-1)
                  tmp(1)=abs(eigc2(i))
                  do 650 j=1,n3    
                     ij=ii+j
                     work(ij,3)=vij(2,3)
 650              continue
 640           continue
               call v3od(v,work,n,n1,n2,n3)    
            else
                call lnkerr('error in number of dimensions')
            endif
         endif            
      elseif(pot.eq.'coulomb') then
         do 700 i=1,numdim
            write(iout,7) i, zij(i,i)
 700     continue   
         do 710 i=1,n1
            work(i,1)=zij(1,1)/eigc1(i) 
 710     continue   
         if(numdim.le.1) then
            call copy(work(1,1),v,n1)
            return
         endif
         do 720 i=1,n2
            work(i,2)=zij(2,2)/eigc2(i) 
 720     continue
         if(numdim.le.2) then
            call v2d(v,work,n,n1,n2)
         elseif(numdim.eq.3) then
            do 730 i=1,n3
               work(i,3)=zij(3,3)/eigc3(i) 
 730        continue
            call v3d(v,work,n,n1,n2,n3)               
         else
            call lnkerr('error in number of dimensions')
         endif               
      elseif(pot.eq.'inverse-r4') then
         write(iout,8)
         do 800 i=1,n1
            work(i,1)=one/( 1.d0 + eigc1(i) )**4
 800     continue   
         if(numdim.le.1) then
            call copy(work(1,1),v,n1)
            return
         endif
         do 810 i=1,n2
            work(i,2)=one/( 1.d0 + eigc2(i) )**4
 810     continue
         if(numdim.le.2) then
            call v2d(v,work,n,n1,n2)
         elseif(numdim.eq.3) then
            do 820 i=1,n3
               work(i,3)=one/( 1.d0 + eigc3(i) )**4 
 820        continue
            call v3d(v,work,n,n1,n2,n3)               
         else
            call lnkerr('error in number of dimensions')
         endif               
      endif
      return
 1    format(/,5x,'diagonal frequency:',/,5x,
     1               '                  i = ',i1,/,5x,
     2               '                  omegaii = ',e15.8)
 2    format(/,5x,'off diagonal frequency:',/,5x,
     1           '                      i = ',i1,' j = ',i1,/,5x,
     2           '                      omegaij = ',e15.8)
 3    format(/,5x,'diagonal exponential constant:',/,5x,
     1            '                             i = ',i1,/,5x,
     2            '                             aii = ',e15.8,/,5x,
     3            '                             sii = ',e15.8)
 4    format(/,5x,'off diagonal exponential constants:',/,5x,
     1            '                                  i = ',i1,' j = '
     2                                                    ,i1,/,5x,
     3            '                                   aij = ',e15.8,
     4                                                       /,5x,
     5            '                                   sij = ',e15.8)
 5    format(/,5x,'diagonal well constant:',/,5x,
     1            '                      i = ',i1,/,5x,
     2            '                      vii = ',e15.8)
 6    format(/,5x,'off diagonal well constants:',/,5x,
     1            '                           i = ',i1,' j = ',i1,
     2                                                       /,5x, 
     3            '                           vij = ',e15.8)
 7    format(/,5x,'diagonal coulomb charges:',/,5x,
     1            '                        i = ',i1,/,5x,
     2            '                        zii = ',e15.8)
 8    format(/,5x,'diagonal 1/r4')
      end       
