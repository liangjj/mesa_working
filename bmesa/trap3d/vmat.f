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
      subroutine vmat(eigc1,eigc2,eigc3,v,work,hbar,mass,scale,n,nd,
     1                numdim,prnt,search,ops)
      implicit integer (a-z)
      real*8 eigc1, eigc2, eigc3, v, work
      real*8 hbar, mass, scale
      real*8 pi, zero, half, one, two, three, ten
      real*8 timau
      real*8 aij, sij, omegij, vij
      real*8 facij, tmp, fpkey
      character*32 pot1, chrkey
      character*8 cpass
      character*320 card
      character*(*) ops, search
      logical prnt, posinp, diag, odiag, toau, useau, logkey
      dimension nd(3)
      dimension eigc1(nd(1)), eigc2(nd(2)), eigc3(nd(3)), v(n)
      dimension work(n,3)
      dimension aij(3,3), sij(3,3), omegij(3,3), vij(3,3), facij(3,3)
      dimension tmp(3)
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
c        if the model potential is an constant the form is:
c        
c             v = - constant
c
c
      if(search.eq.'look for v0') then
         if( posinp('$v0',cpass) ) then
             call cardin(card)      
             pot1=chrkey(card,'unperturbed-potential-type','none',' ')
             scale=one
             diag=.true.
             odiag=.true.
             vij(1,1)=fpkey(card,'v11',zero,' ')
             omegij(1,1)=fpkey(card,'omega-11',zero,' ')
             omegij(1,1)=omegij(1,1)*two*pi
             if(pot1.eq.'harmonic-oscillator') then
                scale=one/(omegij(1,1)*hbar)
             endif   
             aij(1,1)=fpkey(card,'a11',zero,' ')
             sij(1,1)=fpkey(card,'s11',zero,' ')
         endif
      elseif(search.eq.'look for v1') then
         if( posinp('$v1',cpass) ) then
             call cardin(card)
             pot1=chrkey(card,'interaction-potential-type','none',' ')
             diag=logkey(card,'diagonal-potential',.false.,' ')
             odiag=logkey(card,'no-off-diagonal-potential',.false.,' ')
             vij(1,1)=fpkey(card,'v11',zero,' ')
             vij(1,2)=fpkey(card,'v12',one,' ')
             vij(1,3)=fpkey(card,'v13',zero,' ')
             vij(2,1)=vij(1,2)
             vij(2,2)=fpkey(card,'v22',zero,' ')
             vij(2,3)=fpkey(card,'v23',zero,' ')
             vij(3,1)=vij(1,3)
             vij(3,2)=vij(2,3)
             vij(3,3)=fpkey(card,'v33',zero,' ')
             omegij(1,1)=fpkey(card,'omega-11',zero,' ')
             omegij(1,1)=omegij(1,1)*two*pi
             omegij(1,2)=fpkey(card,'omega-12',ten,' ')
             omegij(1,2)=omegij(1,2)*two*pi
             omegij(1,3)=fpkey(card,'omega-13',ten,' ')
             omegij(1,3)=omegij(1,3)*two*pi
             omegij(2,1)=omegij(1,2)
             omegij(2,2)=fpkey(card,'omega-22',zero,' ')
             omegij(2,2)=omegij(2,2)*two*pi
             omegij(2,3)=fpkey(card,'omega-23',zero,' ')
             omegij(2,3)=omegij(2,3)*two*pi
             omegij(3,1)=omegij(1,3)
             omegij(3,2)=omegij(2,3)
             omegij(3,3)=fpkey(card,'omega-33',zero,' ')
             omegij(3,3)=omegij(3,3)*two*pi
             scale=one
             if(pot1.eq.'harmonic-oscillator') then
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
             aij(1,1)=fpkey(card,'a11',zero,' ')
             aij(1,2)=fpkey(card,'a12',one,' ')
             aij(1,3)=fpkey(card,'a13',one,' ')
             aij(2,1)=aij(1,2)
             aij(2,2)=fpkey(card,'a22',zero,' ')
             aij(2,3)=fpkey(card,'a23',zero,' ')
             aij(3,1)=aij(1,3)
             aij(3,2)=aij(2,3)
             aij(3,3)=fpkey(card,'a33',zero,' ')
             sij(1,1)=fpkey(card,'s11',zero,' ')
             sij(1,2)=fpkey(card,'s12',one,' ')
             sij(1,3)=fpkey(card,'s13',one,' ')
             sij(2,1)=sij(1,2)
             sij(2,2)=fpkey(card,'s22',zero,' ')
             sij(2,3)=fpkey(card,'s23',zero,' ')
             sij(3,1)=sij(1,3)
             sij(3,2)=sij(2,3)
             sij(3,3)=fpkey(card,'s33',zero,' ')
         endif      
      else
         call lnkerr('error in directive for potential')
      endif
      toau=logkey(ops,'to-atomic-units',.false.,' ')
      useau=logkey(ops,'use-atomic-units',.false.,' ')
      if(toau) then
         do 10 i=1,numdim
            do 20 j=1,i
               omegij(i,j)=omegij(i,j)*timau
               omegij(j,i)=omegij(i,j)
 20         continue
 10      continue   
      endif
      if(useau) then
         do 30 i=1,numdim
            do 40 j=1,i
               omegij(i,j)=one
               omegij(j,i)=omegij(i,j)
 40         continue
 30      continue   
      endif
      call rzero(v,n)
      if(pot1.eq.'harmonic-oscillator') then
         do 50 i=1,numdim
            do 60 j=1,i
               facij(i,j)=mass*omegij(i,j)*omegij(i,j)*half
               facij(j,i)=facij(i,j)
 60         continue   
 50      continue   
         do 70 i=1,numdim
            if(diag) then
               write(iout,1) i, omegij(i,i)
            endif
            if(.not.odiag) then
               do 80 j=i+1,numdim
                  write(iout,2) i, j, omegij(i,j)
 80            continue   
            endif
 70      continue   
         if(diag) then    
            do 90 i=1,nd(1)
               work(i,1)=facij(1,1)*eigc1(i)*eigc1(i)
 90         continue
            if(numdim.le.1) then
               call copy(work(1,1),v,nd(1))
               return
            endif
            do 100 i=1,nd(2)
               work(i,2)=facij(2,2)*eigc2(i)*eigc2(i)      
 100        continue
            if(numdim.le.2) then
                call v2d(v,work,n,nd)
            elseif(numdim.eq.3) then
                do 110 i=1,nd(3)
                   work(i,3)=facij(3,3)*eigc3(i)*eigc3(i)      
 110            continue
                call v3d(v,work,n,nd)
            else
                call lnkerr('error in number of dimensions')
            endif
         endif
         if(.not.odiag) then    
            do 120 i=1,nd(1)
               ii=nd(2)*(i-1)
               do 130 j=1,nd(2)
                  ij=ii+j
                  work(ij,1) = ( eigc1(i) - eigc2(j) ) *
     1                         ( eigc1(i) - eigc2(j) ) * facij(1,2)
 130           continue
 120        continue
            if(numdim.le.2) then
               call v2od(v,work,n,nd)
            elseif(numdim.eq.3) then
               do 140 i=1,nd(1)
                  ii=nd(3)*(i-1)
                  do 150 j=1,nd(3)
                     ij=ii+j
                     work(ij,2) = (eigc1(i) - eigc3(j)) * 
     1                            (eigc1(i) - eigc3(j)) * facij(1,3) 
 150              continue   
 140           continue
               do 160 i=1,nd(2)
                  ii=nd(3)*(i-1)
                  do 170 j=1,nd(3)
                     ij=ii+j
                     work(ij,3) = (eigc2(i) - eigc3(j)) *
     1                            (eigc2(i) - eigc3(j)) *facij(2,3)
 170              continue
 160           continue
               call v3od(v,work,n,nd)
            else
                call lnkerr('error in number of dimensions')
            endif               
         endif   
      elseif(pot1.eq.'exponential') then
         do 200 i=1,numdim
            if(diag) then
               write(iout,3) i, aij(i,i), sij(i,i)
            endif
            if(.not.odiag) then
               do 210 j=i+1,numdim
                  write(iout,4) i, j, aij(i,j), sij(i,j)
 210           continue   
            endif
 200     continue   
         if(diag) then
            do 220 i=1,nd(1)
               work(i,1)=aij(1,1)*exp(-sij(1,1)*abs(eigc1(i))) 
 220        continue   
            if(numdim.le.1) then
               call copy(work(1,1),v,nd(1))
               return
            endif
            do 230 i=1,nd(2)
               work(i,2)=aij(2,2)*exp(-sij(2,2)*abs(eigc2(i))) 
 230        continue
            if(numdim.le.2) then
                call v2d(v,work,n,nd)
            elseif(numdim.eq.3) then
                do 240 i=1,nd(3)
                   work(i,3)=aij(3,3)*exp(-sij(3,3)*abs(eigc3(i))) 
 240            continue               
            else
                call lnkerr('error in number of dimensions')
            endif               
         endif
         if(.not.odiag) then
            do 250 i=1,nd(1)
               ii=nd(2)*(i-1)
               tmp(1)=abs(eigc1(i))
               do 260 j=1,nd(2)    
                  ij=ii+j
                  tmp(2)=abs ( tmp(1) - abs(eigc2(j)) )
                  work(ij,1)=aij(1,2)*exp(-sij(1,2)*tmp(2))
 260           continue
 250        continue
            if(numdim.le.2) then
               call v2od(v,work,n,nd)            
            elseif(numdim.eq.3) then
               do 270 i=1,nd(1)
                  ii=nd(3)*(i-1)
                  tmp(1)=abs(eigc1(i))
                  do 280 j=1,nd(3)    
                     ij=ii+j
                     tmp(2)=abs ( tmp(1) - abs(eigc3(j)) )
                     work(ij,2)=aij(1,3)*exp(-sij(1,3)*tmp(2))
 280              continue
 270           continue
               do 290 i=1,nd(2)
                  ii=nd(3)*(i-1)
                  tmp(1)=abs(eigc2(i))
                  do 300 j=1,nd(3)    
                     ij=ii+j
                     tmp(2)=abs ( tmp(1) - abs(eigc3(j)) )
                     work(ij,3)=aij(2,3)*exp(-sij(2,3)*tmp(2))
 300              continue
 290           continue
               call v3od(v,work,n,nd)    
            else
               call lnkerr('error in number of dimensions')    
            endif            
         endif            
      elseif(pot1.eq.'well') then
         do 400 i=1,numdim
            if(diag) then
               write(iout,5) i, vij(i,i)
            endif
            if(.not.odiag) then
               do 410 j=i+1,numdim
                  write(iout,6) i, j, vij(i,j)
 410           continue   
            endif
 400     continue            
         if(diag) then
            do 420 i=1,nd(1)
               work(i,1) = vij(1,1)
 420        continue                  
            if(numdim.le.1) then
               call copy(work(1,1),v,nd(1))
               return
            endif
            do 430 i=1,nd(2)
               work(i,2) = vij(2,2)
 430        continue
            if(numdim.le.2) then
               call v2d(v,work,n,nd)
            elseif(numdim.eq.3) then
               do 440 i=1,nd(3)
                  work(i,3)= vij(3,3)
 440           continue   
               call v3d(v,work,n,nd)
            else
               call lnkerr('error in number of dimensions')    
            endif
         endif   
         if(.not.odiag) then
            do 500 i=1,nd(1)
               ii=nd(2)*(i-1)
               do 510 j=1,nd(2)
                  ij=ii+j
                  work(ij,1) = vij(1,2)
 510           continue   
 500        continue
            if(numdim.le.2) then
               call v2od(v,work,n,nd)
            elseif(numdim.eq.3) then
               do 520 i=1,nd(1)
                  ii=nd(3)*(i-1)
                  do 530 j=1,nd(3)    
                     ij=ii+j
                     work(ij,2)=vij(1,3)
 530              continue
 520           continue
               do 540 i=1,nd(2)
                  ii=nd(3)*(i-1)
                  tmp(1)=abs(eigc2(i))
                  do 550 j=1,nd(3)    
                     ij=ii+j
                     work(ij,3)=vij(2,3)
 550              continue
 540           continue
               call v3od(v,work,n,nd)    
            else
                call lnkerr('error in number of dimensions')
            endif                
         endif
      endif
      return
 1    format(/,5x,'diagonal frequency:',/,5x,
     1               '                  i = ',i1,/,5x,
     2               '                  omegaii = ',e15.8)
 2       format(/,5x,'off diagonal frequency:',/,5x,
     1               '                      i = ',i1,' j = ',i1,/,5x,
     2               '                      omegaij = ',e15.8)
 3       format(/,5x,'diagonal exponential constant:',/,5x,
     1               '                             i = ',i1,/,5x,
     2               '                             aii = ',e15.8,/,5x,
     3               '                             sii = ',e15.8)
 4       format(/,5x,'off diagonal exponential constants:',/,5x,
     1               '                                  i = ',i1,' j = '
     2                                                       ,i1,/,5x,
     3               '                                   aij = ',e15.8,
     4          /,5x,'                                   sij = ',e15.8)
 5       format(/,5x,'diagonal well constant:',/,5x,
     1               '                      i = ',i1,/,5x,
     2               '                      vii = ',e15.8)
 6       format(/,5x,'off diagonal well constants:',/,5x,
     1               '                           i = ',i1,' j = ',i1,
     2                                                          /,5x, 
     3               '                           vij = ',e15.8)
      end       
