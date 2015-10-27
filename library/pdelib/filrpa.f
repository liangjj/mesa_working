*deck filrpa.f
c***begin prologue     filrpa
c***date written       960615   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           rpa hamiltonian
c***author             schneider, barry (nsf)
c***source             
c***purpose            calculate and store matrix elements
c***                   of rpa hamiltonian in the transformed representation.
c***
c***description        the transformed representation uses two matrices
c***                   [ A + B ] and [ A - B ] where the original eigenvalue
c***                   problem is written as;
c***
c***                  [  A  B ] [X]   [ E  0 ] [X]
c***                  [       ]     = [      ]
c***                  [ -B -A ] [Y]   [ 0  E ] [Y]      
c***                   
c***references         
c
c***routines called    
c***end prologue       filrpa
      subroutine filrpa(apb,amb,hbufa,diaga,ibufa,hbufb,ibufb,
     1                  e,eps,kappa,q,lenbuf,itdiag,nrpa,
     2                  neapb,neamb,prnt)
      implicit integer (a-z)
      real*8 apb, amb, aplusb, aminusb
      real*8 hbufa, diaga, hbufb
      real*8 e, eps, kappa, q
      character*80 title
      logical itdiag, prnt
      dimension apb(nrpa,*), amb(nrpa,*), ibufa(2,*), ibufb(2,*)
      dimension hbufa(*), diaga(*), hbufb(*), q(nrpa), e(nrpa)
      common/io/inp, iout
      if(itdiag) then
         title='calculate and store non-zero matrix elements of rpa '//
     1         'hamiltonian'
         call iosys('create integer "a plus b matrix buffers" on ham',
     1               -1,0,0,' ')
         write(iout,1) title
         count=0
         neapb=0
         do 10 i=1,nrpa
            do 20 j=1,i-1
               aplusb=2.d0*kappa*q(i)*q(j)
               if(aplusb.ne.0.d0) then
                  count=count+1
                  if(count.gt.lenbuf) then
                     neapb=neapb+lenbuf
                     call iosys('write integer "a plus b matrix '//
     1                          'buffers" to ham without rewinding',
     2                           2*lenbuf,ibufa,0,' ')
                     call iosys('write integer "a plus b matrix '//
     1                          'buffers" to ham without rewinding',
     2                           wptoin(lenbuf),hbufa,0,' ')
                     count=1
                  endif
                  ibufa(1,count)=i
                  ibufa(2,count)=j
                  hbufa(count)=aplusb
               endif
 20         continue
            do 30 j=i+1,nrpa
               aplusb=2.d0*kappa*q(i)*q(j) 
               if(aplusb.ne.0.d0) then
                  count=count+1
                  if(count.gt.lenbuf) then
                     neapb=neapb+lenbuf
                     call iosys('write integer "a plus b matrix '//
     1                          'buffers" to ham without rewinding',
     2                           2*lenbuf,ibufa,0,' ')
                     call iosys('write integer "a plus b matrix '//
     1                          'buffers to ham without rewinding',
     2                           wptoin(lenbuf),hbufa,0,' ')
                     count=1
                  endif               
                  ibufa(1,count)=i
                  ibufa(2,count)=j
                  hbufa(count)=aplusb
               endif
 30         continue
            diaga(i) = eps*(i-1) + 2.d0*kappa*q(i)*q(i)
 10      continue
         if(count.gt.0) then
            neapb=neapb+count
            call iosys('write integer "a plus b matrix '//
     1                 'buffers" to ham without rewinding',
     2                  2*count,ibufa,0,' ')
            call iosys('write integer "a plus b matrix '//
     1                 'buffers" to ham without rewinding',
     2                  wptoin(count),hbufa,0,' ')
         endif
         call iosys('endfile "a plus b matrix buffers" on ham',
     1               0,0,0,' ')
         call iosys('write integer "number of a plus b matrix '//
     1              'elements" to ham',1,neapb,0,' ')
         call iosys('write real "a plus b matrix diagonals" to ham',
     1               nrpa,diaga,0,' ')
         write(iout,2) neapb
         if(prnt) then
            title='packed matrix'
            write(iout,3)
            write(iout,4)
            do 1000 i=1,nrpa
               write(iout,5) i, diaga(i)
 1000       continue
            if(count.ne.0) then
               write(iout,6)
               write(iout,7)
               do 2000 i=1,neapb
                  write(iout,8) ibufa(1,i),ibufa(2,i), hbufa(i)
 2000          continue
            endif                        
         endif
         call iosys('create integer "a minus b matrix buffers" on ham',
     1               -1,0,0,' ')
         count=0
         neamb=0
         do 40 i=1,nrpa
            do 50 j=1,nrpa
c               aminusb = kappa*q(i)*q(ii)
               aminusb =0.d0
               if(i.eq.j) then
                  aminusb = eps*(i-1)
               endif   
               if(aminusb.ne.0.d0) then
                  count=count+1
                  if(count.gt.lenbuf) then
                     neamb=neamb+lenbuf
                     call iosys('write integer "a minus b matrix '//
     1                          'buffers" to ham without rewinding',
     2                           2*lenbuf,ibufb,0,' ')
                     call iosys('write integer "a minus b matrix '//
     1                          'buffers to ham without rewinding',
     2                           wptoin(lenbuf),hbufb,0,' ')
                     count=1
                  endif               
                  ibufb(1,count)=i
                  ibufb(2,count)=j
                  hbufb(count)=aminusb
               endif                             
 50         continue
 40      continue
         if(count.gt.0) then
            neamb=neamb+count
            call iosys('write integer "a minus b matrix '//
     1                 'buffers" to ham without rewinding',
     2                  2*count,ibufb,0,' ')
            call iosys('write integer "a minus b matrix '//
     1                 'buffers" to ham without rewinding',
     2                  wptoin(count),hbufb,0,' ')
         endif
         do 60 i=1,nrpa
            e(i)=eps*eps*(i-1) +2.d0*eps*(i-1)*q(i)*q(i)*kappa
 60      continue   
         write(iout,9) neamb
         if(prnt) then
            title='packed matrix'
            if(count.ne.0) then
               write(iout,11)
               write(iout,7)
               do 3000 i=1,count
                  write(iout,8) ibufb(1,i),ibufb(2,i), hbufb(i)
 3000          continue
            endif          
         endif
         call iosys('endfile "a minus b matrix buffers" on ham',
     1               0,0,0,' ')
         call iosys('write integer "number of a minus b matrix '//
     1              'elements" to ham',1,neapb,0,' ')
      else
         title='calculate  rpa hamiltonian explicitly'
         write(iout,1) title
         do 70 i=1,nrpa
            do 80 j=1,nrpa
               apb(i,j) = 2.d0*kappa*q(i)*q(j)
 80         continue
            apb(i,i) = apb(i,i) + eps*(i-1)
 70      continue
         do 90 i=1,nrpa
            do 100 j=1,nrpa
               amb(i,j)=0.d0
               if(i.eq.j) then
                  amb(i,j) = eps*(i-1)
               endif   
 100        continue
 90      continue
         if(prnt) then
            title='A plus B rpa hamiltonian'
            call prntrm(title,apb,nrpa,nrpa,nrpa,nrpa,iout)
            title= 'A minus B rpa hamiltonian'
            call prntrm(title,amb,nrpa,nrpa,nrpa,nrpa,iout)                   
         endif
      endif         
      return
 1    format(a80)     
 2    format(/,1x,'number non-zero,off-diagonal a matrix elements = ',
     1         i6)
 3    format(/,1x,'a plus b diagonal hamiltonian matrix elements')
 4    format(/,5x,'  index ',4x,'  matrix element  ')
 5    format(5x,i5,6x,e15.8)
 6    format(/,1x,'non-zero off diagonal a plus b matrix elements')
 7    format(/,5x,'    i   ',4x,'   j   ',4x,'  matrix element  ')
 8    format(5x,i5,6x,i5,6x,e15.8)
 9    format(/,1x,'number non-zero,off-diagonal a minus b matrix '
     1            'elements = ',i6)
 11   format(/,1x,'non-zero off diagonal a minus b matrix elements')
      end       







