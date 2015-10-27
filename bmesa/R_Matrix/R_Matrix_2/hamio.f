c*deck hamio.f
c***begin prologue     hamio
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            in(out)put hamiltonian information.
c***                   
c***references         
c
c***routines called    
c***end prologue       hamio
      subroutine hamio(pham,type,inout,hon,sym,n,nc,prn)
      implicit integer (a-z)
      integer*8 pham 
      real*8 ham
      character*(*) inout, type, sym
      character*80 title
      logical prn, hon
      common/io/inp, iout
      pointer (pham,ham(1))
      h=1
      eig=h+n*n
      if(inout.eq.'output') then
         write(iout,1) n
         call iosys('write integer "'//type//' size" to ham',
     1               1,n,0,' ')                      
         call iosys('write real "'//type//' eigenvalues" to '//
     1              'ham',n,ham(eig),0,' ')                      
         call iosys('write real "'//type//' eigenvectors" to '//
     1              'ham',n*n,ham(h),0,' ')
         if(prn) then
            title=type//' eigenvalues'
            call prntfm(title,ham(eig),n,1,n,1,iout)
            title=type//' eigenvectors' 
            call prntfm(title,ham(h),n,n,n,n,iout)
         endif
      elseif(inout.eq.'input') then
         call iosys('read integer "'//type//' size" from ham',
     1               1,n,0,' ')                      
         call iosys('read integer "number of channels" from ham',
     1               1,nc,0,' ')   
         write(iout,2) n, nc         
         eig=1
         gamma=eig+n
         eig0=gamma+n*nc 
         need=eig0+nc
         if(hon) then
            h=need
            need=need+n*n
         endif       
         need=wpadti(need)
         call getmem(need,pham,ngot,'hamin',0)
         call iosys('read real "'//type//' eigenvalues" from '//
     1              'ham',n,ham(eig),0,' ')                      
         if(sym.eq.'symmetric') then
            call iosys('read real "channel energies" from ham',
     1                  nc,ham(eig0),0,' ') 
         else
            call iosys('read integer "number of x channels" from ham',
     1                  1,nx,0,' ')
            call iosys('read real "x channel energies" from ham',
     1                  nx,ham(eig0),0,' ')
            call iosys('read integer "number of y channels" from ham',
     1                  1,ny,0,' ')
            call iosys('read real "y channel energies" from ham',
     1                  ny,ham(eig0+nx),0,' ')
         endif
         write(iout,*) 'got here' 
         call iosys('read real "h amplitudes" from ham',
     1               nc*n,ham(gamma),0,' ') 
         if(hon) then
            call iosys('read real "'//type//' eigenvectors" from '//
     1                 'ham',n*n,ham(h),0,' ')
         endif
         if(prn) then
            title=type//' eigenvalues'
            call prntfm(title,ham(eig),n,1,n,1,iout)
            title='r-matrix amplitudes'
            call prntfm(title,ham(gamma),n,nc,n,nc,iout)
            title='channel eigenvalues'
            call prntfm(title,ham(eig0),nc,1,nc,1,iout)
            if(hon) then
               title=type//' eigenvectors' 
               call prntfm(title,ham(h),n,n,n,n,iout)
            endif
         endif
      endif
      return      
 1    format(/,1x,'writing hamiltonian information for n = ',i6,
     1            ' size matrix')
 2    format(/,1x,'reading hamiltonian information for n = ',i6,
     1            ' size matrix',/,1x,'number of channels = ',i3)
      end       






