*deck hampp.f
c***begin prologue     hampp
c***date written       960615   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           hamppiltonian
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            out of core matrix vector multiply
c***                   using non-zero elements of matrix.
c***description                           
c***                   
c***references         
c
c***routines called    
c***end prologue       hampp
      subroutine hampp(ibuf,rbuf,diag,hpp,pvec,spp,energy,
     #                 header,nwksp,nwksq,nwks,npvec,lenbuf,
     #                 nonz,nodsk,nen,prnt)
      implicit integer (a-z)
      real*8 rbuf, diag, hpp, spp, tmp, pvec, energy
      character*80 title
      character*16 fptoc
      character*(*) header(7)
      logical nodsk, prnt
      dimension rbuf(lenbuf), ibuf(2,lenbuf), diag(*)
      dimension hpp(npvec,npvec)
      dimension pvec(nwks,npvec)
      dimension spp(nwksp,npvec)
      dimension energy(nen)
      common/io/inp, iout 
      data null / 1.d-14 /
      call rzero(spp,nwksp*npvec)
      if(nodsk) then
         call rdham(ibuf,rbuf,header(1),lenbuf,nonz) 
         call hmulpp(ibuf,rbuf,pvec,spp,
     #               nwks,nwksp,nwksq,npvec,nonz,lenbuf)
      else
         trips = nonz/lenbuf
         left = nonz - trips*lenbuf
         do 1000 trp=1,trips
            call iosys('read integer '//header(1)//' from '//
     #                 'hamiltonian without rewinding',
     #                  2*lenbuf,ibuf,0,' ')
            call iosys('read integer '//header(1)//' from '//
     #                 'hamiltonian without rewinding',
     #                  wptoin(lenbuf),rbuf,0,' ')
            call hmulpp(ibuf,rbuf,pvec,spp,
     #                   nwks,nwksp,nwksq,npvec,lenbuf,lenbuf)
 1000    continue
         if(left.ne.0) then
            call iosys('read integer '//header(1)//' from '//
     #                 'hamiltonian without rewinding',
     #                  2*left,ibuf,0,' ')
            call iosys('read integer '//header(1)//' from '//
     #                 'hamiltonian without rewinding',
     #                  wptoin(left),rbuf,0,' ')
            call hmulpp(ibuf,rbuf,pvec,spp,
     #                  nwks,nwksp,nwksq,npvec,left,lenbuf)
         end if
      end if
      count=nwksq
      do 2000 i=1,nwksp
         count=count+1
         tmp=diag(count)
         do 3000 j=1,npvec
            spp(i,j) = spp(i,j) + tmp * pvec(count,j)
 3000    continue
 2000 continue
      call ebtcxx(hpp,pvec(nwksq+1,1),spp,
     #            npvec,nwksp,npvec,npvec,nwks,nwksp)
      do 4000 i=1,npvec
         diag(i) = hpp(i,i)
 4000 continue   
      call iosys('write real '//header(2)//' to '//
     #           'hamiltonian',npvec,diag,0,' ')
      call iosys('create integer '//header(3)//
     #           ' on hamiltonian',-1,0,0,' ')
      n=0
      nonzpp=0
      do 5000 i=1,npvec
         do 6000 j=1,i-1
            if(abs(hpp(i,j)).gt.null) then
               n=n+1
               if(n.gt.lenbuf) then
                  nonzpp=nonzpp+lenbuf
                  call iosys('write integer '//header(3)//
     #                       ' to hamiltonian '//
     #                       'without rewinding',2*lenbuf,ibuf,0,' ') 
                  call iosys('write integer '//header(3)//
     #                       ' to hamiltonian '//
     #                       'without rewinding',wptoin(lenbuf),
     #                        rbuf,0,' ')        
                  n=1
               endif
               ibuf(1,n)=i
               ibuf(2,n)=j
               rbuf(n)=hpp(i,j)
            endif
 6000    continue         
 5000 continue
      if (n.gt.0) then
          nonzpp=nonzpp+n
          call iosys('write integer '//header(3)//
     #               ' to hamiltonian without rewinding',
     #                 2*n,ibuf,0,' ') 
          call iosys('write integer '//header(3)//
     #               ' to hamiltonian without rewinding',
     #                 wptoin(n),rbuf,0,' ')        
      endif
      call iosys('endfile '//header(3)//' on hamiltonian',0,0,0,' ')
      call iosys('write integer '//header(4)//' to hamiltonian',
     #            1,nonzpp,0,' ')
      write(iout,1) nonzpp
      if(prnt) then
         title='Contracted H_PP'
         call prntrm(title,hpp,npvec,npvec,npvec,npvec,iout)
      endif
      call iosys('create real '//header(7)//' on hamiltonian',
     #            nen*npvec*npvec,0,0,' ')
      do 7000 ien=1,nen
         do 8000 i=1,npvec
            hpp(i,i)=diag(i)-energy(ien)
 8000    continue   
         call iosys('write real '//header(7)//' to hamiltonian '//
     #              'without rewinding',npvec*npvec,hpp,0,' ')
         if(prnt) then
            title='static-exchange hamiltonian energy = '
     #            //fptoc(energy(ien))
            call prntrm(title,hpp,npvec,npvec,npvec,npvec,iout)
          endif
 7000 continue   
      return
 1    format(/,1x,'number of non-zero contracted p-space '
     1            'matrix elements = ',i5 )
      end       


















