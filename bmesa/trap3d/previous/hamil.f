*deck hamil.f
c***begin prologue     hamil
c***date written       960615   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           hamiltonian
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            calculate and store non-zero matrix elements
c***                   of hamiltonian and their indices in dvr
c***                   representation.
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       hamil
      subroutine hamil(h01,h02,h03,v,hbuf,ibuf,ind,diag,lenbuf,
     1                 nd,dim,n,prnt,ntot,incore,title)
      implicit integer (a-z)
      real*8 h01, h02, h03
      real*8 hbuf, diag, ham, v, dum
      character*(*) title
      character*1 itoc
      logical prnt, incore
      dimension nd(3)
      dimension h01(nd(1),nd(1)), h02(nd(2),nd(2)), h03(nd(3),nd(3))
      dimension hbuf(lenbuf), ibuf(2,lenbuf), diag(n), v(n)
      dimension ind(n,*)
      common/io/inp, iout
      write(iout,1) title
      call iosys('create integer "'//itoc(dim)//'d buffers '
     1                    //'" on ham',-1,0,0,' ')
      call rzero(diag,n)
      ntot=0
      if(dim.eq.1) then
         count=0
         do 10 i=1,n
            do 20 j=1,i-1
               count=count+1
               if(count.gt.lenbuf) then
                  ntot=ntot+lenbuf
                  call iosys('write integer "'//itoc(dim)//
     1                       'd buffers" to ham without rewinding',
     2                        2*lenbuf,ibuf,0,' ')
                  call iosys('write integer "'//itoc(dim)//
     1                       'd buffers" to ham without rewinding',
     2                        wptoin(lenbuf),hbuf,0,' ')
                  count=1
               endif                        
               ibuf(1,count)=i
               ibuf(2,count)=j
               hbuf(count)=h01(i,j)
 20         continue   
            diag(i)=diag(i) + h01(i,i) + v(i)
 10      continue   
      elseif(dim.eq.2) then
         count=0
         do 30 i=1,n
            i1=ind(i,1)
            j1=ind(i,2)
            indi=ind(i,3)
            do 40 j=1,i-1
               ham=0.d0
               i2=ind(j,1)
               j2=ind(j,2)
               indj=ind(j,3)
               di1i2=0
               dj1j2=0
               if(i1.eq.i2) then
                  di1i2=1
               endif
               if(j1.eq.j2) then
                  dj1j2=1
               endif
               ham = ham + h01(i1,i2)*dj1j2 + h02(j1,j2)*di1i2
               if(ham.ne.0.d0) then
                  count=count+1
                  if(count.gt.lenbuf) then
                     ntot=ntot+lenbuf
                     call iosys('write integer "'//itoc(dim)//
     1                          'd buffers" to ham without '//
     2                          'rewinding',2*lenbuf,ibuf,0,' ')
                     call iosys('write integer "'//itoc(dim)//
     1                          'd buffers" to ham without '//
     2                          'rewinding',wptoin(lenbuf),hbuf,0,' ')
                     count=1
                  endif
                  ibuf(1,count)=indi
                  ibuf(2,count)=indj
                  hbuf(count)=ham
               endif                        
 40         continue
            diag(indi) = diag(indi) + h01(i1,i1) + h02(j1,j1) + v(indi)
 30      continue
      elseif(dim.eq.3) then
         count=0
         do 50 i=1,n
            i1=ind(i,1)
            j1=ind(i,2)
            k1=ind(i,3)
            indi=ind(i,4)
            do 60 j=1,i-1
               i2=ind(j,1)
               j2=ind(j,2)
               k2=ind(j,3)
               indj=ind(j,4)
               ham=0.d0
               di1i2=0
               dj1j2=0
               dk1k2=0
               if(i1.eq.i2) then
                  di1i2=1
               endif
               if(j1.eq.j2) then
                  dj1j2=1
               endif
               if(k1.eq.k2) then
                  dk1k2=1
               endif
               ham = ham + h01(i1,i2)*dj1j2*dk1k2
     1                   + h02(j1,j2)*di1i2*dk1k2
     2                   + h03(k1,k2)*di1i2*dj1j2
               if(ham.ne.0.d0) then
                  count=count+1
                  if(count.gt.lenbuf) then
                     ntot=ntot+lenbuf
                     call iosys('write integer "'//itoc(dim)//
     1                          'd buffers" to ham without '//
     2                          'rewinding',2*lenbuf,ibuf,0,' ')
                     call iosys('write integer "'//itoc(dim)//
     1                          'd buffers" to ham without '//
     2                          'rewinding',wptoin(lenbuf),hbuf,0,' ')
                     count=1
                  endif                        
                  ibuf(1,count)=indi
                  ibuf(2,count)=indj
                  hbuf(count)=ham             
               endif
 60         continue
            diag(indi) = diag(indi) + h01(i1,i1)
     1                              + h02(j1,j1)
     2                              + h03(k1,k1)
     3                              + v(indi)
 50      continue
      endif
      if(count.gt.0) then
         ntot=ntot+count
         call iosys('write integer "'//itoc(dim)//
     1              'd buffers" to ham without rewinding',2*count,
     2               ibuf,0,' ')
         call iosys('write integer "'//itoc(dim)//
     1              'd buffers" to ham without rewinding',
     2               wptoin(count),hbuf,0,' ')
      endif
      call iosys('endfile "'//itoc(dim)//'d buffers" on ham',0,
     1            0,0,' ')
      call iosys('write integer "'//itoc(dim)//'d number of elements '
     1           //'" to ham',1,ntot,0,' ')
      call iosys('write real "'//itoc(dim)//'d diagonals" to ham',
     1            n,diag,0,' ')
      call iosys('rewind all on lamdat read-and-write',0,0,0,' ')
      write(iout,2) ntot, n
      incore=.false.
      if(ntot.le.lenbuf) then
         incore=.true.
         write(iout,3)
      endif         
      if(prnt) then
         call rdham(hbuf,ibuf,diag,dum,dim,lenbuf,n,title,.false.)
      endif
      return
 1    format(a80)     
 2    format(/,1x,'number non-zero,off-diagonal matrix elements = ',i6,
     1       /,1x,'number of diagonal matrix elements = ',i5)
 3    format(/,1x,'all non-zero matrix elements and indices kept in'
     1            ' core')     
      end       






