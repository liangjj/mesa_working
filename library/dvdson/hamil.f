*deck hamil.f
c***begin prologue     hamil
c***date written       960615   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           hamiltonian
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            calculate and store non-zero matrix elements
c***                   of general real hamiltonian and their indices.
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       hamil
      subroutine hamil(h01,ch01,h02,ch02,h03,ch03,v,hbuf,chbuf,ibuf,
     1                 ind,diag,cdiag,lenbuf,n1,n2,n3,dim,n,prnt,ntot,
     2                 incore,title,mattyp)
c
c
c                      note that hbuf and chbuf are equivalenced in the 
c                      calling program and therefore the actual size of the
c                      array hbuf must be twice lenbuf to accomodate the
c                      matrix elements.
c                            
      implicit integer (a-z)
      real*8 h01, h02, h03
      complex*16 ch01, ch02, ch03, cham, chbuf, cdiag, cdum
      real*8 hbuf, diag, ham, v, dum
      character*(*) title
      character*1 itoc
      logical prnt, incore
      character*(*) mattyp
      dimension h01(n1,n1), h02(n2,n2), h03(n3,n3)
      dimension ch01(n1,n1), ch02(n2,n2), ch03(n3,n3)
      dimension hbuf(lenbuf), ibuf(2,lenbuf), diag(n), cdiag(n), v(n)
      dimension chbuf(lenbuf), ind(n,*)
      common/io/inp, iout
      write(iout,1) title
      call iosys('create integer "'//itoc(dim)//'d buffers '
     1                    //'" on ham',-1,0,0,' ')
      if(mattyp.eq.'complex'.or.mattyp.eq.'real-unsymmetric') then
         call czero(cdiag,n)  
         ntot=0
         if(dim.eq.1) then
            count=0
            do 10 i=1,n
               do 20 j=1,i-1
                  cham=(0.d0,0.d0)
                  cham = cham + ch01(i,j) 
                  if(cham.ne.(0.d0,0.d0)) then
                     count=count+1
                     if(count.gt.lenbuf) then
                        ntot=ntot+lenbuf
                        call iosys('write integer "'//itoc(dim)//
     1                             'd buffers" to ham without '//
     2                             'rewinding',2*lenbuf,ibuf,0,' ')
                        call iosys('write integer "'//itoc(dim)//
     1                             'd buffers" to ham without '//
     2                             'rewinding',2*wptoin(lenbuf),
     3                              hbuf,0,' ')
                        count=1
                     endif
                     ibuf(1,count)=i
                     ibuf(2,count)=j
                     chbuf(count) = cham
                  endif
 20            continue   
               do 30 j=i+1,n
                  cham=(0.d0,0.d0)
                  cham = cham + ch01(i,j)
                  if(cham.ne.(0.d0,0.d0)) then
                     count=count+1
                     if(count.gt.lenbuf) then
                        ntot=ntot+lenbuf
                        call iosys('write integer "'//itoc(dim)//
     1                             'd buffers" to ham without '//
     2                             'rewinding',2*lenbuf,ibuf,0,' ')
                        call iosys('write integer "'//itoc(dim)//
     1                             'd buffers" to ham without '//
     2                             'rewinding',2*wptoin(lenbuf),
     3                              hbuf,0,' ')
                        count=1
                     endif
                     ibuf(1,count)=i
                     ibuf(2,count)=j
                     chbuf(count) = cham
                  endif
 30            continue   
               cdiag(i)=cdiag(i) 
     1                           + ch01(i,i) 
     2                                       + v(i)
 10         continue
         elseif(dim.eq.2) then
            count=0
            do 40 i=1,n
               i1=ind(i,1)
               j1=ind(i,2)
               indi=ind(i,3)
               do 50 j=1,i-1
                  cham=(0.d0,0.d0)
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
                  cham = cham 
     1                        + ch01(i1,i2)*dj1j2 
     2                                            + ch02(j1,j2)*di1i2 
                  if(cham.ne.(0.d0,0.d0)) then
                     count=count+1
                     if(count.gt.lenbuf) then
                        ntot=ntot+lenbuf
                        call iosys('write integer "'//itoc(dim)//
     1                             'd buffers" to ham without '//
     2                             'rewinding',2*lenbuf,ibuf,0,' ')
                        call iosys('write integer "'//itoc(dim)//
     1                             'd buffers" to ham without '//
     2                             'rewinding',2*wptoin(lenbuf),
     3                              hbuf,0,' ')
                        count=1
                     endif
                     ibuf(1,count)=indi
                     ibuf(2,count)=indj
                     chbuf(count) = cham
                  endif                        
 50            continue
               do 60 j=i+1,n
                  cham=(0.d0,0.d0)
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
                  cham = cham 
     1                        + ch01(i1,i2)*dj1j2 
     2                                           + ch02(j1,j2)*di1i2 
                  if(cham.ne.(0.d0,0.d0)) then
                     count=count+1
                     if(count.gt.lenbuf) then
                        ntot=ntot+lenbuf
                        call iosys('write integer "'//itoc(dim)//
     1                             'd buffers" to ham without '//
     2                             'rewinding',2*lenbuf,ibuf,0,' ')
                        call iosys('write integer "'//itoc(dim)//
     1                             'd buffers" to ham without '//
     2                             'rewinding',2*wptoin(lenbuf),
     3                              hbuf,0,' ')
                        count=1
                     endif
                     ibuf(1,count)=indi
                     ibuf(2,count)=indj
                     chbuf(count) = cham
                  endif                        
 60            continue
               cdiag(i) = cdiag(i) 
     1                              + ch01(i1,i1) 
     2                                            + ch02(j1,j1) 
     3                                                 + v(i)
 40         continue
         elseif(dim.eq.3) then
            count=0
            do 70 i=1,n
               i1=ind(i,1)
               j1=ind(i,2)
               k1=ind(i,3)
               indi=ind(i,4)
               do 80 j=1,i-1
                  i2=ind(j,1)
                  j2=ind(j,2)
                  k2=ind(j,3)
                  indj=ind(j,4)
                  cham=(0.d0,0.d0)
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
                  cham = cham + ch01(i1,i2)*dj1j2*dk1k2
     1                   + ch02(j1,j2)*di1i2*dk1k2
     2                   + ch03(k1,k2)*di1i2*dj1j2
                  if(cham.ne.(0.d0,0.d0)) then
                     count=count+1
                     if(count.gt.lenbuf) then
                        ntot=ntot+lenbuf
                        call iosys('write integer "'//itoc(dim)//
     1                             'd buffers" to ham without '//
     2                             'rewinding',2*lenbuf,ibuf,0,' ')
                        call iosys('write integer "'//itoc(dim)//
     1                             'd buffers" to ham without '//
     2                             'rewinding',2*wptoin(lenbuf),
     3                              hbuf,0,' ')
                        count=1
                     endif                        
                     ibuf(1,count)=indi
                     ibuf(2,count)=indj
                     chbuf(count) = cham             
                  endif
 80            continue
               do 90 j=i+1,n
                  i2=ind(j,1)
                  j2=ind(j,2)
                  k2=ind(j,3)
                  indj=ind(j,4)
                  cham=(0.d0,0.d0)
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
                  cham = cham + ch01(i1,i2)*dj1j2*dk1k2
     1                   + ch02(j1,j2)*di1i2*dk1k2
     2                   + ch03(k1,k2)*di1i2*dj1j2
                  if(cham.ne.(0.d0,0.d0)) then
                     count=count+1
                     if(count.gt.lenbuf) then
                        ntot=ntot+lenbuf
                        call iosys('write integer "'//itoc(dim)//
     1                             'd buffers" to ham without '//
     2                             'rewinding',2*lenbuf,ibuf,0,' ')
                        call iosys('write integer "'//itoc(dim)//
     1                             'd buffers" to ham without '//
     2                             'rewinding',2*wptoin(lenbuf),
     3                              hbuf,0,' ')
                        count=1
                     endif                        
                     ibuf(1,count)=indi
                     ibuf(2,count)=indj
                     chbuf(count) = cham             
                  endif
 90            continue
               cdiag(i) = cdiag(i) + ch01(i1,i1)
     1                                   + ch02(j1,j1)
     2                                   + ch03(k1,k1)
     3                                   + v(i)
 70         continue
         endif
         if(count.gt.0) then
            ntot=ntot+count
            call iosys('write integer "'//itoc(dim)//
     1                 'd buffers" to ham without rewinding',2*count,
     2                  ibuf,0,' ')
            call iosys('write integer "'//itoc(dim)//
     1                 'd buffers" to ham without rewinding',
     2                  2*wptoin(count),hbuf,0,' ')
         endif
         call iosys('endfile "'//itoc(dim)//'d buffers" on ham',0,
     1               0,0,' ')
         call iosys('write integer "'//itoc(dim)//'d number of '//
     1              'elements" to ham',1,ntot,0,' ')
         call iosys('write real "'//itoc(dim)//'d diagonals" to ham',
     1               2*n,cdiag,0,' ')
      else
         call rzero(diag,n)  
         ntot=0
         if(dim.eq.1) then
            count=0
            do 1000 i=1,n
               do 2000 j=1,i-1
                  ham=0.d0
                  ham = ham 
     1                      + h01(i,j) 
                  if(ham.ne.0.0d0) then
                     count=count+1
                     if(count.gt.lenbuf) then
                        ntot=ntot+lenbuf
                        call iosys('write integer "'//itoc(dim)//
     1                             'd buffers" to ham without '//
     2                             'rewinding',2*lenbuf,ibuf,0,' ')
                        call iosys('write integer "'//itoc(dim)//
     1                             'd buffers" to ham without '//
     2                             'rewinding',wptoin(lenbuf),
     3                              hbuf,0,' ')
                        count=1
                     endif
                     ibuf(1,count)=i
                     ibuf(2,count)=j
                     hbuf(count) = ham
                  endif
 2000          continue   
               diag(i)=diag(i) 
     1                         + h01(i,i) 
     2                                    + v(i)
 1000       continue
         elseif(dim.eq.2) then
            count=0
            do 3000 i=1,n
               i1=ind(i,1)
               j1=ind(i,2)
               indi=ind(i,3)
               do 4000 j=1,i-1
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
                  ham = ham 
     1                      + h01(i1,i2)*dj1j2 
     2                                         + h02(j1,j2)*di1i2 
                  if(ham.ne.0.d0) then
                     count=count+1
                     if(count.gt.lenbuf) then
                        ntot=ntot+lenbuf
                        call iosys('write integer "'//itoc(dim)//
     1                             'd buffers" to ham without '//
     2                             'rewinding',2*lenbuf,ibuf,0,' ')
                        call iosys('write integer "'//itoc(dim)//
     1                             'd buffers" to ham without '//
     2                             'rewinding',wptoin(lenbuf),
     3                              hbuf,0,' ')
                        count=1
                     endif
                     ibuf(1,count)=indi
                     ibuf(2,count)=indj
                     hbuf(count) = ham
                  endif                        
 4000          continue
               diag(i) = diag(i) 
     1                            + h01(i1,i1) 
     2                                         + h02(j1,j1) 
     3                                                  + v(i)
 3000       continue
         elseif(dim.eq.3) then
            count=0
            do 5000 i=1,n
               i1=ind(i,1)
               j1=ind(i,2)
               k1=ind(i,3)
               indi=ind(i,4)
               do 6000 j=1,i-1
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
     1                            + h02(j1,j2)*di1i2*dk1k2
     2                                  + h03(k1,k2)*di1i2*dj1j2
                  if(ham.ne.0.d0) then
                     count=count+1
                     if(count.gt.lenbuf) then
                        ntot=ntot+lenbuf
                        call iosys('write integer "'//itoc(dim)//
     1                             'd buffers" to ham without '//
     2                             'rewinding',2*lenbuf,ibuf,0,' ')
                        call iosys('write integer "'//itoc(dim)//
     1                             'd buffers" to ham without '//
     2                             'rewinding',wptoin(lenbuf),
     3                              hbuf,0,' ')
                        count=1
                     endif                        
                     ibuf(1,count)=indi
                     ibuf(2,count)=indj
                     hbuf(count) = ham             
                  endif
 6000          continue
               diag(i) = diag(i) + h01(i1,i1)
     1                                        + h02(j1,j1)
     2                                                    + h03(k1,k1)
     3                                        + v(i)
 5000       continue
         endif
         if(count.gt.0) then
            ntot=ntot+count
            call iosys('write integer "'//itoc(dim)//
     1                 'd buffers" to ham without rewinding',2*count,
     2                  ibuf,0,' ')
            call iosys('write integer "'//itoc(dim)//
     1                 'd buffers" to ham without rewinding',
     2                  wptoin(count),hbuf,0,' ')
         endif
         call iosys('endfile "'//itoc(dim)//'d buffers" on ham',0,
     1               0,0,' ')
         call iosys('write integer "'//itoc(dim)//'d number of '//
     1              'elements" to ham',1,ntot,0,' ')
         call iosys('write real "'//itoc(dim)//'d diagonals" to ham',
     1               n,diag,0,' ')
      endif         
      call iosys('rewind all on ham read-and-write',0,0,0,' ')  
      write(iout,2) ntot, n
      incore=.false.
      if(ntot.le.lenbuf) then
         incore=.true.
         write(iout,3)
      endif
      if(prnt) then
         if(mattyp.eq.'complex'.or.mattyp.eq.'real-unsymmetric') then
            call cdham(chbuf,ibuf,cdiag,cdum,dim,lenbuf,n,title,.false.)
         else
            call rdham(hbuf,ibuf,diag,dum,dim,lenbuf,n,title,.false.)
         endif
      endif
      return
 1    format(a80)     
 2    format(/,1x,'number non-zero,off-diagonal matrix elements = ',i6,
     1       /,1x,'number of diagonal matrix elements = ',i5)
 3    format(/,1x,'all non-zero matrix elements and indices kept in'
     1            ' core')     
      end       






