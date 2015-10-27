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
      subroutine hamil(eig1,eig2,eig3,t1,t2,t3,q3,hbuf,ibuf,ind,
     1                 diag,lenbuf,nply,dim,matdim,pottyp,prnt,
     2                 title,ntot,incore)
      implicit integer (a-z)
      real*8 t1, t2, t3, eig1, eig2, eig3, q3
      real*8 hbuf, diag, ham, tmp1, tmp2
      character*(*) title
      character*(*) pottyp
      character*1 itoc
      logical prnt, incore
      dimension eig1(nply), eig2(nply), eig3(nply)
      dimension t1(nply,nply), t2(nply,nply), t3(nply,nply)
      dimension hbuf(lenbuf), ibuf(2,lenbuf), diag(matdim)
      dimension ind(matdim,*)
      common/io/inp, iout
      write(iout,1) title
      call iosys('create integer "'//itoc(dim)//'d buffers '
     1                    //title//'" on ham',-1,0,0,' ')
      call rzero(diag,matdim)
      ntot=0
      if(dim.eq.1) then
         count=0
         do 10 i=1,nply
            do 20 j=1,i-1
               count=count+1
               if(count.gt.lenbuf) then
                  ntot=ntot+lenbuf
                  call iosys('write integer "'//itoc(dim)//'d buffers '
     1                       //title//'" to ham without rewinding',
     2                       2*lenbuf,ibuf,0,' ')
                  call iosys('write integer "'//itoc(dim)//'d buffers '
     1                       //title//'" to ham without rewinding',
     2                       wptoin(lenbuf),hbuf,0,' ')
                  count=1
               endif                        
               ibuf(1,count)=i
               ibuf(2,count)=j
               hbuf(count)=t1(i,j)
 20         continue   
            diag(i)=diag(i) + t1(i,i)
 10      continue   
c     now add in the potential energy contribution.
         if(pottyp.eq.'exponential') then
            do 30 i=1,nply
               diag(i) = diag(i) + exp(-eig1(i))          
 30         continue
         elseif(pottyp.eq.'one') then   
            do 40 i=1,nply
               diag(i) = diag(i) - 1.d0
 40        continue
         elseif(pottyp.eq.'coulomb') then
            do 50 i=1,matdim
               diag(i) = diag(i) - 1.d0/sqrt(eig1(i))
 50         continue
         endif
      elseif(dim.eq.2) then
         count=0
         do 60 i=1,matdim
            i1=ind(i,1)
            j1=ind(i,2)
            indi=ind(i,3)
            do 70 j=1,i-1
               ham=0.d0
               i2=ind(j,1)
               j2=ind(j,2)
               indj=ind(j,3)
               if(j1.eq.j2) then
                  ham = ham + t1(i1,i2)
                  if(i1.eq.i2) then
                     ham = ham + t2(j1,j1)
                  endif
               else
                  if(i1.eq.i2) then
                     ham = ham + t2(j1,j2)
                  endif
               endif
               if(ham.ne.0.d0) then
                  count=count+1
                  if(count.gt.lenbuf) then
                     ntot=ntot+lenbuf
                     call iosys('write integer "'//itoc(dim)//
     1                          'd buffers '//title//'" to ham '//
     2                          'without rewinding',2*lenbuf,
     3                           ibuf,0,' ')
                     call iosys('write integer "'//itoc(dim)//
     1                          'd buffers '//title//'" to ham '//
     2                          'without rewinding',wptoin(lenbuf),
     3                           hbuf,0,' ')
                     count=1
                  endif
                  ibuf(1,count)=indi
                  ibuf(2,count)=indj
                  hbuf(count)=ham
               endif                        
 70         continue
            diag(indi) = diag(indi) + t1(i1,i1) + t2(j1,j1)
 60      continue
c     now add in the potential energy contribution.
c     the third coordinate is held constant at q3 on this particular face.
         if(pottyp.eq.'exponential') then
            tmp1=q3*q3
            do 100 i=1,matdim
               i1=ind(i,1)
               j1=ind(i,2)
               indi=ind(i,3)
               tmp2=sqrt( eig1(i1)*eig1(i1) + eig2(j1)*eig2(j1) +tmp1 )
               diag(indi) = diag(indi) + exp(-tmp2)          
 100        continue
         elseif(pottyp.eq.'one') then   
            do 110 i=1,matdim
               i1=ind(i,1)
               j1=ind(i,2)
               indi=ind(i,3)
               diag(indi) = diag(indi) - 1.d0
 110        continue
         elseif(pottyp.eq.'coulomb') then
            tmp1=q3*q3
            do 120 i=1,matdim
               i1=ind(i,1)
               j1=ind(i,2)
               indi=ind(i,3)
               tmp2=eig1(i1)*eig1(i1) + eig2(j1)*eig2(j1) + tmp1
               diag(indi) = diag(indi) - 1.d0/sqrt(tmp2)
 120        continue
         endif
      elseif(dim.eq.3) then
         count=0
         do 200 i=1,matdim
            i1=ind(i,1)
            j1=ind(i,2)
            k1=ind(i,3)
            indi=ind(i,4)
            do 210 j=1,i-1
               i2=ind(j,1)
               j2=ind(j,2)
               k2=ind(j,3)
               indj=ind(j,4)
               ham=0.d0
               if(k1.eq.k2) then
                  if(j1.eq.j2) then
                     ham = ham + t1(i1,i2)
                     if(i1.eq.i2) then
                        ham = ham + t2(j1,j1) + t3(k1,k1)
                     endif
                  else
                     if(i1.eq.i2) then
                        ham = ham + t2(j1,j2)
                     endif
                  endif
               else
                  if(i1.eq.i2.and.j1.eq.j2) then
                     ham = ham + t3(k1,k2)
                  endif
               endif
               if(ham.ne.0.d0) then
                  count=count+1
                  if(count.gt.lenbuf) then
                     ntot=ntot+lenbuf
                     call iosys('write integer "'//itoc(dim)//
     1                          'd buffers '//title//'" to ham '//
     2                          'without rewinding',2*lenbuf,
     3                           ibuf,0,' ')
                     call iosys('write integer "'//itoc(dim)//
     1                          'd buffers '//title//'" to ham '//
     2                          'without rewinding',wptoin(lenbuf),
     3                           hbuf,0,' ')
                     count=1
                  endif                        
                  ibuf(1,count)=indi
                  ibuf(2,count)=indj
                  hbuf(count)=ham             
               endif
 210        continue
            diag(indi) = diag(indi) + t1(i1,i1) + t2(j1,j1) + t3(k1,k1)
 200     continue
c     now add in the potential energy contribution.
         if(pottyp.eq.'exponential') then
            do 300 i=1,matdim
               i1=ind(i,1)
               j1=ind(i,2)
               k1=ind(i,3)
               indi=ind(i,4)
               tmp1=sqrt ( eig1(i1)*eig1(i1) + eig2(j1)*eig2(j1) 
     1                                       + eig3(k1)*eig3(k1) )
               diag(indi) = diag(indi) + exp(-tmp1)
 300        continue
         elseif(pottyp.eq.'one') then   
            do 400 i=1,matdim
               i1=ind(i,1)
               j1=ind(i,2)
               k1=ind(i,3)
               indi=ind(i,4)
               diag(indi) = diag(indi) - 1.d0
 400        continue
         elseif(pottyp.eq.'coulomb') then
            do 500 i=1,matdim
               i1=ind(i,1)
               j1=ind(i,2)
               k1=ind(i,3)
               indi=ind(i,4)
               tmp1=eig1(i1)*eig1(i1) + eig2(j1)*eig2(j1)
     1                                + eig3(j1)*eig3(j1)
               diag(indi) = diag(indi) - 1.d0/sqrt(tmp1)
 500        continue
         endif
      endif
      if(count.gt.0) then
         ntot=ntot+count
         call iosys('write integer "'//itoc(dim)//
     1              'd buffers '//title//'" to ham '//
     2              'without rewinding',2*count,
     3               ibuf,0,' ')
         call iosys('write integer "'//itoc(dim)//
     1              'd buffers '//title//'" to ham '//
     2              'without rewinding',wptoin(count),
     3               hbuf,0,' ')
      endif
      call iosys('endfile "'//itoc(dim)//'d buffers '//title//
     1           '" on ham',0,0,0,' ')
      call iosys('write integer "'//itoc(dim)//'d number of elements '
     1           //title//'" to ham',1,ntot,0,' ')
      call iosys('write real "'//itoc(dim)//'d diagonals '
     1          //title//'" to ham',matdim,diag,0,' ')
      call iosys('rewind all on ham read-and-write',0,0,0,' ')
      write(iout,2) ntot, matdim
      incore=.false.
      if(ntot.le.lenbuf) then
         incore=.true.
         write(iout,3)
      endif         
      if(prnt) then
         call rdham(hbuf,ibuf,diag,dim,lenbuf,matdim,title)
      endif
      return
 1    format(a80)     
 2    format(/,1x,'number non-zero,off-diagonal matrix elements = ',i6,
     1       /,1x,'number of diagonal matrix elements = ',i5)
 3    format(/,1x,'all non-zero matrix elements and indices kept in'
     1            ' core')     
      end       






