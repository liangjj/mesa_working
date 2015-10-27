*deck ihamxt.f
c***begin prologue     ihamxt
c***date written       960615   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           hamiltonian
c***author             schneider, barry (nsf)
c***source             
c***purpose            calculate and store non-zero matrix elements
c***                   of hamiltonian and their indices in dvr
c***                   representation.
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       ihamxt
      subroutine ihamxt(ham,hbuf,ibuf,ind,diag,lenbuf,n,prnt,
     1                  ntot,incore)
      implicit integer (a-z)
      complex*16 ham, hbuf, diag
      character*1 itoc
      logical prnt, incore
      dimension ham(n,n), hbuf(lenbuf), ibuf(2,lenbuf), diag(n)
      dimension ind(n)
      common/io/inp, iout
      dim=1
      call iosys('create integer "'//itoc(dim)//'d buffers" '//
     1           'on lamdat',-1,0,0,' ')
      call czero(diag,n)
      ntot=0
      count=0
      do 10 i=1,n
         ii=ind(i)
         do 20 j=1,n
            if(i.ne.j) then         
               jj=ind(j)
               if(ham(i,j).ne.(0.d0,0.d0)) then
                  count=count+1
                  if(count.gt.lenbuf) then
                     ntot=ntot+lenbuf
                     call iosys('write integer "'//itoc(dim)//
     1                          'd buffers" to lamdat without '//
     2                          'rewinding',2*lenbuf,ibuf,0,' ')
                     call iosys('write integer "'//itoc(dim)//
     1                          'd buffers" to lamdat without '//
     2                          'rewinding',wptoin(2*lenbuf),hbuf,
     3                           0,' ')
                     count=1
                  endif
                  ibuf(1,count)=ii
                  ibuf(2,count)=jj
                  hbuf(count)=ham(i,j)
               endif                        
            else
               diag(ii) = ham(i,i)
            endif 
 20      continue
 10   continue
      if(count.gt.0) then
         ntot=ntot+count
         call iosys('write integer "'//itoc(dim)//
     1              'd buffers" to lamdat without rewinding',2*count,
     2               ibuf,0,' ')
         call iosys('write integer "'//itoc(dim)//
     1              'd buffers" to lamdat without rewinding',
     2               wptoin(2*count),hbuf,0,' ')
      endif
      call iosys('endfile "'//itoc(dim)//'d buffers" on lamdat',0,
     1            0,0,' ')
      call iosys('write integer "'//itoc(dim)//'d number of elements '
     1           //'" to lamdat',1,ntot,0,' ')
      call iosys('write real "'//itoc(dim)//'d diagonals '
     1          //'" to lamdat',2*n,diag,0,' ')
      call iosys('rewind all on lamdat read-and-write',0,0,0,' ')
      write(iout,2) ntot, n
      incore=.false.
      if(ntot.le.lenbuf) then
         incore=.true.
         write(iout,3)
      endif         
      if(prnt) then
         call rdham(hbuf,ibuf,diag,dim,lenbuf,n)
      endif
      return
 1    format(a80)     
 2    format(/,1x,'number non-zero,off-diagonal matrix elements = ',i6,
     1       /,1x,'number of diagonal matrix elements = ',i5)
 3    format(/,1x,'all non-zero matrix elements and indices kept in'
     1            ' core')     
      end       






