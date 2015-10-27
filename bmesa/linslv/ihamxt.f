*deck ihamxt.f
c***begin prologue     ihamxt
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
c***end prologue       ihamxt
      subroutine ihamxt(ham,hbuf,ibuf,ind,diag,lenbuf,n,prnt,ntot,
     1                  incore)
      implicit integer (a-z)
      real*8 ham, hbuf, diag, tmp
      character*1 itoc
      logical prnt, incore
      dimension ham(n,n), hbuf(lenbuf), ibuf(2,lenbuf), diag(n)
      dimension ind(n)
      common/io/inp, iout
      dim=1
      call iosys('create integer "'//itoc(dim)//'d buffers" '//
     1           'on lamdat',-1,0,0,' ')
      call rzero(diag,n)
      ntot=0
      count=0
      do 10 i=1,n
            ii=ind(i)
            do 20 j=1,i-1
               jj=ind(j)
               tmp=ham(i,j)
               if(tmp.ne.0.d0) then
                  count=count+1
                  if(count.gt.lenbuf) then
                     ntot=ntot+lenbuf
                     call iosys('write integer "'//itoc(dim)//
     1                          'd buffers" to lamdat without '//
     2                          'rewinding',2*lenbuf,ibuf,0,' ')
                     call iosys('write integer "'//itoc(dim)//
     1                          'd buffers" to lamdat without '//
     2                          'rewinding',wptoin(lenbuf),hbuf,0,' ')
                     count=1
                  endif
                  ibuf(1,count)=ii
                  ibuf(2,count)=jj
                  hbuf(count)=tmp
               endif                        
 20         continue
            diag(i) = diag(i) + ham(i,i)
 10      continue
      if(count.gt.0) then
         ntot=ntot+count
         call iosys('write integer "'//itoc(dim)//
     1              'd buffers" to lamdat without rewinding',2*count,
     2               ibuf,0,' ')
         call iosys('write integer "'//itoc(dim)//
     1              'd buffers" to lamdat without rewinding',
     2               wptoin(count),hbuf,0,' ')
      endif
      call iosys('endfile "'//itoc(dim)//'d buffers" on lamdat',0,
     1            0,0,' ')
      call iosys('write integer "'//itoc(dim)//'d number of elements '
     1           //'" to lamdat',1,ntot,0,' ')
      call iosys('write real "'//itoc(dim)//'d diagonals '
     1          //'" to lamdat',n,diag,0,' ')
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






