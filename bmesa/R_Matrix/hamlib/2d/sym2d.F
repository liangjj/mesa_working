*deck sym2d.f
c***begin prologue     sym2d
c***date written       000710   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           two-dim
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            non zero hamiltonian elements and indices for 
c***                   two dimensional, symmetrized hamiltonian. 
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       sym2d
      subroutine sym2d(hbuf,buf,diag,v,ind,hi,dim,len,nonz,n,incore)
      implicit integer (a-z)
      real*8 hbuf, diag, v, hi, sqrt2
      logical incore
      dimension n(dim+1)
      dimension hbuf(len), buf(len,2), diag(n(dim+1))
      dimension ind(n(1),n(1)), v(n(dim+1))
      dimension hi(n(1),n(1))
      common/io/inp, iout
      sqrt2=sqrt(2.d0)
c
c 
      nonz=0
      if(incore) then
         do 10 i=1,n(1)
            do 20 j=1,i-1 
c
c              <I| = <ij|   i>j
c
               twoi=ind(i,j)
c          
c              the |J> allowed are either
c                     |ia> , |ja> , |ai> , |aj>
c              the first index must be larger than the second and 
c              |I> must be larger than |j>
c
c              first |ia> and |ai>
c
               do 30 a=1,j-1
                 twoj=ind(i,a)
                 if(hi(j,a).ne.0.d0) then
                    nonz=nonz+1  
                    buf(nonz,1)=twoi
                    buf(nonz,2)=twoj
                    hbuf(nonz)=hi(j,a)  
                 endif
 30            continue
c               
c              now |ja> and |aj>
c          
               do 40 a=j+1,i-1
                  twoj=ind(a,j)
                  if(hi(i,a).ne.0.d0) then
                     nonz=nonz+1
                     buf(nonz,1)=twoi
                     buf(nonz,2)=twoj
                     hbuf(nonz)=hi(i,a)
                  endif
 40            continue   
               do 50 a=1,j-1
                  twoj=ind(j,a)
                  if(hi(i,a).ne.0.d0) then
                     nonz=nonz+1
                     buf(nonz,1)=twoi
                     buf(nonz,2)=twoj
                     hbuf(nonz)=hi(i,a)
                  endif
 50            continue   
 20         continue
 10      continue      
         do 70 i=1,n(1)
            ii=ind(i,i)
            do 80 a=1,n(1)
               if(hi(i,a).ne.0.d0) then
                  if(i.gt.a) then  
                     jj=ind(i,a)
                     twoi=ii
                     twoj=jj
                     nonz=nonz+1
                     buf(nonz,1)=twoi
                     buf(nonz,2)=twoj
                     hbuf(nonz)=sqrt2*hi(i,a)
                  elseif(i.lt.a) then
                     jj=ind(a,i)
                     nonz=nonz+1
                     twoi=jj
                     twoj=ii 
                     buf(nonz,1)=twoi
                     buf(nonz,2)=twoj
                     hbuf(nonz)=sqrt2*hi(i,a)
                  endif 
               endif
 80         continue
 70      continue   
      else
         count=0
         do 100 i=1,n(1)
            do 110 j=1,i-1
               twoi=ind(i,j)
               do 120 a=1,j-1
                  if(hi(j,a).ne.0.d0) then
                     twoj=ind(i,a)
                     count=count+1
                     if(count.gt.len) then
                        nonz=nonz+len
                        call iosys('write integer '//
     1                             '"hamiltonian buffers" to ham '//
     2                             'without rewinding',2*len,buf,0,' ')
                        call iosys('write integer '//
     1                             '"hamiltonian buffers" to ham '//
     2                             'without rewinding',wptoin(len),
     3                              hbuf,0,' ')
                        count=1
                     endif        
                     buf(count,1)=twoi
                     buf(count,2)=twoj
                     hbuf(count)=hi(j,a)  
                  endif
 120           continue
               do 130 a=j+1,i-1
                  if(hi(i,a).ne.0.d0) then
                     twoj=ind(a,j)
                     count=count+1
                     if(count.gt.len) then
                        nonz=nonz+len
                        call iosys('write integer '//
     1                             '"hamiltonian buffers" to ham '//
     2                             'without rewinding',2*len,buf,0,' ')
                        call iosys('write integer '//
     1                             '"hamiltonian buffers" to ham '//
     2                             'without rewinding',wptoin(len),
     3                              hbuf,0,' ')
                        count=1
                     endif        
                     buf(count,1)=twoi
                     buf(count,2)=twoj
                     hbuf(count)=hi(i,a)
                  endif
 130           continue      
               do 150 a=1,j-1
                  if(hi(i,a).ne.0.d0) then
                     twoj=ind(j,a)
                     count=count+1
                     if(count.gt.len) then
                        nonz=nonz+len
                        call iosys('write integer '//
     1                             '"hamiltonian buffers" to ham '//
     2                             'without rewinding',2*len,buf,0,' ')
                        call iosys('write integer '//
     1                             '"hamiltonian buffers" to ham '//
     2                             'without rewinding',wptoin(len),
     3                              hbuf,0,' ')
                        count=1
                     endif        
                     buf(count,1)=twoi
                     buf(count,2)=twoj
                     hbuf(count)=hi(i,a)
                  endif
 150           continue   
 110        continue
 100     continue      
         do 200 i=1,n(1)
            ii=ind(i,i)
            do 210 a=1,n(1)
               if(hi(i,a).ne.0.d0) then
                  if(i.gt.a) then  
                     jj=ind(i,a)
                     twoi=ii
                     twoj=jj
                     count=count+1
                     if(count.gt.len) then
                        nonz=nonz+len
                        call iosys('write integer '//
     1                             '"hamiltonian buffers" to ham '//
     2                             'without rewinding',2*len,buf,0,' ')
                        call iosys('write integer '//
     1                             '"hamiltonian buffers" to ham '//
     2                             'without rewinding',wptoin(len),
     3                              hbuf,0,' ')
                        count=1
                     endif        
                     buf(count,1)=twoi
                     buf(count,2)=twoj
                     hbuf(count)=sqrt2*hi(i,a)
                  elseif(i.lt.a) then
                     jj=ind(a,i)
                     count=count+1
                     twoi=jj
                     twoj=ii 
                     if(count.gt.len) then
                        nonz=nonz+len
                        call iosys('write integer '//
     1                             '"hamiltonian buffers" to ham '//
     2                             'without rewinding',2*len,buf,0,' ')
                        call iosys('write integer '//
     1                             '"hamiltonian buffers" to ham '//
     2                             'without rewinding',wptoin(len),
     3                              hbuf,0,' ')
                        count=1
                     endif        
                     buf(count,1)=twoi
                     buf(count,2)=twoj
                     hbuf(count)=sqrt2*hi(i,a)
                  endif 
               endif
 210        continue
 200     continue   
         if(count.ne.0) then
            nonz=nonz+count
            call iosys('write integer "hamiltonian buffers" to ham '//
     1                 'without rewinding',2*count,buf,0,' ')
            call iosys('write integer "hamiltonian buffers" to ham '//
     1                 'without rewinding',wptoin(count),hbuf,0,' ')
         endif
         call iosys('endfile "hamiltonian buffers" on ham',0,0,0,' ')
         call iosys('write integer "number of non-zero matrix '//
     1              'elements" to ham',1,nonz,0,' ')
         call iosys('rewind all on ham read-and-write',0,0,0,' ')  
      endif
c
c     do the diagonals
c
      do 300 i=1,n(1)
         do 310 j=1,i
            twoi=ind(i,j)
            diag(twoi) = hi(i,i) + hi(j,j) + v(twoi)
 310     continue
 300  continue   
      return
      end       










