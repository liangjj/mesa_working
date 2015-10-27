*deck hamfl2.f
c***begin prologue     hamfl2
c***date written       000710   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           two-dim
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            non zero hamiltonian elements and indices for 
c***                   two dimensional packed hamiltonian. 
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       hamfl2
      subroutine hamfl2(hbuf,buf,diag,v,ind,len,nonz,
     1                  hbufx,bufx,diagx,lenx,nonzx,
     2                  hbufy,bufy,diagy,leny,nonzy,
     3                  nx,ny,n,incore)
      implicit integer (a-z)
      real*8 hbuf, diag, v, hbufx, diagx, hbufy, diagy
      logical incore
      dimension hbuf(len), buf(len,2), diag(n), ind(ny,nx)
      dimension v(n)
      dimension hbufx(lenx), bufx(lenx,2), diagx(nx)
      dimension hbufy(leny), bufy(leny,2), diagy(ny)
      common/io/inp, iout
c
c 
      if(incore) then
         nonz=0
         do 10 xi=1,nx
            do 20 yij=1,nonzy 
               yi=bufy(yij,1)
               yj=bufy(yij,2)
               i=ind(yi,xi)
               j=ind(yj,xi)
               nonz=nonz+1
               buf(nonz,1)=i
               buf(nonz,2)=j
               hbuf(nonz)=hbufy(yij)  
 20         continue
 10      continue      
         do 30 yi=1,ny
            do 40 xij=1,nonzx
               xi=bufx(xij,1) 
               xj=bufx(xij,2) 
               i=ind(yi,xi)
               j=ind(yi,xj)
               nonz=nonz+1 
               buf(nonz,1)=i
               buf(nonz,2)=j
               hbuf(nonz)=hbufx(xij)  
 40         continue   
 30      continue
      else
         nonz=0
         count=0
         do 100 xi=1,nx
            do 110 yij=1,nonzy
               yi=bufy(yij,1)
               yj=bufy(yij,2)
               i=ind(yi,xi)
               j=ind(yj,xi)
               count=count+1
               if(count.gt.len) then
                  nonz=nonz+len
                  call iosys('write integer "hamiltonian buffers" '//
     1                       'to ham without rewinding',
     2                        2*len,buf,0,' ')
                  call iosys('write integer "hamiltonian buffers" '//
     1                       'to ham without rewinding',
     2                        wptoin(len),hbuf,0,' ')
                  count=1
               endif        
               buf(count,1)=i
               buf(count,2)=j
               hbuf(count)=hbufy(yij)  
 110        continue
 100     continue      
         do 200 yi=1,ny
            do 210 xij=1,nonzx
               xi=bufx(xij,1) 
               xj=bufx(xij,2) 
               i=ind(yi,xi)
               j=ind(yi,xj)
               count=count+1
               if(count.gt.len) then
                  nonz=nonz+len
                  call iosys('write integer "hamiltonian buffers" '//
     1                       'to ham without rewinding',
     2                        2*len,buf,0,' ')
                  call iosys('write integer "hamiltonian buffers" '//
     1                       'to ham without rewinding',
     2                        wptoin(len),hbuf,0,' ')
                  count=1
               endif        
               buf(count,1)=i
               buf(count,2)=j
               hbuf(count)=hbufx(xij)  
 210        continue   
 200     continue
         if(count.gt.0) then
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
      cnt=0
      do 300 xi=1,nx
         do 310 yi=1,ny
            cnt=cnt+1
            diag(cnt) = diagx(xi) + diagy(yi) + v(cnt)
 310     continue
 300  continue   
      return
      end       










