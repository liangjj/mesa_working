*deck conh2.f
c***begin prologue     conh2
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
c***end prologue       conh2
      subroutine conh2(hbuf,buf,diag,v,ind,len,nonz,
     1                  hx,hy,nx,ny,n,incore)
      implicit integer (a-z)
      real*8 hbuf, diag, v, hx, hy
      logical incore
      dimension hbuf(len), buf(len,2), diag(n), ind(ny,nx)
      dimension v(n)
      dimension hx(nx,nx), hy(ny,ny)
      common/io/inp, iout
c
c 
      if(incore) then
         nonz=0
         do 10 xi=1,nx
            do 20 yi=1,ny 
               do 30 yj=1,yi-1
                  if(hy(yi,yj).ne.0.d0) then
                     i=ind(yi,xi)
                     j=ind(yj,xi)
                     nonz=nonz+1
                     buf(nonz,1)=i
                     buf(nonz,2)=j
                     hbuf(nonz)=hy(yi,yj)  
                  endif
 30            continue
 20         continue
 10      continue      
         do 40 yi=1,ny
            do 50 xi=1,nx
               do 60 xj=1,xi-1
                  if(hx(xi,xj).ne.0.d0) then                
                     i=ind(yi,xi)
                     j=ind(yi,xj)
                     nonz=nonz+1 
                     buf(nonz,1)=i
                     buf(nonz,2)=j
                     hbuf(nonz)=hx(xi,xj)  
                  endif
 60            continue   
 50         continue   
 40      continue
      else
         nonz=0
         count=0
         do 100 xi=1,nx
            do 110 yi=1,ny
               do 120 yj=1,yi-1
                  if(hy(yi,yj).ne.0.d0) then
                     i=ind(yi,xi)
                     j=ind(yj,xi)
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
                     buf(count,1)=i
                     buf(count,2)=j
                     hbuf(count)=hy(yi,yj)  
                  endif
 120           continue   
 110        continue
 100     continue      
         do 200 yi=1,ny
            do 210 xi=1,nx
               do 220 xj=1,xi-1
                  if(hx(xi,xj).ne.0.d0) then                
                     i=ind(yi,xi)
                     j=ind(yi,xj)
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
                     buf(count,1)=i
                     buf(count,2)=j
                     hbuf(count)=hx(xi,xj)  
                  endif
 220           continue   
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
            diag(cnt) = hx(xi,xi) + hy(yi,yi) + v(cnt)
 310     continue
 300  continue   
      return
      end       










