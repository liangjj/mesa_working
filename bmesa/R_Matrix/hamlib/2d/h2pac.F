*deck h2pac.f
c***begin prologue     h2pac
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
c***end prologue       h2pac
      subroutine h2pac(hbuf,buf,diag,v,ind,hbufx,bufx,diagx,
     2                 hbufy,bufy,diagy,dim,len,nonz,n,incore)
      implicit integer (a-z)
      real*8 hbuf, diag, v, hbufx, diagx, hbufy, diagy
      logical incore
      dimension n(dim+1), len(dim+1), nonz(dim+1)
      dimension hbuf(len(dim+1)), buf(len(dim+1),2), diag(n(dim+1))
      dimension ind(n(2),n(1)), v(n(dim+1))
      dimension hbufx(len(1)), bufx(len(1),2), diagx(n(1))
      dimension hbufy(len(2)), bufy(len(2),2), diagy(n(2))
      common/io/inp, iout
c
c
      if(incore) then
         nonz(dim+1)=0
         do 10 xi=1,n(1)
            do 20 yij=1,nonz(2) 
               yi=bufy(yij,1)
               yj=bufy(yij,2)
               i=ind(yi,xi)
               j=ind(yj,xi)
               nonz(dim+1)=nonz(dim+1)+1
               buf(nonz(dim+1),1)=i
               buf(nonz(dim+1),2)=j
               hbuf(nonz(dim+1))=hbufy(yij)  
 20         continue
 10      continue      
         do 30 yi=1,n(2)
            do 40 xij=1,nonz(1)
               xi=bufx(xij,1) 
               xj=bufx(xij,2) 
               i=ind(yi,xi)
               j=ind(yi,xj)
               nonz(dim+1)=nonz(dim+1)+1 
               buf(nonz(dim+1),1)=i
               buf(nonz(dim+1),2)=j
               hbuf(nonz(dim+1))=hbufx(xij)  
 40         continue   
 30      continue
      else
         nonz(dim+1)=0
         count=0
         do 100 xi=1,n(1)
            do 110 yij=1,nonz(2)
               yi=bufy(yij,1)
               yj=bufy(yij,2)
               i=ind(yi,xi)
               j=ind(yj,xi)
               count=count+1
               if(count.gt.len(dim+1)) then
                  nonz(dim+1)=nonz(dim+1)+len(dim+1)
                  call iosys('write integer "hamiltonian buffers" '//
     1                       'to ham without rewinding',
     2                        2*len(dim+1),buf,0,' ')
                  call iosys('write integer "hamiltonian buffers" '//
     1                       'to ham without rewinding',
     2                        wptoin(len(dim+1)),hbuf,0,' ')
                  count=1
               endif        
               buf(count,1)=i
               buf(count,2)=j
               hbuf(count)=hbufy(yij)  
 110        continue
 100     continue      
         do 200 yi=1,n(2)
            do 210 xij=1,nonz(1)
               xi=bufx(xij,1) 
               xj=bufx(xij,2) 
               i=ind(yi,xi)
               j=ind(yi,xj)
               count=count+1
               if(count.gt.len(dim+1)) then
                  nonz(dim+1)=nonz(dim+1)+len(dim+1)
                  call iosys('write integer "hamiltonian buffers" '//
     1                       'to ham without rewinding',
     2                        2*len(dim+1),buf,0,' ')
                  call iosys('write integer "hamiltonian buffers" '//
     1                       'to ham without rewinding',
     2                        wptoin(len(dim+1)),hbuf,0,' ')
                  count=1
               endif        
               buf(count,1)=i
               buf(count,2)=j
               hbuf(count)=hbufx(xij)  
 210        continue   
 200     continue
         if(count.gt.0) then
            nonz(dim+1)=nonz(dim+1)+count
            call iosys('write integer "hamiltonian buffers" to ham '//
     1                 'without rewinding',2*count,buf,0,' ')
            call iosys('write integer "hamiltonian buffers" to ham '//
     1                 'without rewinding',wptoin(count),hbuf,0,' ')
         endif
         call iosys('endfile "hamiltonian buffers" on ham',0,0,0,' ')
         call iosys('write integer "number of non-zero matrix '//
     1              'elements" to ham',1,nonz(dim+1),0,' ')
         call iosys('rewind all on ham read-and-write',0,0,0,' ')  
      endif
c
c     do the diagonals
c
      cnt=0
      do 300 xi=1,n(1)
         do 310 yi=1,n(2)
            cnt=cnt+1
            diag(cnt) = diagx(xi) + diagy(yi) + v(cnt)
 310     continue
 300  continue   
      return
      end       










