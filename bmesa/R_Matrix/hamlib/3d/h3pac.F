*deck h3pac.f
c***begin prologue     h3pac
c***date written       000710   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           one-dim
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            non zero hamiltonian matrix element and indices 
c***                   for three dimensional hamiltonian. 
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       h3pac
      subroutine h3pac(hbuf,buf,diag,v,ind,
     1                 hbufx,bufx,diagx,
     2                 hbufy,bufy,diagy,
     3                 hbufz,bufz,diagz,dim,len,nonz,
     4                 n,incore)
      implicit integer (a-z)
      real*8 hbuf, diag, v, hbufx, diagx, hbufy, diagy, hbufz, diagz
      logical incore
      dimension n(dim+1), len(dim+1), nonz(dim+1)
      dimension hbuf(len(dim+1)), buf(len(dim+1),2), diag(n(dim+1))
      dimension ind(n(3),n(2),n(1)), v(n(dim+1))
      dimension hbufx(len(1)), bufx(len(1),2), diagx(n(1))
      dimension hbufy(len(2)), bufy(len(2),2), diagy(n(2))
      dimension hbufz(len(3)), bufz(len(3),2), diagz(n(3))
      common/io/inp, iout
      if(incore) then
         nonz(dim+1)=0
         do 10 zij=1,nonz(3)
            zi=bufz(zij,1)
            zj=bufz(zij,2)
            do 20 xi=1,n(1)
               do 30 yi=1,n(2)
                  i=ind(zi,yi,xi)
                  j=ind(zj,yi,xi)
                  nonz(dim+1)=nonz(dim+1)+1
                  buf(nonz(dim+1),1)=i
                  buf(nonz(dim+1),2)=j
                  hbuf(nonz(dim+1))=hbufz(zij)
 30            continue
 20         continue
 10      continue   
         do 40 yij=1,nonz(2)
            yi=bufy(yij,1)
            yj=bufy(yij,2)
            do 50 xi=1,n(1)
               do 60 zi=1,n(3)
                  i=ind(zi,yi,xi)
                  j=ind(zi,yj,xi)
                  nonz(dim+1)=nonz(dim+1)+1
                  buf(nonz(dim+1),1)=i
                  buf(nonz(dim+1),2)=j
                  hbuf(nonz(dim+1))=hbufy(yij)  
 60            continue
 50         continue
 40      continue
         do 70 xij=1,nonz(1)
            xi=bufy(xij,1)
            xj=bufy(xij,2)
            do 80 yi=1,n(2)
               do 90 zi=1,n(3)
                  i=ind(zi,yi,xi)
                  j=ind(zi,yi,xj)
                  nonz(dim+1)=nonz(dim+1)+1
                  buf(nonz(dim+1),1)=i
                  buf(nonz(dim+1),2)=j
                  hbuf(nonz(dim+1))=hbufx(xij)  
 90            continue
 80         continue
 70      continue   
      else
         nonz(dim+1)=0
         count=0
         do 100 zij=1,nonz(3)
            zi=bufz(zij,1)
            zj=bufz(zij,2)
            do 110 xi=1,n(1)
               do 120 yi=1,n(2)
                  i=ind(zi,yi,xi)
                  j=ind(zj,yi,xi)
                  count=count+1
                  if(count.gt.len(dim+1)) then
                     nonz(dim+1)=nonz(dim+1)+len(dim+1)
                     call iosys('write integer "hamiltonian buffers" '//
     1                          'to ham without rewinding',
     2                           2*len(dim+1),buf,0,' ')
                     call iosys('write integer "hamiltonian buffers" '//
     1                          'to ham without rewinding',
     2                           wptoin(len(dim+1)),hbuf,0,' ')
                     count=1
                  endif        
                  buf(count,1)=i
                  buf(count,2)=j
                  hbuf(count)=hbufz(zij)  
 120           continue
 110        continue
 100     continue   
         do 130 yij=1,nonz(2)
            yi=bufy(yij,1)
            yj=bufy(yij,2)
            do 140 xi=1,n(1)
               do 150 zi=1,n(3)
                  i=ind(zi,yi,xi)
                  j=ind(zi,yj,xi)
                  count=count+1
                  if(count.gt.len(dim+1)) then
                     nonz(dim+1)=nonz(dim+1)+len(dim+1)
                     call iosys('write integer "hamiltonian buffers" '//
     1                          'to ham without rewinding',
     2                           2*len(dim+1),buf,0,' ')
                     call iosys('write integer "hamiltonian buffers" '//
     1                          'to ham without rewinding',
     2                           wptoin(len(dim+1)),hbuf,0,' ')
                     count=1
                  endif        
                  buf(count,1)=i
                  buf(count,2)=j
                  hbuf(count)=hbufy(yij)  
 150           continue
 140        continue
 130     continue
         do 160 xij=1,nonz(1)
            xi=bufy(xij,1)
            xj=bufy(xij,2)
            do 170 yi=1,n(2)
               do 180 zi=1,n(3)
                  i=ind(zi,yi,xi)
                  j=ind(zi,yi,xj)
                  count=count+1
                  if(count.gt.len(dim+1)) then
                     nonz(dim+1)=nonz(dim+1)+len(dim+1)
                     call iosys('write integer "hamiltonian buffers" '//
     1                          'to ham without rewinding',
     2                           2*len(dim+1),buf,0,' ')
                     call iosys('write integer "hamiltonian buffers" '//
     1                          'to ham without rewinding',
     2                           wptoin(len(dim+1)),hbuf,0,' ')
                     count=1
                  endif        
                  buf(count,1)=i
                  buf(count,2)=j
                  hbuf(count)=hbufx(xij)  
 180           continue
 170        continue
 160     continue   
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
      cnt=0
      do 200 xi=1,n(1)
         do 210 yi=1,n(2)
            do 220 zi=1,n(3)
               cnt=cnt+1
               diag(cnt) = diagx(xi) + diagy(yi) + diagz(zi) + v(cnt)
 220        continue
 210     continue   
 200  continue   
      return
      end       


