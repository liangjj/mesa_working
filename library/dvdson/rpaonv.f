*deck rpaonv.f
c***begin prologue     rpaonv
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***purpose            multiply the rpa hamiltonian matrix on a vector.
c***references         
c
c***routines called    
c***end prologue       rpaonv
      subroutine rpaonv(hbufa,hbufb,ibufa,ibufb,diaga,veco,vecn,
     1                  t,n,nvc,lenbuf,neapb,neamb,incore,
     2                  title,iter,prnt)
      implicit integer (a-z)
      real*8 hbufa, hbufb, diaga
      real*8 veco, vecn, t
      character*(*) title
      logical incore, prnt
      dimension hbufa(lenbuf), hbufb(lenbuf), t(n,nvc)
      dimension ibufa(2,lenbuf), ibufb(2,lenbuf)
      dimension veco(n,nvc), vecn(n,nvc), diaga(n)
      common/io/inp, iout
      call rzero(t,n*nvc)
      call rzero(vecn,n*nvc)
      if(incore) then
         do 10 i=1,neapb
            ii=ibufa(1,i)
            jj=ibufa(2,i)
            do 20 j=1,nvc
               t(ii,j) = t(ii,j) + hbufa(i)*veco(jj,j)
 20         continue
 10      continue
         do 30 i=1,n
            do 40 j=1,nvc
               t(i,j) = t(i,j) + diaga(i)*veco(i,j)
 40         continue
 30      continue
         do 50 i=1,neamb
            ii=ibufb(1,i)
            jj=ibufb(2,i)
            do 60 j=1,nvc
               vecn(ii,j) = vecn(ii,j) + hbufb(i)*t(jj,j)
 60         continue
 50      continue                  
      else
         trips=neapb/lenbuf
         left=neapb-trips*lenbuf
         do 100 nt=1,trips
            call iosys('read integer "a plus b matrix '//
     1                 'buffers" from ham without rewinding',
     2                  2*lenbuf,ibufa,0,' ')
            call iosys('read integer "a plus b matrix '//
     1                 'buffers from ham without rewinding',
     2                  wptoin(lenbuf),hbufa,0,' ')            
            do 200 i=1,lenbuf
               ii=ibufa(1,i)
               jj=ibufa(2,i)
               do 300 j=1,nvc
                  t(ii,j) = t(ii,j) + hbufa(i)*veco(jj,j)
 300           continue
 200        continue
 100     continue
         if(left.ne.0) then
            call iosys('read integer "a plus b matrix '//
     1                 'buffers" from ham without rewinding',
     2                  2*left,ibufa,0,' ')
            call iosys('read integer "a plus b matrix '//
     1                 'buffers from ham without rewinding',
     2                  wptoin(left),hbufa,0,' ')   
            do 400 i=1,left
               ii=ibufa(1,i)
               jj=ibufa(2,i)
               do 500 j=1,nvc
                  t(ii,j) = t(ii,j) + hbufa(i)*veco(jj,j)
 500           continue
 400        continue
         endif
         do 600 i=1,n
            do 700 j=1,nvc
               t(i,j) = t(i,j) + diaga(i)*veco(i,j)
 700        continue
 600     continue
         trips=neamb/lenbuf
         left=neamb-trips*lenbuf
         do 1000 nt=1,trips
            call iosys('read integer "a minus b matrix '//
     1                 'buffers" from ham without rewinding',
     2                  2*lenbuf,ibufb,0,' ')
            call iosys('read integer "a minus b matrix '//
     1                 'buffers from ham without rewinding',
     2                  wptoin(lenbuf),hbufb,0,' ')            
            do 2000 i=1,lenbuf
               ii=ibufb(1,i)
               jj=ibufb(2,i)
               do 3000 j=1,nvc
                  vecn(ii,j) = vecn(ii,j) + hbufb(i)*t(jj,j)
 3000          continue
 2000       continue
 1000    continue
         if(left.ne.0) then
            call iosys('read integer "a minus b matrix '//
     1                 'buffers" from ham without rewinding',
     2                  2*left,ibufb,0,' ')
            call iosys('read integer "a minus b matrix '//
     1                 'buffers from ham without rewinding',
     2                  wptoin(left),hbufb,0,' ')   
            do 4000 i=1,left
               ii=ibufb(1,i)
               jj=ibufb(2,i)
               do 5000 j=1,nvc
                  vecn(ii,j) = vecn(ii,j) + hbufb(i)*t(jj,j)
 5000          continue
 4000       continue
         endif 
      endif          
      if(prnt) then
         call prntrm(title,vecn,n,nvc,n,nvc,iout)
      endif             
      return
      end       
