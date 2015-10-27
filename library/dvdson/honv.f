*deck honv.f
c***begin prologue     honv
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***purpose            multiply the hamiltonian matrix on a vector.
c***references         
c
c***routines called    
c***end prologue       honv
      subroutine honv(hbuf,chbuf,ibuf,diag,cdiag,veco,cveco,vecn,cvecn,
     1                n,nvc,lenbuf,nel,it,incore,mattyp,
     2                title,iter,prnt)
      implicit integer (a-z)
      real*8 hbuf, diag
      complex*16 chbuf, cdiag
      real*8 veco, vecn
      complex*16 cveco, cvecn
      character*1 it
      character*(*) title
      logical incore, prnt
      character*(*) mattyp
      dimension hbuf(lenbuf), ibuf(2,lenbuf)
      dimension chbuf(lenbuf)
      dimension veco(n,nvc), vecn(n,nvc), diag(n)
      dimension cveco(n,nvc), cvecn(n,nvc), cdiag(n)
      common/io/inp, iout
      if(mattyp.eq.'complex'.or.mattyp.eq.'real-unsymmetric') then
         call czero(cvecn,n*nvc)
         if(incore) then
            do 10 i=1,nel
               ii=ibuf(1,i)
               jj=ibuf(2,i)
               do 20 j=1,nvc
                  cvecn(ii,j) = cvecn(ii,j) + chbuf(i)*cveco(jj,j)
 20            continue
 10         continue
         else
            trips=nel/lenbuf
            left=nel-trips*lenbuf
            do 30 nt=1,trips
               call iosys('read integer "'//it//
     1                    'd buffers" from ham without rewinding',
     2                     2*lenbuf,ibuf,0,' ')
               call iosys('read integer "'//it//
     1                    'd buffers" from ham without rewinding',
     2                     2*wptoin(lenbuf),hbuf,0,' ')
               do 40 i=1,lenbuf
                  ii=ibuf(1,i)
                  jj=ibuf(2,i)
                  do 50 j=1,nvc
                     cvecn(ii,j) = cvecn(ii,j) + chbuf(i)*cveco(jj,j)
 50               continue
 40            continue
 30         continue
            if(left.ne.0) then
               call iosys('read integer "'//it//
     1                    'd buffers" from ham without rewinding',
     2                     2*left,ibuf,0,' ')
               call iosys('read integer "'//it// 
     1                    'd buffers" from ham without rewinding',
     2                     2*wptoin(left),hbuf,0,' ')
               do 60 i=1,left
                  ii=ibuf(1,i)
                  jj=ibuf(2,i)
                  do 70 j=1,nvc
                     cvecn(ii,j) = cvecn(ii,j) + chbuf(i)*cveco(jj,j)
 70               continue
 60            continue
            endif   
         endif          
         do 80 i=1,n
            do 90 j=1,nvc
               cvecn(i,j) = cvecn(i,j) + cdiag(i)*cveco(i,j)
 90         continue
 80      continue
         if(prnt) then
            call prntcm(title,cvecn,n,nvc,n,nvc,iout)
         endif             
      else
         call rzero(vecn,n*nvc)
         if(incore) then
            do 1000 i=1,nel
               ii=ibuf(1,i)
               jj=ibuf(2,i)
               do 2000 j=1,nvc
                  vecn(ii,j) = vecn(ii,j) + hbuf(i)*veco(jj,j)
                  vecn(jj,j) = vecn(jj,j) + hbuf(i)*veco(ii,j)
 2000          continue
 1000       continue
         else
            trips=nel/lenbuf
            left=nel-trips*lenbuf
            do 3000 nt=1,trips
               call iosys('read integer "'//it//
     1                    'd buffers" from ham without rewinding',
     2                     2*lenbuf,ibuf,0,' ')
               call iosys('read integer "'//it//
     1                    'd buffers" from ham without rewinding',
     2                     wptoin(lenbuf),hbuf,0,' ')
               do 4000 i=1,lenbuf
                  ii=ibuf(1,i)
                  jj=ibuf(2,i)
                  do 5000 j=1,nvc
                     vecn(ii,j) = vecn(ii,j) + hbuf(i)*veco(jj,j)
 5000             continue
 4000          continue
 3000       continue
            if(left.ne.0) then
               call iosys('read integer "'//it//
     1                    'd buffers" from ham without rewinding',
     2                     2*left,ibuf,0,' ')
               call iosys('read integer "'//it// 
     1                    'd buffers" from ham without rewinding',
     2                     wptoin(left),hbuf,0,' ')
               do 6000 i=1,left
                  ii=ibuf(1,i)
                  jj=ibuf(2,i)
                  do 7000 j=1,nvc
                     vecn(ii,j) = vecn(ii,j) + hbuf(i)*veco(jj,j)
 7000             continue
 6000          continue
            endif   
         endif          
         do 8000 i=1,n
            do 9000 j=1,nvc
               vecn(i,j) = vecn(i,j) + diag(i)*veco(i,j)
 9000       continue
 8000    continue
         if(prnt) then
            call prntrm(title,vecn,n,nvc,n,nvc,iout)
         endif             
      endif     
      return
      end       
