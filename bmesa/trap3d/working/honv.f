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
      subroutine honv(hbuf,ibuf,diag,vecold,vecnew,n,nvc,lenbuf,
     1                nel,it,incore)
      implicit integer (a-z)
      real*8 hbuf, diag, vecold, vecnew
      character*1 it
      logical incore
      dimension hbuf(lenbuf), ibuf(2,lenbuf)
      dimension vecold(n,nvc), vecnew(n,nvc), diag(n)
      common/io/inp, iout
      call rzero(vecnew,n*nvc) 
      if(incore) then
         do 10 i=1,nel
            ii=ibuf(1,i)
            jj=ibuf(2,i)
            do 20 j=1,nvc
               vecnew(ii,j) = vecnew(ii,j) + hbuf(i)*vecold(jj,j)
               vecnew(jj,j) = vecnew(jj,j) + hbuf(i)*vecold(ii,j)
   20       continue
   10    continue
      else
         trips=nel/lenbuf
         left=nel-trips*lenbuf
         do 30 nt=1,trips
            call iosys('read integer "'//it//
     1                 'd buffers" from ham without rewinding',
     2                  2*lenbuf,ibuf,0,' ')
            call iosys('read integer "'//it//
     1                 'd buffers" from ham without rewinding',
     2                  wptoin(lenbuf),hbuf,0,' ')
            do 40 i=1,lenbuf
               ii=ibuf(1,i)
               jj=ibuf(2,i)
               do 50 j=1,nvc
                  vecnew(ii,j) = vecnew(ii,j) + hbuf(i)*vecold(jj,j)
                  vecnew(jj,j) = vecnew(jj,j) + hbuf(i)*vecold(ii,j)
 50            continue
 40         continue
 30      continue
         if(left.ne.0) then
            call iosys('read integer "'//it//
     1                 'd buffers" from ham without rewinding',2*left,
     2                    ibuf,0,' ')
            call iosys('read integer "'//it// 
     1                 'd buffers" from ham without rewinding',
     2                  wptoin(left),hbuf,0,' ')
            do 60 i=1,left
               ii=ibuf(1,i)
               jj=ibuf(2,i)
               do 70 j=1,nvc
                  vecnew(ii,j) = vecnew(ii,j) + hbuf(i)*vecold(jj,j)
                  vecnew(jj,j) = vecnew(jj,j) + hbuf(i)*vecold(ii,j)
 70            continue
 60         continue
         endif   
      endif          
      do 80 i=1,n
         do 90 j=1,nvc
            vecnew(i,j) = vecnew(i,j) + diag(i)*vecold(i,j)
 90      continue
 80   continue   
      return
      end       
