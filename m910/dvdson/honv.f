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
      subroutine honv(ibuf,hbuf,diag,veco,vecn,n,nvc,headr,
     1                lenbuf,nel,incore,title,prnt)
      implicit integer (a-z)
      real*8 hbuf, diag
      real*8 veco, vecn
      character*(*) title, headr
      logical incore, prnt
      dimension hbuf(lenbuf), ibuf(2,lenbuf)
      dimension veco(n,nvc), vecn(n,nvc), diag(n)
      common/io/inp, iout
      call rzero(vecn,n*nvc)
      if(incore) then
         do 1000 i=1,nel
            ii=ibuf(1,i)
            jj=ibuf(2,i)
            do 2000 j=1,nvc
               vecn(ii,j) = vecn(ii,j) + hbuf(i)*veco(jj,j)
               vecn(jj,j) = vecn(jj,j) + hbuf(i)*veco(ii,j)
 2000       continue
 1000    continue
      else
         trips=nel/lenbuf
         left=nel-trips*lenbuf
         do 3000 nt=1,trips
            call iosys('read integer '//headr//' from hamiltonian '//
     1                 'without rewinding',2*lenbuf,ibuf,0,' ')
            call iosys('read integer '//headr//' from hamiltonian '//
     1                 'without rewinding',wptoin(lenbuf),hbuf,0,' ')
            do 4000 i=1,lenbuf
               ii=ibuf(1,i)
               jj=ibuf(2,i)
               do 5000 j=1,nvc
                  vecn(ii,j) = vecn(ii,j) + hbuf(i)*veco(jj,j)
                  vecn(jj,j) = vecn(jj,j) + hbuf(i)*veco(ii,j)
 5000          continue
 4000       continue
 3000    continue
         if(left.ne.0) then
            call iosys('read integer '//headr//' from hamiltonian '//
     1                 'without rewinding',2*left,ibuf,0,' ')
            call iosys('read integer '//headr//' from hamiltonian '//
     1                 'without rewinding',wptoin(left),hbuf,0,' ')
            do 6000 i=1,left
               ii=ibuf(1,i)
               jj=ibuf(2,i)
               do 7000 j=1,nvc
                  vecn(ii,j) = vecn(ii,j) + hbuf(i)*veco(jj,j)
                  vecn(jj,j) = vecn(jj,j) + hbuf(i)*veco(ii,j)
 7000          continue
 6000       continue
         endif   
      endif          
      do 8000 i=1,n
         do 9000 j=1,nvc
            vecn(i,j) = vecn(i,j) + diag(i)*veco(i,j)
 9000    continue
 8000 continue
      if(prnt) then
         call prntrm(title,vecn,n,nvc,n,nvc,iout)
      endif             
      return
      end       
