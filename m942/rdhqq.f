*deck rdhqq.f
c***begin prologue     rdhqq
c***date written       960615   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           rdhqqiltonian
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            non-zero elements of hqq
c***                   
c***description                           
c***                   
c***references         
c
c***routines called    
c***end prologue       rdhqq
      subroutine rdhqq(ibuf,rbuf,hqqf,hqqt,diagqq,header,n,ntri,lenbuf,
     #                 nonz,nodsk,type,prnt)
      implicit integer (a-z)
      real*8 rbuf, hqqf, hqqt, diagqq
      logical nodsk, prnt
      character*(*) type, header
      dimension rbuf(lenbuf), ibuf(2,lenbuf)
      dimension hqqf(n,n), hqqt(ntri), diagqq(n)
      common/io/inp, iout 
      if(type.eq.'triangle') then
         call rzero(hqqt,ntri)
      else
         call rzero(hqqf,n*n)
      endif
      call dagfil(hqqf,hqqt,diagqq,n,type)
      if(nonz.gt.0) then
         if(nodsk) then
            call iosys('read integer '//header//' from '//
     #                 'hamiltonian without rewinding',
     #                  2*nonz,ibuf,0,' ')
            call iosys('read integer '//header//' from '//
     #                 'hamiltonian without rewinding',
     #                  wptoin(nonz),rbuf,0,' ')
            call hfill(ibuf,rbuf,hqqf,hqqt,n,nonz,type)
         else
            trips = nonz/lenbuf
            left = nonz - trips*lenbuf
            do 1000 trp=1,trips
               call iosys('read integer '//header//' from '//
     #                    'hamiltonian without rewinding',
     #                     2*lenbuf,ibuf,0,' ')
               call iosys('read integer '//header//' from '//
     #                    'hamiltonian without rewinding',
     #                     wptoin(lenbuf),rbuf,0,' ')
               call hfill(ibuf,rbuf,hqqf,hqqt,n,lenbuf,type)
 1000       continue   
            if(left.ne.0) then
               call iosys('read integer '//header//' from '//
     #                    'hamiltonian without rewinding',
     #                     2*left,ibuf,0,' ')
               call iosys('read integer '//header//' from '//
     #                    'hamiltonian without rewinding',
     #                     wptoin(left),rbuf,0,' ')
               call hfill(ibuf,rbuf,hqqf,hqqt,n,left,type)
            end if
         endif
      endif
      if(type.eq.'full') then
         do 2000 i=1,n
            do 3000 j=1,i
               hqqf(j,i)=hqqf(i,j)
 3000       continue
 2000    continue
         call iosys('write real H_QQ to hamiltonian',n*n,hqqf,0,' ')
      else
         call iosys('write real H_QQ to hamiltonian',ntri,
     $               hqqt,0,' ')
      endif   
      return
      end       
