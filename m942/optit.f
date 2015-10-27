*deck @(#)optit.f	1.4  8/3/91
c
c***begin prologue     optit
c***date written       000811   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           ci, hamiltonian matrix, contracted
c***                   basis
c***author             schneider, barry(nsf)
c***source             @(#)m942.f	1.4   8/3/91
c
c***purpose            to compute the optical potential 
c***                   when H_QQ cannot be stored in core.
c
c***references
c
c***routines called    (none)
c
c***end prologue       optit
c
c
      subroutine optit(diag,hqp,opt,buf,ibuf,energy,vec,hvec,
     #                 resid,mat,mattmp,b,btmp,t,list,cnverg,thresh,n,
     #                 nrhs,nen,nonzqq,nonzqp,maxit,maxvec,lenbuf,
     #                 incore,header,prnt,prntdvd)
      implicit integer (a-z)
c
      character*(*) header(10)
      character*4 itoc
      real*8 diag, hqp, opt, buf, energy
      real*8 vec, hvec, resid, mat, mattmp, b, btmp, t
      real*8 cnverg, thresh
      logical incore, prnt, prndvd
      dimension diag(n), hqp(n,nrhs), opt(nrhs,nrhs)
      dimension buf(lenbuf), ibuf(2,lenbuf), energy(nen)
      dimension vec(n,maxvec), hvec(n,maxvec), resid(n,nrhs)
      dimension mat(maxvec,maxvec), mattmp(maxvec,maxvec)
      dimension b(maxvec,nrhs), btmp(maxvec,nrhs), t(n,nrhs)
      dimension list(nrhs), prnt(3), prndvd(11)
      common/io/inp,iout
c
c     rewind the optical potential file which now holds
c     the 'static-exchange" part of the hamiltonian.  that
c     is (H_PP - Energy)
c
      call iosys('rewind '//header(7)//' on hamiltonian',
     #            0,0,0,' ')
c
c     read in the H_QQ diagonals
c
      call iosys('read real '//header(8)//' from hamiltonian',n,diag,
     #            0,' ')
c
      if(nonzqq.ne.0) then 
         if(incore) then
            call iosys('rewind '//header(9)//' on hamiltonian '//
     #                 'read-and-write',0,0,0,' ')
            call iosys('read integer '//header(9)//' from '//
     #                 'hamiltonian without rewinding',
     #                  2*nonzqq,ibuf,0,' ')
            call iosys('read integer '//header(9)//' from '//
     #                 'hamiltonian without rewinding',
     #                  wptoin(nonzqq),buf,0,' ')
         endif
      endif
c
c     read in the right hand sides, H_QP
c
      call rdhqp(ibuf,buf,hqp,header(5),n,nrhs,lenbuf,
     #           nonzqp,incore,prnt(2))
c
      do 10 i=1,nen
c
c        solve the linear equations, (H_QQ _ Energy) X_QP = H_QP
c
         call lindvd(ibuf,rbuf,diag,energy(i),vec,hvec,
     #               hqp,mat,mattmp,b,btmp,resid,t,list,cnverg,thresh,
     #               n,nrhs,maxit,maxvec,lenbuf,incore,header(8),
     #               prndvd)
c
c        read in the solutions into t
c
           do 20 j=1,nrhs
              call iosys('read real "solution for right hand side = '
     #                    //itoc(list(j))//'" from hamiltonian',n,
     #                      t(1,j),0,' ')
 20        continue   
c                                                 -1
c        form (H_PP - Energy + H_PQ (Energy - H_QQ)  H_QP
c
         call frmopt(opt,t,hqp,header(7),energy(i),n,
     #               nrhs,prnt(3))

 10   continue   
      return
      end
