*deck @(#)optdir.f	1.4  8/3/91
c
c***begin prologue     optdir
c***date written       000811   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           ci, hamiltonian matrix, contracted
c***                   basis
c***author             schneider, barry(nsf)
c***source             @(#)m942.f	1.4   8/3/91
c
c***purpose            to compute the optical potential 
c***                   when H_QQ and H_QP are in core.
c
c***references
c
c***routines called    (none)
c
c***end prologue       optdir
c
c
      subroutine optdir(diag,hqq,hqp,xqp,opt,tmp,buf,ibuf,energy,ipvt,
     #                  header,lenbuf,n,npvec,ntri,nen,nonzqq,nonzqp,
     #                  nodsk,prnt)
      implicit integer (a-z)
c
      character*(*) header(10)
      real*8 hqq, hqp, buf, opt, diag, xqp, tmp
      real*8 energy
      logical nodsk, prnt
      dimension diag(n), hqq(ntri), hqp(n,npvec), buf(lenbuf)
      dimension ibuf(2,lenbuf), xqp(n,npvec), energy(nen)
      dimension opt(npvec,npvec), ipvt(npvec), prnt(3)
      dimension tmp(ntri)
      common/io/inp,iout
c
c     rewind the optical potential file which now holds
c     the 'static-exchange" part of the hamiltonian.  that
c     is (H_PP - Energy)
c
      call iosys('rewind '//header(7)//' on hamiltonian',
     #            0,0,0,' ')
c
c
c     read in the hqq diagonals
c
      call iosys('read real '//header(8)//
     #           ' from hamiltonian',n,diag,0,' ')
c
c
c     fill all of hqq
c
      call rdhqq(ibuf,buf,tmp,tmp,diag,header(9),
     #           n,ntri,lenbuf,nonzqq,nodsk,
     #           'triangle',prnt(2))
c
c     read in the right hand sides, hqp
c
      call rdhqp(ibuf,buf,hqp,header(5),n,npvec,lenbuf,
     #           nonzqp,nodsk,prnt(1))
      do 10 i=1,nen
         call copy(tmp,hqq,ntri)
c
c        make Energy - H_QQ
c
         call eshft(hqq,hqq,energy(i),n,ntri,'triangle')
c
c        copy the right hand side into xqp
c
         call copy(hqp,xqp,npvec*n)
c
c        solve the linear equations, (Energy -H_QQ) X_QP = H_QP
c
         call dspsv('u',n,npvec,hqq,ipvt,xqp,n,info)
c
c                                                 -1
c        form (H_PP - Energy + H_PQ (Energy - H_QQ)  H_QP
c
         call frmopt(opt,xqp,hqp,header(7),energy(i),n,
     #               npvec,prnt(3))
 10   continue   
      return
      end
