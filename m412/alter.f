*deck @(#)alter.f	5.1  11/6/94
      subroutine alter(c,eigval,scr,num)
c***begin prologue     alter
c***date written       850601  yymmdd
c***revision date      yymmdd  yymmdd
c***keywords           $alter, alter, guess
c***author             martin, richard (lanl)
c***source             @(#)alter.f	5.1   11/6/94
c***purpose            alters initial guess eigenvectors and eigenvalues.
c***description
c     call alter(c,eigval,scr,num)
c       c       eigenvector array(num,num) ordered by columns.
c       eigval  eigenvalue array(num).
c       scr     scratch array(num).
c       num     number of orbitals.
c***references
c***routines called    positn(chr), lnkerr(mdutil), ffnext(chr),
c                      ctoi(chr), vmove(math)
c***end prologue       alter
      implicit integer(a-z)
      real*8 c(num,num),eigval(num),scr(num)
      real*8 scrval
      character card*80,ffnext*16,found*16
      logical positn
c
      common/io/inp,iout
c
 1000 format(a)
 1010 format(5x,'orbital',i4,' switched with orbital',i4)
 1020 format(' alter pair:',2i4,'   number of orbitals:',i4)
c
c     alter the order of the eigenvectors.
c
c     abort if alteration section not present.
      if(.not.positn('$alter',card,inp))
     $   call lnkerr('no $alter section found.')
c
   20 read(inp,1000,end=9001) card
         if(index(card,'$').ne.0) return
         pos=0
   25    found=ffnext(card,pos,start,end)
            if(found.eq.'eos') goto 20
c
c
            if(found.ne.'integer')
     $         call lnkerr('alter switches must be integer. i found:'
     $                      //card)
            v1=ctoi(card(start:end))
c
            found=ffnext(card,pos,start,end)
            if(found.eq.'eos')
     $         call lnkerr('alter switches must be in pairs. i found:'
     $                     //card)
            if(found.ne.'integer')
     $         call lnkerr('alter switches must be integer. i found:'
     $                      //card)
            v2=ctoi(card(start:end))
c
c           check for switches out of range.
            if(v1.gt.num.or.v2.gt.num) then
               write(iout,1020) v1,v2,num
               call lnkerr('orbital index in an alter pair is '
     $             //'greater than the number of orbitals.')
            endif
c
c           move the vectors.
            call vmove(scr,c(1,v1),num)
            call vmove(c(1,v1),c(1,v2),num)
            call vmove(c(1,v2),scr,num)
c
c           move the eigenvalues.
            scrval=eigval(v1)
            eigval(v1)=eigval(v2)
            eigval(v2)=scrval
            write(iout,1010) v1,v2
            goto 25
c
c
 9001 call lnkerr('no terminus found for the $alter section.')
c
      end
