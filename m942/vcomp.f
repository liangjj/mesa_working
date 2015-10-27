*deck @(#)vcomp
      subroutine vcomp(cvec,pvec,s,eig,work,scr,nwks,nwksp,
     #                 nwksq,npvec,ncomp,prnt)
      implicit integer(a-z)
      real*8 cvec(nwks,ncomp), pvec(nwks,npvec), s(nwksp,nwksp)
      real*8 eig(nwksp), work(5*nwksp), scr(nwks,nwksp)
      real*8 fac, nrzero
      character*80 title
      logical prnt
      common/io/inp, iout
      data nrzero / 1.d-08 /
c
      if(prnt) then
         title='input vectors'
         call prntrm(title,pvec,nwks,npvec,nwks,npvec,iout)
      endif
      call rzero(scr,nwks*nwksp)
      do 10 i=1,nwksp
         scr(nwksq+i,i) = 1.d0
 10   continue
      call ambctxx(scr(nwksq+1,1),pvec(nwksq+1,1),pvec(nwksq+1,1),
     #             nwksp,npvec,nwksp,nwks,nwks,nwks)
      if(prnt) then
         title='projected vectors'
         call prntrm(title,scr,nwks,nwksp,nwks,nwksp,iout)
      endif
      call rzero(s,nwksp*nwksp)
      do 20 i=1,nwksp
         do 30 j=1,i
            do 40 k=1,nwksp
               s(i,j) = s(i,j) + scr(nwksq+k,i)*scr(nwksq+k,j)
 40         continue
            s(j,i) = s(i,j)
 30      continue
 20   continue
      if(prnt) then
         title='overlap of projected vectors'
         CALL prntrm(title,s,nwksp,nwksp,nwksp,nwksp,iout)
      endif
      CALL dsyev('v','l',nwksp,s,nwksp,eig,work,5*nwksp,info)
      title='eigenvalues of overlap matrix'
      CALL prntrm(title,eig,nwksp,1,nwksp,1,iout)
c      if(prnt) then
c         title='eigenvectors'
c         CALL prntrm(title,s,nwksp,nwksp,nwksp,nwksp,iout)
c      endif
      count=0
      do 50 i=1,nwksp
         if(abs(eig(i)).gt.nrzero) then
            count=count+1
            fac=1.d0/sqrt(eig(i))
            do 60 j=1,nwksp
               s(j,count) = fac * s(j,i)
 60         continue
         endif
 50   continue
      if(prnt) then
         title='eigenvectors'
         CALL prntrm(title,s,nwksp,count,nwksp,count,iout)
      endif
c      call rzero(cvec,nwks*ncomp)
      call ebcxx(cvec(nwksq+1,1),scr(nwksq+1,1),s,nwksp,nwksp,count,
     #           nwks,nwks,nwksp)
      if(prnt) then
         title='eigenvectors'
         CALL prntrm(title,cvec,nwks,count,nwks,count,iout)
      endif
      return
      end
