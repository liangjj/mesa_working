*deck orth.f
c***begin prologue     orth
c***date written       950704   (yymmdd)
c***revision date               (yymmdd)
c***keywords           orthogonality
c***author             schneider, barry(nsf)
c***source             
c***purpose            check orthogonality of numerical eigenvectors
c***routines called
c***end prologue     orth
      subroutine orth(s,vec,wts,r,scr,nroots,n)
      implicit integer (a-z)
      dimension s(nroots,nroots), vec(n,nroots) ,wts(n), scr(n,nroots)
      dimension r(n)
      real*8 s, vec, wts, anorm, sdot, pi, k, r, scr, thresh
      real*8 snorm, sum
      character*80 title
      logical nowts
      common /io/ inp, iout
      parameter (thresh=1.d-08 )
      nowts=.false.
      title='input vectors in orth'
      call prntrm(title,vec,n,nroots,n,nroots,iout)
c----------------------------------------------------------------------c
c              schmidt orthonormalization of vectors                   c
c----------------------------------------------------------------------c
      nvec=1
      anorm=1.d0/sqrt(snorm(vec(1,1),vec(1,1),wts,n,nowts))
      call sscal(n,anorm,vec(1,1),1)
      do 10 i=2,nroots
         anorm=snorm(vec(1,i),vec(1,i),wts,n,nowts)
         if (anorm.ne.0.d+00) then
             anorm=1.d+00/sqrt(anorm)
         endif
         call sscal(n,anorm,vec(1,i),1)
         do 20 j=1,i
            if (i.ne.j) then
                sum=-snorm(vec(1,j),vec(1,i),wts,n,nowts)
                call saxpy (n,sum,vec(1,j),1,vec(1,i),1)
            endif
   20    continue
         anorm=snorm(vec(1,i),vec(1,i),wts,n,nowts)
         if (anorm.gt.thresh) then
             nvec=nvec+1
             anorm=1.d+00/sqrt(anorm)
             do 30 j=1,n
                vec(j,nvec)=anorm*vec(j,i)
   30        continue
         else
             write(iout,1) i, anorm 
         endif
   10 continue
      title='vectors'
      call prntrm(title,vec,n,nroots,n,nroots,iout)
      call copy(vec,scr,n*nroots)
      call vmmul(wts,scr,scr,n,nroots)
      call ebtc(s,vec,scr,nroots,n,nroots)
      title='overlap matrix in orth'
      call prntrm(title,s,nroots,nroots,nroots,nroots,iout)
      pi=3.141592653d0
      do 40 i=1,nroots
         k=(2*i-1)*pi/2.d0
         do 50 j=1,n
            vec(j,i)=sin(k*r(j))
 50      continue
 40   continue   
      do 60 i=1,nroots
         anorm=1.d0/sqrt(snorm(vec(1,i),vec(1,i),wts,n,nowts))
         call sscal(n,anorm,vec(1,i),1)
 60   continue   
      call copy(vec,scr,n*nroots)
      call vmmul(wts,scr,scr,n,nroots)
      call ebtc(s,scr,vec,nroots,n,nroots)
      title='vectors'
      call prntrm(title,vec,n,nroots,n,nroots,iout)
      title='overlap matrix'
      call prntrm(title,s,nroots,nroots,nroots,nroots,iout)
      return
 1    format(/,1x,'vector = ',i3,' below threshold. norm = ',e15.8)
      end

