*deck @(#)jmatrix.f	1.1  4/25/95
      subroutine jmatrix(values,d,nbf,nnp,jmat,ncoul,
     $                   ntriang,ndmat)
c
c***begin prologue     jmatrix.f
c***date written       870521   (yymmdd)  
c***revision date      11/6/94      
c
c***keywords           j matrix, coulomb matrix 
c***author             saxe, paul and martin, richard (lanl)
c***source             @(#)jmatrix.f	1.1   4/25/95
c***purpose            to form the j matrix given integrals and
c                      the density matrices. 
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       jmatrix.f
      implicit none
c     --- input variables ---
      integer nbf,nnp,ncoul,ndmat,ntriang
c     --- input arrays (unmodified) ---
      real*8 d(nnp,ndmat)
c     --- input arrays (scratch) ---
      real*8 values(nnp,ntriang)
c     --- output arrays ---
      real*8 jmat(nnp,ncoul)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer maxkl,minkl,kl,k,k1,l,l1
      integer lnread,shell
      real*8 zero,half,one,two
      logical called,read
c
      parameter (zero=0.0d+00,half=0.5d+00,one=1.0d+00,two=2.0d+00)
c
      data called /.false./
      save called
c
      common /io/ inp,iout
c
c     --- rewind the integrals, initalize the j matrix, and let's go ---
      called=.false.
      call iosys('rewind "sorted ao integrals" on ints',0,0,0,' ')
      call rzero(jmat,nnp*ncoul)
c
c     --- read through integrals if first iteration, or if
c           cannot hold all the integrals in core ---
      if (.not.called.or.ntriang.lt.nnp) then
         maxkl=0
         read=.true.
      else
         minkl=1
         maxkl=nnp
         read=.false.
      end if
c
      kl=0
      do 100 k=1,nbf
         k1=k*(k-1)/2+1
         do 90 l=1,k
            l1=l*(l-1)/2+1
            kl=kl+1
c
c           --- check that this triangle of integrals is in core ---
            if (kl.gt.maxkl) then
               minkl=maxkl+1
               maxkl=min(nnp,maxkl+ntriang)
               lnread=(maxkl-minkl+1)*nnp
               call iosys('read real "sorted ao integrals" from ints '//
     $                     'without rewinding',lnread,values,0,' ')
            end if
c
c            --- form coulomb terms ---
            if (k.eq.l.and.read) then
               call sscal(nnp,half,values(1,kl-minkl+1),1)
            end if
c
            if (kl.eq.maxkl) then
               do 80 shell=1,ncoul
                  call sgemv('n',nnp,maxkl-minkl+1,one,values,nnp,
     $                       d(minkl,shell),1,one,jmat(1,shell),1)
   80          continue
            end if
   90    continue
  100 continue
c
c     --- double coulomb terms to account for [ij;kl]
c                   and [ij;lk] ---
      do 110 shell=1,ncoul
         call sscal(nnp,two,jmat(1,shell),1)
  110 continue
c
c
      return
      end
