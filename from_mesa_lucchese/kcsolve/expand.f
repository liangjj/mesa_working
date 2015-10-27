      subroutine expand(htopp,store,tq,nbig,nsmall,nofree,notot,
     $  nqvecs,nq)
      implicit real*8 (a-h,o-z)
      complex*16 htopp(nbig,nsmall),store(nqvecs,nofree),cterm
      real*8 tq(nq,nqvecs)
      do i=1,nofree
         do j=1,nq
            cterm=0.d0
            do k=1,nqvecs
               cterm=cterm+store(k,i)*tq(j,k)
            enddo
            htopp(j+notot,i)=cterm
         enddo
      enddo
      return
      end
