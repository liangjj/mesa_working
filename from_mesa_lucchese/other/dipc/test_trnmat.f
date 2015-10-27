      subroutine trnmat(hpvbx,hpvby,hpvbz,hpvbtx,hpvbty,hpvbtz,
     $hpvb,hpvbt,hbx,hby,hbz,hbtx,tr,nchan2,nlm,
     $ngauss,nmotot,nchnl,nbfmax,lm2,iene,lmtop)
      implicit real*8(a-h,o-z)
      real*8 hpvbx(lm2,nbfmax,nchnl)
     $,hpvbz(lm2,nbfmax,nchnl),hpvbtx(lm2,nbfmax,nchnl)
     $,hpvbtz(lm2,nbfmax,nchnl)
     $,hpvb(lm2,nbfmax,nchnl),hpvbt(lm2,nbfmax,nchnl)
      complex*16 hpvby(lmtop,nbfmax,nchnl), hpvbty(lmtop,nbfmax,nchnl)
     $,cterm
      real*8 hbx(nbfmax,nbfmax),hby(nbfmax,nbfmax),
     $hbz(nbfmax,nbfmax),hbtx(nbfmax,nbfmax),tr(nbfmax,nbfmax)
      integer nlm(nchnl)
      do 1 ic=1,nchan2
      nl=nlm(ic)
      call mxma(hpvbx(1,1,ic),2,lm2,tr,1,nbfmax,hpvbtx(1,1,ic),
     $          2,lm2,nl,ngauss,nmotot)
      call mxma(hpvbx(2,1,ic),2,lm2,tr,1,nbfmax,hpvbtx(2,1,ic),
     $          2,lm2,nl,ngauss,nmotot)
c      call mxma(hpvby(1,1,ic),2,lm2,tr,1,nbfmax,hpvbty(1,1,ic),
c     $          2,lm2,nl,ngauss,nmotot)
c      call mxma(hpvby(2,1,ic),2,lm2,tr,1,nbfmax,hpvbty(2,1,ic),
c     $          2,lm2,nl,ngauss,nmotot)
      do i=1,nl
         do j=1,nmotot
            cterm=0.d0
            do k=1,ngauss
               cterm=cterm+hpvby(i,k,ic)*tr(k,j)
            enddo
            hpvbty(i,j,ic)=cterm
         enddo
      enddo
      call mxma(hpvbz(1,1,ic),2,lm2,tr,1,nbfmax,hpvbtz(1,1,ic),
     $          2,lm2,nl,ngauss,nmotot)
      call mxma(hpvbz(2,1,ic),2,lm2,tr,1,nbfmax,hpvbtz(2,1,ic),
     $          2,lm2,nl,ngauss,nmotot)
      call mxma(hpvb(1,1,ic),2,lm2,tr,1,nbfmax,hpvbt(1,1,ic),
     $          2,lm2,nl,ngauss,nmotot)
      call mxma(hpvb(2,1,ic),2,lm2,tr,1,nbfmax,hpvbt(2,1,ic),
     $          2,lm2,nl,ngauss,nmotot)
 1    continue
      if(iene.eq.1)then
         do i=1,ngauss
            do j=1,nmotot
               term=0.d0
               do k=1,ngauss
                  term=term+hbx(i,k)*tr(k,j)
               enddo
               hbtx(i,j)=term
            enddo
         enddo
         do i=1,nmotot
            do j=1,nmotot
               term=0.d0
               do k=1,ngauss
                  term=term+tr(k,i)*hbtx(k,j)
               enddo
               hbx(i,j)=term
            enddo
         enddo
         do i=1,ngauss
            do j=1,nmotot
               term=0.d0
               do k=1,ngauss
                  term=term+hby(i,k)*tr(k,j)
               enddo
               hbtx(i,j)=term
            enddo
         enddo
         do i=1,nmotot
            do j=1,nmotot
               term=0.d0
               do k=1,ngauss
                  term=term+tr(k,i)*hbtx(k,j)
               enddo
               hby(i,j)=term
            enddo
         enddo
         do i=1,ngauss
            do j=1,nmotot
               term=0.d0
               do k=1,ngauss
                  term=term+hbz(i,k)*tr(k,j)
               enddo
               hbtx(i,j)=term
            enddo
         enddo
         do i=1,nmotot
            do j=1,nmotot
               term=0.d0
               do k=1,ngauss
                  term=term+tr(k,i)*hbtx(k,j)
               enddo
               hbz(i,j)=term
            enddo
         enddo
      endif
      return
      end
