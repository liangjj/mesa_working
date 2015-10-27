*deck @(#)sethpp.f	5.1  11/6/94
      subroutine sethpp(hpp,hss,grad,diag,t,mdim,nrhs,c0,npvec,
     $                  npdim,scrtch,incore,lnbuf,energy,
     $                  nbins,binsiz,nblks,blksiz,tsize,mxblk,lsize,
     $                  bptr,xbin,ibin,tvec,iu,iv,iprt,tstopt,freeze,
     $                  hqq,filtyp)
c
      implicit real*8(a-h,o-z)
      real*8 hss(*),diag(mdim),t(mdim,nrhs),tvec(*)
      real*8 c0(mdim,npvec),grad(mdim,nrhs),scrtch(*)
      real*8 hpp(npvec,npvec),xbin(*)
      character*(*) filtyp
      integer ibin(*),bptr(*)
      integer binsiz,blksiz,tsize,binsz2
      integer tstopt,freeze,twalks,hqq
      logical debug
c
      parameter (debug=.false.)
      data small/1.0d-9/
      save small
c
      common /io/inp,iout
c
c
      binsz2=binsiz+2
      twalks=mdim-npdim
      if(incore.eq.1) then
         call rdham(bptr,nbins,xbin,ibin,hss,binsiz,
     $              binsz2,tblks,nblks,lsize,blksiz,mxblk,iu,iprt,
     $              diag,mdim,energy,filtyp)
      else
         call tdiag(bptr,nbins,xbin,ibin,hss,binsiz,
     $              binsz2,tblks,nblks,lsize,blksiz,mxblk,iu,iprt,
     $              diag,mdim,energy,filtyp)
      end if
c
      mtot=mdim*npvec
c     
      if(tstopt.ne.0) then
         call iosys('read real "p-space vectors" from '//filtyp,
     $               mtot,c0,0,' ')
      else
         call rzero(c0,mtot)
         do 1111 i=1,npvec
            c0(i,i)=1.d0
 1111    continue
      endif
c
c     ----- construct hpp and start hpq -----
c
      call tmult(bptr,nbins,xbin,ibin,hss,binsiz,
     $            binsz2,tblks,nblks,lsize,blksiz,mxblk,iu,iprt,
     $            c0,grad,tvec,mdim,npvec,incore,energy,filtyp)
c
      call ebtc(hpp,c0,grad,npvec,mdim,npvec)
c
c     ----- add projection operator contributions to the diagonals -----
c
      if (hqq.eq.0) then
          if(mdim.ne.npvec) then
             do 8 i=1,npvec
                do 9 j=1,mdim
                   diag(j)=diag(j)-2.d0*grad(j,i)*c0(j,i)
   9            continue
   8         continue
          end if
c
c         -- t(j,i) = c0(j,k)*hpp(k,i)  .. skipping zeros in c0 --
c
         call mxma(c0,1,mdim,hpp,1,npvec,t,1,mdim,
     $             mdim,npvec,npvec)
c
c        -- finish projection operator contributions to the diagonal --
c
          if(npvec.ne.mdim) then
             do 11 j=1,npvec
                do 10 i=1,mdim
                   diag(i)=diag(i)+c0(i,j)*t(i,j)
  10            continue
  11         continue
             do 12 i=1,mdim
                if(abs(diag(i)).lt.small) diag(i)=1.d0
  12         continue
          end if
      endif
c
c     -- finish hpq --
c
      if(freeze.eq.0) then
         call project(grad,mdim,nrhs,c0,npvec,scrtch,t)
      else
         call zap(grad,mdim,twalks,nrhs)
      end if
c
c
      if(debug) then
         write(iout,20)
         call matout(hss,mdim,mdim,mdim,mdim,iout)
         write(iout,21)
         call matout(c0,mdim,nrhs,npdim,nrhs,iout)
         write(iout,22)
         call matout(grad,mdim,nrhs,mdim,nrhs,iout)
      endif
c
 20   format(//,'  hss ',/)
 21   format(//,'  c0 ',/)
 22   format(//,'  rhs ',/)
c
      return
      end
