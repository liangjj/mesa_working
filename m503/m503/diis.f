*deck  @(#)diis.f	5.1 11/6/94
      subroutine diis(f,t1,t2,t3,t4,t5,triang,error,
     $                diissv,energs,bmatrx,cicoef,
     $                focksv,c,cdiis,diiser,energy,stdiis,iprint,
     $                nnp,nbf,mxds,mxiter,iter,calc,
     $                ndiis,mndiis,ptr,dsptr,
     $                prnt,first,ops)
c
c     module to perform pulay's diis extrapolation.
c
      implicit integer (a-z)
c
      character*(*) ops
      character*(*) calc
      real*8 f(nnp),t1(nbf,nbf),t2(nbf,nbf),t3(nbf,nbf),t4(nbf,nbf)
      real*8 t5(nbf,nbf),triang(nnp),error(nbf,nbf,mxds)
      real*8 diissv(mxiter),energs(mxiter)
      real*8 bmatrx(*),cicoef(2)
      real*8 focksv(nnp,mxds)
      real*8 c(nbf,nbf),cdiis(nbf,nbf)
      real*8 diiser,energy,sdot,scalar,ten
      integer ptr(mxiter),dsptr(mxds)
      logical prnt,logkey,first
c
      common /io/     inp,iout
c
      parameter (ten=10.0d+00)
c
c
 1000 format (5x,'diis extrapolation in use')
c
c     ----- note t1 and bmatrx are implicitly equivalenced.
c
c
c     ----- compute the diis error vector fds-sdf -----
c
      call trtosq(t2,triang,nbf,nnp)
      call ebc(t4,t1,t2,nbf,nbf,nbf)
      call ebc(t5,t4,t3,nbf,nbf,nbf)
      call ebc(t4,t3,t2,nbf,nbf,nbf)
      call ambc(t5,t4,t1,nbf,nbf,nbf)
c
c
      if (logkey(ops,'scf=print=diis-error',.false.,' ')) then
            write (iout,190)
  190       format (//,' the diis vector in the a. o. basis',/)
            call matout(t5,nbf,nbf,nbf,nbf,iout)
      end if
c
c     ----- find maximum element of error vector -----
c
      diiser=0.0d+00
      do 12 i=1,nbf
         do 11 j=1,i
            diiser=max(diiser,abs(t5(i,j)))
   11    continue
   12 continue
c
c
      diissv(iter)=diiser
c
c     ----- print results if needed each iteration -----
c
      if (prnt.and.logkey(ops,'print=scf=energy',.false.,' ')) then
         if (calc.eq.'gvb') then
            write (iout,5001) iter,energy,diiser,cicoef(1),cicoef(2)
 5001       format (5x,i4,2(5x,f15.9),2x,2f11.4)
         else
            write (iout,31) iter,energy,diiser
   31       format(5x,i4,2(5x,f15.9))
         end if
      end if
c
c     ----- check if error small enough for diis procedure -----
c
      if (diiser.lt.ten**(-stdiis)) then
         if (first) then
            write(iout,1000)
            first=.false.
c
c           ----- save the current mo's in 'cdiis' as a basis for diis
c
            call vmove(cdiis,c,nbf**2)
         end if
c
c        ----- transform to the ao basis -----
c
c        call trtosq(t1,triang,nbf,nnp)
c        call ebct(t2,t1,c,nbf,nbf,nbf)
c        call ebc(t5,c,t2,nbf,nbf,nbf)
c
         if (logkey(ops,'scf=print=ao-diis-error',.false.,' ')) then
            write (iout,110)
  110       format (//,' the diis vector in the a.o. basis',/)
            call matout(t5,nbf,nbf,nbf,nbf,iout)
         end if
c
c        ----- transform the error vector to original mo basis -----
c
         call ebc(t2,t5,cdiis,nbf,nbf,nbf)
         call ebtc(t3,cdiis,t2,nbf,nbf,nbf)
c
         if (logkey(ops,'scf=print=mo-diis-error',.false.,' ')) then
            write (iout,111)
  111       format (//,' the diis error vector in the original ',
     $              'mo basis',/)
            call matout(t3,nbf,nbf,nbf,nbf,iout)
         end if
c
         ndiis=ndiis+1
c
c        ----- find an available diis storage spot -----
c
         do 701 i=1,mxds
            if (dsptr(i).le.0) then
               dsptr(i)=ndiis
               ptr(ndiis)=i
               go to 702
            end if
  701    continue
c
c        ---- none available, so use oldest one -----
c
         i=ptr(mndiis)
         ptr(mndiis)=-9999998
         dsptr(i)=ndiis
         ptr(ndiis)=i
         mndiis=mndiis+1
c
  702    continue
c
         call vmove(error(1,1,ptr(ndiis)),t3,nbf**2)
         call vmove(focksv(1,ptr(ndiis)),f,nnp)
         if (ndiis.ge.2) then
c
c     ----- form small matrix -----
c
  100       continue
            bmatrx(1)=0.0d+00
            n=1
            do 13 i=mndiis,ndiis
               n=n+1
               bmatrx(n)=-1.0d+00
   13       continue
            scalar=sdot(nbf**2,error(1,1,ptr(mndiis)),1,
     $                 error(1,1,ptr(mndiis)),1)
            do 15 i=mndiis,ndiis
               n=n+1
               bmatrx(n)=-1.0d+00
               do 14 j=mndiis,ndiis
                  n=n+1
                  bmatrx(n)=sdot(nbf**2,error(1,1,ptr(i)),1,
     $                           error(1,1,ptr(j)),1)/scalar
   14          continue
   15       continue
            n=n+1
            bmatrx(n)=-1.0d+00
            do 16 i=mndiis,ndiis
               n=n+1
               bmatrx(n)=0.0d+00
   16       continue
c
c     ----- solve small set of equations -----
c
            if (logkey(ops,'scf=print=diis-matrix',.false.,' ')) then
               write (iout,112)
  112          format (//,' small matrix for diis ',/)
               call matout(bmatrx,ndiis-mndiis+2,ndiis-mndiis+3,
     $                     ndiis-mndiis+2,ndiis-mndiis+3,iout)
            end if
c
            call flin(bmatrx,ndiis-mndiis+2,ndiis-mndiis+2,1,scalar)
c
c           write (iout,17) ndiis,scalar
c   17      format (9x,' diis iteration',i4,', the determinant is ',
c     #             e9.1)
c
c     ----- check if problems solving the equations -----
c
            if (abs(scalar).lt.1.0d-15) then
               i=ptr(mndiis)
               dsptr(i)=-1
               ptr(mndiis)=-9999998
               mndiis=mndiis+1
               if (mndiis.eq.ndiis-1) go to 101
               go to 100
            end if
c
c     ----- form new fock matrix -----
c
            call rzero(f,nnp)
            n=(ndiis-mndiis+2)**2+1
c
c           if (btest(iprint,14)) then
c              write (iout,114) (iq,ptr(iq),bmatrx(n+ptr(iq)),
c     #               iq=mndiis,ndiis)
c  114         format ( //,' solution vector for diis ',/,
c     #               (1x,2i4,f15.9))
c           end if
c
            do 18 i=mndiis,ndiis
               n=n+1
               call saxpy(nnp,bmatrx(n),focksv(1,ptr(i)),1,f,1)
   18       continue
c
c           if (btest(iprint,15)) then
c              write (iout,115)
c  115         format (//,' fock matrix after diis ',/)
c              call print(f,nnp,nbf,iout)
c           end if
         end if
      end if
c
  101 continue
c
c
c
c
c
      return
      end
