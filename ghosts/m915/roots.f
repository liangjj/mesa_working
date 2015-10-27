*deck @(#)roots.f	2.1  10/10/91
      subroutine roots(ints,v,s,root,dvdvec,dvdmat,ndvdmx,prtflg)
c
      implicit integer (a-z)
c
      real*8 ints(*),v(nwks,mxvc),s(nwks,mxvc),root(nroots)
      real*8 dvdvec(*),dvdmat(*)
      real*8 rep,fzcore,eguess,eci,cnverg,sqcdif,czero,edav
      character*20 status,itoc*4
      character*8 prtflg
c
      common /io/     inp,iout
      common /dims/   nbf,nsym,norbs,nrows,nrows4,nwks,nwks2,nlevs,
     #                nrowoc,levfrm,nwksmx,nlwkmx,nuwkmx,bmax,nroots,
     #                orbfrm
      common /d901/   rep,fzcore,eguess,eci,refwlk,mxiter,cnverg,icnvg,
     #                iter,sqcdif,czero,nroot
      common /nvectr/ nvc,mxvc
c
c
      write(iout,*)'  roots:  nuclear repulsion energy ',rep
c
   10 continue
         call bliu('solve',status,v,s,nwks,mxiter,0,0,0,dvdmat,root,
     #             dvdvec)
         if (status.eq.'converged') go to 100
c
         nvc=0
   20    continue
            nvc=nvc+1
            call bliu('new trial',status,v(1,nvc),s(1,nvc),nwks,
     #                mxiter,0,0,0,0,root,dvdvec)
         if (status.ne.'done'.and.nvc.lt.mxvc) go to 20
            if (status.eq.'done'.and.nvc.eq.1) go to 10
            if (status.eq.'done') nvc=nvc-1
c
            call rzero(s,nwks*nvc)
            call loopy
            do 30 ivc=1,nvc
c     write (iout,83) ivc,(v(iq,ivc),iq=1,nwks)
c  83 format ('1'//,' vector ivc=',i5,/,(1x,10f12.6))
c     write (iout,84) ivc,(s(iq,ivc),iq=1,nwks)
c  84 format (//,'  sigma, ivc=',i5,/,(1x,10f12.6))
               call bliu('with vectors',0,v(1,ivc),s(1,ivc),nwks,
     #                   mxiter,0,0,0,0,0,0)
   30       continue
            if (status.eq.'done') go to 10
            nvc=0
         go to 20
c
  100 continue
c
c     ----- recover each root's ci vector, find the most important
c           configuration, and write vector to rwf
c
      call iosys('read real diagonals from bliu',-1,s,0,' ')
c
      if (prtflg.ne.'minimum') then
      write(iout,*)'  roots:  nuclear repulsion energy ',rep
      write (iout,90)
   90 format (///,' root   reference  guess energy    ci energy  ',
     #        '  davidson energy  c(0)')
      end if
c
      do 200 iroot=1,nroots
         call bliu('get vector',status,v,0,nwks,0,iroot,0,0,0,0,0)
         if (status.ne.'ok') go to 200
         call iosys('write real "ci root '//itoc(iroot)//'" to rwf',
     #               nwks,v,0,' ')
c
         czero=0.0d+00
         refwlk=0
         do 110 i=1,nwks
            if (abs(v(i,1)).gt.abs(czero)) then
               czero=v(i,1)
               refwlk=i
            end if
  110    continue
c
         eguess=s(refwlk,1)+rep+fzcore
         eci=root(iroot)+rep+fzcore
         edav=eci+(eci-eguess)*(1.0d+00-czero**2)
         if (iroot.eq.1) then
            call iosys('write real energy to rwf',1,eci,0,' ')
         end if
c
      if (prtflg.ne.'minimum') then
         write (iout,120) iroot,refwlk,eguess,eci,edav,czero
  120    format (1x,i3,i10,3g18.9,f8.4)
      end if
c
  200 continue
c
      call bliu('finish',0,0,0,0,0, 0,0,0,0,0,0)
c
c
      return
      end
