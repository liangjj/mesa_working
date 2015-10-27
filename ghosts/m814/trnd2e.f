*deck @(#)trnd2e.f	1.1  11/30/90
      subroutine trnd2e(c,values,t1,t2,mo,lab,nbf,nnp,nactiv,ncore,
     #                  nnpact,ntriang,nbfna,asort,lnsort,nder,third)
c
c   4 december 1986   pws at lanl
c      changing iosys open to character unit name.
c
      implicit integer (a-z)
c
      real*8 c(nbf,nbf),values(nnp,ntriang),t1(nbf,nbf),t2(nbf,nbf)
      real*8 third(nbfna,nnpact),asort(lnsort),mo(nnpact,nnpact)
      integer lab(nnp)
      logical prnt
c
      common /io/ inp,iout
c
      prnt=.false.
c
c     ----- open a scratch file -----
c
      junk=nnp**2/50
      call iosys('open scratch as scratch on ssd',junk,0,0,' ')
c
c     ----- loop through degrees of freedom, transforming as we go
c
      call iosys('rewind "sorted ao derivative integrals"'//
     #' on ints',0,0,0,' ')
      call iosys('create real der_i(xbcd) on ints',
     #            nbf*nactiv*nnpact*nder,0,0,' ')
      call iosys('create real der_g on ints',nnpact**2*nder,0,0,' ')
c
      do 1000 der=1,nder
c
c        ----- initialize the sorting routines -----
c
         call sorter('start',asort,asort,lnsort,nnp*nnpact,5120,0,0,0,
     #               'half','scratch',prnt)
c
c        ----- initialize the labels -----
c
         do 10 i=1,nnpact
            lab(i)=(i-1)*nnp
   10    continue
c
         maxkl=0
         kl=0
         do 200 k=1,nbf
            do 100 l=1,k
               kl=kl+1
c
c              ----- check that these integrals are in core -----
c
               if (kl.gt.maxkl) then
                  minkl=maxkl+1
                  maxkl=min(nnp,maxkl+ntriang)
                  lnread=(maxkl-minkl+1)*nnp
      call iosys('read real "sorted ao derivative integrals"'//
     #' from ints without rewinding',lnread,values,0,' ')
               end if
c
c              ----- first half-transformation -----
c
               call trtosq(t1,values(1,kl-minkl+1),nbf,nnp)
c
c..bhl
c               if(der.eq.1) then
c                write(iout,*)' der_g ao  k,l ',k,l
c                call matout(t1,nbf,nbf,nbf,nbf,iout)
c               endif
c..bhl
c
               call ebtc(t2,c(1,ncore+1),t1,nactiv,nbf,nbf)
               call ebc (t1,t2,c(1,ncore+1),nactiv,nbf,nactiv)
c
c              ----- increment labels and pass integrals to sort -----
c
               do 20 i=1,nnpact
                  lab(i)=lab(i)+1
   20          continue
               call sqtotr(t2,t1,nactiv,nnpact)
c
               call sorter('with bin',asort,asort,0,nnpact,lab,t1,t2,
     #                      0,0,0,prnt)
c
  100       continue
  200    continue
c
c        ----- finish the sorting of the half-transformed integrals
c              to [kl;ab]
c
         call sorter('end',asort,asort,0,0,0,0,0,0,0,0,prnt)
c
c        ----- and do third quarter-transformation -----
c
         call iosys('rewind half on scratch',0,0,0,' ')
         maxij=0
         ij=0
         do 400 i=1,nactiv
            do 300 j=1,i
               ij=ij+1
c
c              ----- check that this triangle is in core -----
c
               if (ij.gt.maxij) then
                  minij=maxij+1
                  maxij=min(nnpact,maxij+ntriang)
                  lnread=(maxij-minij+1)*nnp
                  call iosys('read real half from scratch '//
     #                       'without rewinding',lnread,values,0,' ')
               end if
c
c              ----- third quarter-transformation -----
c
               call trtosq(t1,values(1,ij-minij+1),nbf,nnp)
               call ebc(third(1,ij),t1,c(1,ncore+1),nbf,nbf,nactiv)
c
  300       continue
  400    continue
c
c        ----- and write out these transformed integrals -----
c
         call iosys('write real der_i(xbcd) to ints without rewinding',
     #               nbf*nactiv*nnpact,third,0,' ')
c
c        ----- fourth quarter-transformation -----
c
         do 500 ij=1,nnpact
            call ebtc(t1,c(1,ncore+1),third(1,ij),nactiv,nbf,nactiv)
            call sqtotr(mo(1,ij),t1,nactiv,nnpact)
c..bhl
c            if(der.eq.1) then
c             write(iout,*)' ndf=1 der_g mo ij= ',ij
c             call matout(t1,nactiv,nactiv,nactiv,nactiv,iout)
c            endif
c..bhl
  500    continue
c
         call iosys('write real der_g to ints without rewinding',
     #               nnpact**2,mo,0,' ')
c
 1000 continue
c
      return
      end
