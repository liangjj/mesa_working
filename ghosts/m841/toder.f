*deck @(#)toder.f	1.1  11/30/90
      subroutine toder(dm,val,lab,bin,rsort,isort,lnsort,
     $     bftgrp,bftcmp,grpsiz,gptij,gptkl,gklsiz,nbf,nnp,ngrp,
     $     nnpgrp,ntriang,n2pdm,ci,mcscf)
c
      implicit integer (a-z)
c
      real*8 dm(nnp,ntriang)
      real*8 val(nnp)
      integer lab(nnp)
      integer bin(nnp)
      integer isort(lnsort)
      real*8 rsort(*)
      integer bftgrp(nbf)
      integer bftcmp(nbf)
      integer grpsiz(ngrp)
      integer gptij(nnpgrp)
      integer gptkl(nnpgrp)
      integer gklsiz(nnpgrp)
      logical ci
      logical mcscf
      character*32 file
c
      common /io/ inp,iout
c
      ioff(i,j)=i*(i-1)/2+j
c
c      if (mcscf) then
c.m840         file='mcscf ao 2pdm'
c      else
c.m840         file='ci ao 2pdm'
c      end if
c
       file='sorted ao integrals'
c
c     ----- initalize the sorting routines -----
c
      call sorter('start',isort,rsort,lnsort,n2pdm,0,0,0,0,
     #             'group ordered 2e-ints','rwf',.false.)
c
c     ----- loop over buffer loads of ao 2e-ints -----
c
      call iosys('rewind "'//file//'" on ints',0,0,0,' ')
      maxij=0
      ij=0
      do 4 i=1,nbf
         igrp=bftgrp(i)
         igrpsz=grpsiz(igrp)
         icmp=bftcmp(i)
         do 3 j=1,i
            jgrp=bftgrp(j)
            jgrpsz=grpsiz(jgrp)
            jcmp=bftcmp(j)
            nnij=ioff(igrp,jgrp)
            ij=ij+1
c
c     ----- check that this triangle of integrals is in core -----
c
            if (ij.gt.maxij) then
               minij=maxij+1
               maxij=min(nnp,maxij+ntriang)
               lnread=(maxij-minij+1)*nnp
               call iosys('read real "'//file//'" from ints '
     #                    //'without rewinding',lnread,dm,0,' ')
            end if
c
            kl=0
            do 2 k=1,i
               kgrp=bftgrp(k)
               kgrpsz=grpsiz(kgrp)
               kcmp=bftcmp(k)
               if (i.eq.k) then
                  lmax=j
               else
                  lmax=k
               end if
               do 1 l=1,lmax
                  kl=kl+1
                  lgrp=bftgrp(l)
                  lgrpsz=grpsiz(lgrp)
                  lcmp=bftcmp(l)
                  nnkl=ioff(kgrp,lgrp)
c
                  val(kl)=dm(kl,ij-minij+1)
c
                  if (nnij.ge.nnkl) then
                     lab(kl)=gptij(nnij)+gptkl(nnkl)*gklsiz(nnij)+
     $                    icmp+igrpsz*(jcmp-1+jgrpsz*(kcmp-1+kgrpsz*
     $                    (lcmp-1)))
                  else
                     lab(kl)=gptij(nnkl)+gptkl(nnij)*gklsiz(nnkl)+
     $                    kcmp+kgrpsz*(lcmp-1+lgrpsz*(icmp-1+igrpsz*
     $                    (jcmp-1)))
                  end if
c
c
 1             continue
 2          continue
c
            call sorter('with bin',isort,rsort,0,kl,lab,bin,
     #                   val,0,0,0,.false.)
    3    continue
    4 continue
c
      call sorter('end',isort,rsort,0,0,0,0,0,0,0,0,.false.)
c
c
      return
      end
