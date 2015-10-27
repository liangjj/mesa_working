*deck @(#)dfock.f	5.1  11/6/94
      subroutine dfock(values,corden,cfock,nnp,nbf,
     #                 fsq,isq,ntriang,nder)
c
      implicit integer (a-z)
c
      real*8 values(nnp,ntriang),corden(nnp),cfock(nnp)
      real*8 fsq(nbf,nbf),isq(nbf,nbf)
      real*8 dkl,djl,djk
c
c     ----- initialize the integral file, etc. ------
c
      call iosys('rewind "sorted ao derivative integrals"'//
     #           ' on scr',0,0,0,' ')
c
c     ----- loop over degrees of freedom -----
c
         call rzero(fsq,nbf**2)
         maxkl=0
         kl=0
c
         do 6 k=1,nbf
            do 5 l=1,k
               kl=kl+1
c
c              ----- check that this triangle is in core -----
c
               if (kl.gt.maxkl) then
                  minkl=maxkl+1
                  maxkl=min(nnp,maxkl+ntriang)
                  lnread=(maxkl-minkl+1)*nnp
c
      call iosys('read real "sorted ao derivative integrals"'//
     #' from scr without rewinding',lnread,values,0,' ')
               end if
c
c              ----- coulomb term -----
c
               dkl=corden(kl)*2.0d+00
               if (k.eq.l) then
                  do 1 iq=1,nnp
                     values(iq,kl-minkl+1)=values(iq,kl-minkl+1)*0.5d+00
    1             continue
               end if
c
               call saxpy(nnp,dkl,values(1,kl-minkl+1),1,cfock,1)
c
c              ----- exchange term -----
c
               call trtosq(isq,values(1,kl-minkl+1),nbf,nnp)
c
               do 4 j=1,nbf
                  junk=max(j,l)
                  djl=corden(junk*(junk-1)/2+min(j,l))
                  junk=max(j,k)
                  djk=corden(junk*(junk-1)/2+min(j,k))
                  do 2 i=k,nbf
                     fsq(i,k)=fsq(i,k)+isq(i,j)*djl
    2             continue
                  do 3 i=l,nbf
                     fsq(i,l)=fsq(i,l)+isq(i,j)*djk
    3             continue
    4          continue
    5       continue
    6    continue
c
c        ----- combine coulomb and exchange to form fock-matrix -----
c
         ij=0
         do 8 i=1,nbf
            do 7 j=1,i
               ij=ij+1
               cfock(ij)=cfock(ij)-fsq(i,j)*0.5d+00
    7       continue
    8    continue
c
c
c
      return
      end
