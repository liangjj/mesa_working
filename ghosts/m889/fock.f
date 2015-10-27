*deck @(#)fock.f	1.1  11/30/90
      subroutine fock(values,d,f,nnp,num,fsq,isq,ntriang)
c
c***purpose:to form the two-electron contribution to the fock matrix
c           for closed-shell rhf using an ordered integral list.
c
c paul saxe                  27 july 1984                       lanl
c
      implicit integer (a-z)
c
      real*8 values(nnp,ntriang),d(nnp),f(nnp),fsq(num,num)
      real*8 isq(num,num),djl,djk,dkl,dil,dik
c
      data call /0/
c
      save start
c
c
      call=call+1
      call iosys('rewind "sorted ao integrals" on ints',0,0,0,' ')
      call rzero(f,nnp)
      call rzero(fsq,num**2)
c
      maxkl=0
      kl=0
      do 6 k=1,num
         do 5 l=1,k
            kl=kl+1
c
c     ----- check that this triangle of integrals is in core -----
c
            if (kl.gt.maxkl) then
               minkl=maxkl+1
               maxkl=min(nnp,maxkl+ntriang)
               lnread=(maxkl-minkl+1)*nnp
               call iosys('read real "sorted ao integrals" from ints '
     #                    //'without rewinding',lnread,values,0,' ')
            end if
c
c    ----- form coulomb term -----
c
            dkl=d(kl)*2.0d+00
            if (k.eq.l) then
               do 1 iq=1,nnp
                  values(iq,kl-minkl+1)=0.5d+00*values(iq,kl-minkl+1)
    1          continue
            end if
            call saxpy(nnp,dkl,values(1,kl-minkl+1),1,f,1)
c
c    ----- form exchange term -----
c
            call trtosq(isq,values(1,kl-minkl+1),num,nnp)
c
            do 4 j=1,num
               junk=max(j,l)
               djl=d(junk*(junk-1)/2+min(j,l))
               junk=max(j,k)
               djk=d(junk*(junk-1)/2+min(j,k))
               do 2 i=k,num
                  fsq(i,k)=fsq(i,k)+isq(i,j)*djl
    2          continue
               do 3 i=l,num
                  fsq(i,l)=fsq(i,l)+isq(i,j)*djk
    3          continue
    4       continue
    5    continue
    6 continue
c
c     ----- combine coulomb and exchange to form fock-matrix -----
c
      ij=0
      do 8 i=1,num
         do 7 j=1,i
            ij=ij+1
            f(ij)=f(ij)-fsq(i,j)*0.5d+00
    7    continue
    8 continue
c
c
      return
      end
