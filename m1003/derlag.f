*deck @(#)derlag.f	5.1  11/6/94
      subroutine derlag(talag,scr,grad,buf,lbufso,mixhes,
     1                  nco,nao,nob,nmix,ldf)
c
c   add ta_der_lagrangian to mo_der_lagrangian
c
      implicit integer(a-z)
      real*8 talag(*),scr(*),grad(*),buf(*)
      integer mixhes(nob,*)
      logical debug
c
      parameter (debug=.false.)
c
      common/io/inp,iout
c
c
      noc=nco+nao
      nvo=nob-noc
      nocnob=noc*nob
c
      if (debug) then
         write(iout,*)' mixhes '
         do 1 i=1,nob
            write(iout,9001)(mixhes(i,j),j=1,noc)
    1    continue
      end if

 9001 format(1x,7(1x,i4))
c
c
      ipass=lbufso/nocnob
      ipass=min(ldf,ipass)
      if(ipass.lt.1) then
         call lnkerr(' m1001: derlag increase lbufso')
      end if
      npass=(ldf-1)/ipass+1
c
      call iosys('rewind mo_der_lagrangian on rwf',0,0,0,' ')
      kx=1
      ix=0
c
      do 10 i=1,ldf,ipass
         kpass=min(ipass,ldf-i+1)
         nread=kpass*nocnob
c
         call iosys('read real mo_der_lagrangian from rwf '//
     1              'without rewinding',nread,buf,0,' ')
c
         jx=0
         do 6 m=1,kpass
            do 5 j=1,nocnob
               scr(j)=talag(ix+j)+buf(jx+j)*.5
    5       continue
c
            if (debug) then
               write(iout,*)' mo_der_lagr'
               call matout(buf(jx+1),nob,noc,nob,noc,6)
            end if
c
            ix=ix+nocnob
            jx=jx+nocnob
c
            call mcgrd(grad(kx),scr,nob,nco,nao,nvo,lmix,lmixt,
     1                 mixhes,nob)
c
            if (debug) then
               write(iout,*)' rhs of orb. gradient'
               write(iout,9901)(grad(kx-1+mm),mm=1,nmix)
            end if
 9901 format(5(1x,f12.8))
c
            kx=kx+nmix
c
    6    continue
c
 10   continue
c
c
c
      return
      end
