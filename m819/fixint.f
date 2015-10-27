*deck  @(#)fixint.f	5.1 11/6/94
      subroutine fixint(s,t,v,values,nnp,nbf,nsmall,ntriang,ops)
      implicit integer(a-z)
      character*(*) ops
      real*8 s(nnp),t(nnp),v(nnp),values(*)
      real*8 ezero,epsilon,fpkey 
      logical logkey
      common /io/ inp,iout 
      character*128 tints,zints 
c
      call iosys('read character "transformed integral filename" '//
     $           'from rwf',0,0,0,tints)
      call iosys('read character "zeroed integral filename" '//
     $           'from rwf',0,0,0,zints)
c 
      len=nnp**2+nnp+200000
      len=min(len,5000000)
      if(logkey(ops,'unit=ssd=zint',.false.,' ')) then
         call iosys('open zints as new on ssd',len,0,0,zints)
      else
         call iosys('open zints as new',len,0,0,zints)
      end if
c
c
      call iosys('open tints as old',0,0,0,tints)
c
c     ----- read in s, t and v one-electron integrals -----
c
      call iosys('read real "mo one-electron integrals" from tints',
     $            nnp,v,0,' ')
c
      ezero=0.0d+00
      epsilon=fpkey(ops,'effective-integral-zero', 0.0d+00,' ')
      iq=nsmall*(nsmall+1)/2+1
      do 1 i=nsmall+1,nbf
         call rzero(v(iq),i-1)
         iq=iq+i
         v(iq-1)=ezero
         ezero=ezero+epsilon  
1     continue
c      
c
      call iosys('write real "mo one-electron integrals" to zints',
     $            nnp,v,0,' ')
c
      call iosys('create real file "mo two-electron integrals"'//
     $           ' on zints',-1,0,0,' ')
c
      call iosys('rewind "mo two-electron integrals" on tints',
     $            0,0,0,' ')
c
c
c     ----- read through integrals if first iteration, or if
c           cannot hold all the integrals in core -----
c
      maxkl=0
c
c -- read the first block of integrals
c
      ip=1
      minkl=maxkl+1
      maxkl=min(nnp,maxkl+ntriang)
      lnread=(maxkl-minkl+1)*nnp
      call iosys('read real "mo two-electron integrals" from tints '//
     $           'without rewinding',lnread,values,0,' ')
      lastrd=lnread
c
      ioff=nsmall*(nsmall+1)/2
      lzero=nnp-ioff
c
      kl=0
      do 6 k=1,nsmall
         do 5 l=1,k
            kl=kl+1
c
c           ----- check that this triangle of integrals is in core -----
c
            if (kl.gt.maxkl) then
                ip=1
                minkl=maxkl+1
                maxkl=min(nnp,maxkl+ntriang)
                lnread=(maxkl-minkl+1)*nnp
                call iosys('write real "mo two-electron integrals" to'//
     $                     ' zints without rewinding',lastrd,values,
     $                      0,' ')
                call iosys('read real "mo two-electron integrals" '//
     $                     'from tints without rewinding',lnread,
     $                      values,0,' ')
                lastrd=lnread
            end if
            call rzero(values(ip+ioff),lzero)
            ip=ip+nnp
    5    continue  
    6 continue
c
      do 8 k=nsmall+1,nbf
         do 7 l=1,k
            kl=kl+1
c
c           ----- check that this triangle of integrals is in core -----
c
            if (kl.gt.maxkl) then
                ip=1
                minkl=maxkl+1
                maxkl=min(nnp,maxkl+ntriang)
                lnread=(maxkl-minkl+1)*nnp
                call iosys('write real "mo two-electron integrals" to'//
     $                      ' zints without rewinding',lastrd,values,
     $                      0,' ')
                call iosys('read real "mo two-electron integrals" '//
     $                     'from tints without rewinding',lnread,values,
     $                     0,' ')
                lastrd=lnread
            end if
            call rzero(values(ip),nnp)
            ip=ip+nnp
    7    continue  
    8 continue
c
      call iosys('write real "mo two-electron integrals" to zints '//
     $            'without rewinding',lastrd,values,0,' ')
c
      call iosys('endfile "mo two-electron integrals" on zints',
     $            0,0,0,' ')
c
c
      return
      end
