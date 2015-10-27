c $Header: trnsf.f,v 1.2 93/01/01 12:43:28 bis Exp $
*deck @(#)trnsf.f	1.1 9/9/91
      subroutine trnsf(rpt,thpt,phpt,wtr,wtth,wtph,f,volume,yukawa,
     1                 intgrd,ns,nrmax,nthmax,nphmax,which)
      parameter (numshl=50 , nocen=10 )
      implicit real *8 (a-h,o-z)
      dimension rpt(nrmax,ns), thpt(nthmax,ns), phpt(nphmax,ns)
      dimension wtr(nrmax,ns), wtth(nthmax,ns), wtph(nphmax,ns)
      dimension f(*), dumcen(3,nocen)
      character*(*) intgrd, which
      common/spheri/ nshell, nr(numshl), nthell(numshl), nphell(numshl),
     1               ntheta(numshl,numshl), nphi(numshl,numshl)
      common /centsi/ncent
      common /quadpt/ nsub(numshl), nth(numshl), nph(numshl)
      common /centsr/cc(3),a(3,nocen),eta(nocen),
     1               znuc(nocen)
      common/io/ inp,iout
      if (which.eq.'all') then
          call copy(a,dumcen,3*nocen)
      else
          call rzero(dumcen,3*nocen)
      endif
      volume = 0.0d+00
      yukawa = 0.d0
      do 10 is=1,ns
         do 20 i=1,nr(is)
            do 30 j=1,nth(is)
               do 40 k=1,nph(is)
                  volume=volume+wtr(i,is)*wtth(j,is)*wtph(k,is)
   40          continue
   30       continue
   20    continue
   10 continue     
c
c
c  example is a yukawa potential at location a(1 to 3,i))
c
c
c
c  evaluate integrand in untransformed coordinates
c
      if (which.eq.'all') then
          call iosys ('create real "multicenter yukawa potential" on '//
     1                'lamdat',ns*nrmax*nthmax*nphmax,0,0,' ')
      else
          call iosys ('create real "yukawa potential on '//
     1                'center-'//which//'" on lamdat',
     2                 ns*nrmax*nthmax*nphmax,0,0,' ')
      endif
      do 50 is=1,ns
         call fillf(f,rpt(1,is),thpt(1,is),phpt(1,is),wtr(1,is),
     1              wtth(1,is),wtph(1,is),dumcen,yukawa,eta,nr(is),
     2              nth(is),nph(is),nrmax,nthmax,nphmax,ncent,nocen,
     3              intgrd,which)
         if (which.eq.'all') then
             call iosys ('write real "multicenter yukawa potential" '//
     1                   'to lamdat without rewinding',
     2                    nr(is)*nth(is)*nph(is),f,0,' ')
         else
             call iosys ('write real "yukawa potential on '//
     1                   'center-'//which//'" to lamdat without '//
     2                   'rewinding',nr(is)*nth(is)*nph(is),f,0,' ')
         endif
   50 continue
      return
      end


