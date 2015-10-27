*deck @(#)nrmlz.f
c***begin prologue     nrmlz
c***date written       920417   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           nrmlz, link 7000
c***author             schneider, barry (nsf)
c***source             m7000
c***purpose            normalize functions
c
c***references       
c
c***routines called
c***end prologue       nrmlz
      subroutine nrmlz(fns,fnsc,ddfns,ddfnsc,type,nbf,npts)
      implicit integer (a-z)
      real *8 fns, ddfns, sdot, norm
      complex *16 fnsc, ddfnsc, cdotu, cnorm
      character *(*) type
      dimension fns(npts,nbf), fnsc(npts,nbf), ddfns(npts,nbf)
      dimension ddfnsc(npts,nbf)
      if (type.eq.'real') then
          do 10 i=1,nbf
             norm=sdot(npts,fns(1,i),1,fns(1,i),1)
             norm=1.d0/sqrt(norm)
             do 20 j=1,npts
                fns(j,i)=norm*fns(j,i)
                ddfns(j,i)=norm*ddfns(j,i)
   20        continue
   10     continue
      elseif (type.eq.'complex') then
          do 30 i=1,nbf
             cnorm=cdotu(npts,fnsc(1,i),1,fnsc(1,i),1)
             cnorm=1.d0/sqrt(cnorm)
             do 40 j=1,npts
                fnsc(j,i)=cnorm*fnsc(j,i)
                ddfnsc(j,i)=cnorm*ddfnsc(j,i)
   40        continue
   30     continue
      endif      
      return
      end







