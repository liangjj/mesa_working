*deck @(#)scalfn.f
c***begin prologue     scalfn
c***date written       920417   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           scalfn, link 7000
c***author             schneider, barry (nsf)
c***source             m7000
c***purpose            scale functions
c
c***references       
c
c***routines called
c***end prologue       scalfn
      subroutine scalfn(fns,fnsc,vec,type,nbf,npts)
      implicit integer (a-z)
      real *8 fns, vec
      complex *16 fnsc
      character *(*) type
      dimension fns(npts,nbf), fnsc(npts,nbf), vec(npts)
      if (type.eq.'real') then
          do 10 i=1,nbf
              call vmul(fns(1,i),fns(1,i),vec,npts)        
   10     continue
      elseif (type.eq.'complex') then
          do 20 i=1,nbf
              call cvmul(fnsc(1,i),fnsc(1,i),vec,vec,npts,'real')        
   20     continue
      endif      
      return
      end







