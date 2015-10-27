*deck @(#)vsarr.f	1.1 9/8/91
c***begin prologue     vsarr
c***date written       920417   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m7001, vectors
c***author             schneider, barry (nsf)
c***source             m7000
c***purpose            rearrange vectors
c***routines called    iosys, util and mdutil
c***end prologue       vsarr
      subroutine vsarr(fns,fnsc,ddfns,ddfnsc,work,nsym,lb,mb,nout,
     1                 nbf,n,dimsym,type)
      implicit integer (a-z)
      real *8  fns, ddfns, work
      complex *16 fnsc, ddfnsc
      character *(*) type
      dimension fns(n,nbf), ddfns(n,nbf), nsym(dimsym), work(nbf)
      dimension lb (nbf), mb(nbf)
      dimension fnsc(n,nbf), ddfnsc(n,nbf)
      ii=0
      nout=0
      do 10 i=1,dimsym
         if (nsym(i).ne.0) then
             newsym=0
             do 20 j=1,nsym(i)
                if (work(ii+j).ne.0.d0) then
                    newsym=newsym+1
                    nout=nout+1
                    lb(nout)=lb(ii+j)
                    mb(nout)=mb(ii+j)
                    if (type.eq.'real') then
                        call copy(fns(1,ii+j),fns(1,nout),n)
                        call copy(ddfns(1,ii+j),ddfns(1,nout),n)
                    elseif (type.eq.'complex') then
                        call cc2opy(fnsc(1,ii+j),fnsc(1,nout),n)
                        call cc2opy(ddfnsc(1,ii+j),ddfnsc(1,nout),n)
                    endif
                endif
   20        continue
             ii=ii+nsym(i)
             nsym(i)=newsym 
         endif
   10 continue   
      return
      end

