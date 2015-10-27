*deck @(#)vearr.f	1.1 9/8/91
c***begin prologue     vearr
c***date written       920417   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m7001, vectors
c***author             schneider, barry (nsf)
c***source             m7000
c***purpose            rearrange vectors
c***routines called    iosys, util and mdutil
c***end prologue       vearr
      subroutine vearr(s,sc,eig,eigc,tol,nsym,lb,mb,nout,nbf,
     1                 dimsym,type)
      implicit integer (a-z)
      real *8  s, eig, eg, tol, teig
      complex *16 sc, eigc, cdotu, tmp
      character *(*) type
      dimension s(nbf,nbf), nsym(dimsym), eig(nbf), lb (nbf), mb(nbf)
      dimension sc(nbf,nbf), eigc(nbf)
      ii=0
      nout=0
      do 10 i=1,dimsym
         if (nsym(i).ne.0) then
             newsym=0
             do 20 j=1,nsym(i)
                if (type.eq.'real') then
                    teig=abs(eig(ii+j))
                elseif (type.eq.'complex') then
                    teig=abs(eigc(ii+j))
                endif
                if (teig.gt.tol) then                    
                    newsym=newsym+1
                    nout=nout+1
                    lb(nout)=lb(ii+j)
                    mb(nout)=mb(ii+j)
                    if (type.eq.'real') then
                        eg=1.d0/sqrt(eig(ii+j)) 
                        call smul(s(1,nout),s(1,ii+j),eg,nbf)
                    elseif (type.eq.'complex') then
                        tmp=cdotu(nbf,sc(1,ii+j),1,sc(1,ii+j),1)
                        tmp=1.d0/sqrt(tmp*eigc(ii+j))     
                        call csmul(sc(1,nout),sc(1,ii+j),tmp,nbf)
                    endif
                endif
   20        continue
             ii=ii+nsym(i)
             nsym(i)=newsym 
         endif
   10 continue
      return
      end
