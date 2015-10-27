*deck @(#)sph.f	1.1 9/8/91
      subroutine sph(plm,cphi,sphi,cmp,cmphi,smphi,ylmp,ylmm,dfct,n,
     1               mu,lmax,maxfac,pass)
      implicit integer (a-z)
      common /io/ inp, iout
      dimension plm(n,0:lmax), ylmp(n,0:lmax), ylmm(n,0:lmax)
      dimension dfct(0:maxfac)
      dimension cphi(n), sphi(n), cmp(n), cmphi(n), smphi(n)
      real *8 plm, cphi, sphi, cmphi, smphi, const, ylmp, ylmm
      real *8 sq2, dfct, pi, fourpi
      complex*16 cmp
      data pi, fourpi /3.14159265358979323846d+00,
     1                12.56637061435917295384d+00/
      sq2=sqrt(2.d+00)
      if(mu.eq.0) then
         do 10 l=0,lmax
            const=sqrt((2*l+1)/fourpi)
            do 20 i=1,n
               ylmp(i,l)=plm(i,l)*const
               plm(i,l)=ylmp(i,l)
   20       continue
   10    continue
         words=n*(lmax+1)
         call iosys ('write real ylm to ylms without rewinding',words,
     1               ylmp,0,' ')
         write (iout,100) pass,mu,words
         return
      else
         do 30 i=1,n
            cmp(i)=(dcmplx(cphi(i),sphi(i)))**mu
            cmphi(i)=real(cmp(i))
            smphi(i)=imag(cmp(i))
   30    continue
         do 40 l=mu,lmax
         const=sqrt((2*l+1)/fourpi*dfct(l-mu)/dfct(l+mu))
c  sqrt(2) factor added to normalize "real valued" ylms
         const=const*sq2
            do 50 i=1,n
               ylmp(i,l)=plm(i,l)*cmphi(i)*const
               ylmm(i,l)=plm(i,l)*smphi(i)*const
   50       continue
   40    continue
         words=n*(lmax+1)
         call iosys ('write real ylm to ylms without rewinding',words,
     1                  ylmp(1,0),0,' ')
         call iosys ('write real ylm to ylms without rewinding',words,
     1                  ylmm(1,0),0,' ')
         words=words+words
         write (iout,100) pass,mu,words
      endif
  100 format(/,5x,'pass',1x,i3,2x,'m value',1x,i3,2x,'words written',1x,
     1             i8)
      return
      end
