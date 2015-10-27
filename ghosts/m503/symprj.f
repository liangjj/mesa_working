*deck %W%  %G%
      subroutine symprj(c,scr,ct,sympt,nbf,nsym,nobs)
      implicit integer(a-z)
      real*8  c(nbf,nbf),ct(nbf,nbf),scr(nbf,nbf)
      real*8 thrshv,thrshe,thrshd
      integer sympt(nbf),nobs(nsym)
c
      common /io/inp,iout
c
      data thrshv/1.d-4/,thrshe/1.d-3/,thrshd/1.d-5/
c
c         write(iout,*)'  SortBl ',thrshe
c
c
c
      call iosys('read integer sympt from rwf',nbf,sympt,
     #   0,' ')
cc
c         write(iout,*)' sympt '
c         write(iout,101)(sympt(i),i=1,nbf)
c  101    format(2x,5i5)
cc
cc
       iof=0
       do 5 is=1,nsym
       nob=nobs(is)
c
c   symmetry modulo 4 for d2h case
c
       iss=is-4*((is-1)/4)
c
       do 4 i=iof+1,iof+nob
        do 3 j=1,nbf
         if(sympt(j).ne.iss) then
           c(j,i)=0.0
         endif
  3     continue
  4    continue
       iof=iof+nob
  5    continue
c
c   d2h projection
c
      if(nsym.eq.8) then
c
       n2=nbf/2
c
       do 40 i=1,nbf
         do 30 j=1,nbf
         ct(j,i)=abs(c(j,i))
  30     continue
  40   continue
       do 41 i=1,nbf
        do 31 j=1,n2
         scr(j,i)=(ct(j,i)+ct(j+n2,i))*.5
  31    continue
  41   continue
       do 42 i=1,nbf
        do 32 j=1,n2
         c(j,i)=sign(scr(j,i),c(j,i))
  32    continue
  42   continue
       do 43 i=1,nbf
        do 33 j=1,n2
         c(j+n2,i)=sign(scr(j,i),c(j+n2,i))
  33    continue
  43   continue
      end if
c
c
       return
       end
