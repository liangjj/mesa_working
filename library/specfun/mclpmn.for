        program mclpmn
c
c       ============================================================
c       purpose: this program computes the associated legendre 
c                functions pmn(z) and their derivatives pmn'(z) for
c                a complex argument using subroutine clpmn
c       input :  x --- real part of z
c                y --- imaginary part of z
c                m --- order of pmn(z),  m = 0,1,2,...,n
c                n --- degree of pmn(z), n = 0,1,2,...,n
c       output:  cpm(m,n) --- pmn(z)
c                cpd(m,n) --- pmn'(z)
c       examples:
c                n = 5, x = 0.5, y = 0.2
c
c       m     re[pmn(z)]    im[pmn(z)]    re[pmn'(z)]   im[pmn'(z)]
c      -------------------------------------------------------------
c       0    .252594d+00  -.530293d+00   -.347606d+01  -.194250d+01
c       1    .333071d+01   .135206d+01    .117643d+02  -.144329d+02
c       2   -.102769d+02   .125759d+02    .765713d+02   .598500d+02
c       3   -.572879d+02  -.522744d+02   -.343414d+03   .147389d+03
c       4    .335711d+03  -.389151d+02   -.226328d+03  -.737100d+03
c       5   -.461125d+03   .329122d+03    .187180d+04   .160494d+02
c
c                n = 5, x = 2.5, y = 1.0
c
c       m     re[pmn(z)]    im[pmn(z)]    re[pmn'(z)]   im[pmn'(z)]
c      -------------------------------------------------------------
c       0   -.429395d+03   .900336d+03   -.350391d+02   .193594d+04
c       1   -.216303d+04   .446358d+04   -.208935d+03   .964685d+04
c       2   -.883477d+04   .174005d+05   -.123703d+04   .381938d+05
c       3   -.273211d+05   .499684d+05   -.568080d+04   .112614d+06
c       4   -.565523d+05   .938503d+05   -.167147d+05   .219713d+06
c       5   -.584268d+05   .863328d+05   -.233002d+05   .212595d+06
c       ============================================================
c
        implicit double precision (x,y)
        implicit complex*16 (c,z)
        dimension cpm(0:40,0:40),cpd(0:40,0:40)
        write(*,*)'  please enter m, n, x and y '
        read(*,*) m,n,x,y
        write(*,30) m,n,x,y
        call clpmn(40,m,n,x,y,cpm,cpd)
        write(*,*)'   m   n    re[pmn(z)]    im[pmn(z)]    ',
     &            're[pmn''(z)]   im[pmn''(z)]'
        write(*,*)' -----------------------------------',
     &            '-------------------------------'
        do 10 j=0,n
10         write(*,20)m,j,cpm(m,j),cpd(m,j)
20      format(1x,2i4,1x,2d14.6,1x,2d14.6)
30      format(1x,'m =',i2,', ','n =',i2,', ','x =',f5.1,
     &         ', ','y =',f5.1)
        end


        subroutine clpmn(mm,m,n,x,y,cpm,cpd)
c
c       =========================================================
c       purpose: compute the associated legendre functions pmn(z)   
c                and their derivatives pmn'(z) for a complex 
c                argument
c       input :  x  --- real part of z
c                y  --- imaginary part of z
c                m  --- order of pmn(z),  m = 0,1,2,...,n
c                n  --- degree of pmn(z), n = 0,1,2,...,n
c                mm --- physical dimension of cpm and cpd
c       output:  cpm(m,n) --- pmn(z)
c                cpd(m,n) --- pmn'(z)
c       =========================================================
c
        implicit double precision (x,y)
        implicit complex*16 (c,z)
        dimension cpm(0:mm,0:n),cpd(0:mm,0:n)
        z=cmplx(x,y)
        do 10 i=0,n
        do 10 j=0,m
           cpm(j,i)=(0.0d0,0.0d0)
10         cpd(j,i)=(0.0d0,0.0d0)
        cpm(0,0)=(1.0d0,0.0d0)
        if (dabs(x).eq.1.0d0.and.y.eq.0.0d0) then
           do 15 i=1,n
              cpm(0,i)=x**i
15            cpd(0,i)=0.5d0*i*(i+1)*x**(i+1)
           do 20 j=1,n
           do 20 i=1,m
              if (i.eq.1) then
                 cpd(i,j)=(1.0d+300,0.0d0)
              else if (i.eq.2) then
                 cpd(i,j)=-0.25d0*(j+2)*(j+1)*j*(j-1)*x**(j+1)
              endif
20         continue
           return
        endif
        ls=1
        if (cdabs(z).gt.1.0d0) ls=-1
        zq=cdsqrt(ls*(1.0d0-z*z))
        zs=ls*(1.0d0-z*z)
        do 25 i=1,m
25         cpm(i,i)=-ls*(2.0d0*i-1.0d0)*zq*cpm(i-1,i-1)
        do 30 i=0,m
30         cpm(i,i+1)=(2.0d0*i+1.0d0)*z*cpm(i,i)
        do 35 i=0,m
        do 35 j=i+2,n
           cpm(i,j)=((2.0d0*j-1.0d0)*z*cpm(i,j-1)-(i+j-
     &              1.0d0)*cpm(i,j-2))/(j-i)
35      continue
        cpd(0,0)=(0.0d0,0.0d0)
        do 40 j=1,n
40         cpd(0,j)=ls*j*(cpm(0,j-1)-z*cpm(0,j))/zs
        do 45 i=1,m
        do 45 j=i,n
           cpd(i,j)=ls*i*z*cpm(i,j)/zs+(j+i)*(j-i+1.0d0)
     &              /zq*cpm(i-1,j)
45      continue
        return
        end
