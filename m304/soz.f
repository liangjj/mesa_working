*deck @(#)soz.f	5.1  11/6/94
      subroutine soz(nprim,nprimi,nprimj,pointi,pointj,exp,xyz,
     #               itype,jtype,
     #               ncarti,ncartj,zso)
c
c***module to form the primitive one-electron spin-orbit integrals
c
c operator lz - integrals are done analytically
c
c matt braunstein                dec. 1991                  lanl
c
      implicit integer (a-z)
c
      real*8 exp(nprim),zso(xyz),facp,facd2,facd4,facd6,f1,f3,f8
c
      common /io/ inp,iout
c
      data facp/2.094395102/,
     $     facd2/0.418879021/,facd4/0.8377580410/,facd6/1.256637061/,
     $     f1/1.795195802/,f3/3.949430764/,f8/0.1196797201/
c
c g functions and above not yet programmed
c
      if(itype.ge.5) then
        call lnkerr ('m304: soint, angular momentum type not 
     #                yet programmed')
      end if
c
c for s type functions, the integral is zero
c
c   zero array
c
        do 5 i=1,xyz
          zso(i)=0.d0
5       continue
      if(itype.eq.1)then
        return
      end if
c
c d,f functions not yet tested
c
      if(itype.eq.3.or.itype.eq.4)then
        write(iout,*)'warning - d,f functions not yet programmed'
      end if
      
c
c p type functions
c
      if (itype.eq.2) then 
        icount=0
        do 10 i=1,ncarti
          do 20 j=1,ncartj
            do 30 k=1,nprimj
              do 40 l=1,nprimi
                icount=icount+1
                if (i.eq.1.and.j.eq.2) then
                  zso(icount)=-facp/
     #                        (exp(pointi+(l-1))+exp(pointj+(k-1)))
                end if
                if (i.eq.2.and.j.eq.1) then
                  zso(icount)=facp/
     #                        (exp(pointi+(l-1))+exp(pointj+(k-1)))
                end if
40            continue
30          continue
20        continue
10      continue
        return
      end if
c
c d type functions
c
      if (itype.eq.3) then
        icount=0
        do 50 i=1,ncarti
          do 60 j=1,ncartj
            do 70 k=1,nprimj
              do 80 l=1,nprimi
                icount=icount+1
                if (i.eq.1.and.j.eq.4) then
                  zso(icount)=-facd6/
     #                        (exp(pointi+(l-1))+exp(pointj+(k-1)))**2
                end if
                if (i.eq.2.and.j.eq.4) then
                  zso(icount)=facd6/
     #                        (exp(pointi+(l-1))+exp(pointj+(k-1)))**2
                end if
                if (i.eq.4.and.j.eq.1) then
                  zso(icount)=facd4/
     #                        (exp(pointi+(l-1))+exp(pointj+(k-1)))**2
                end if
                if (i.eq.4.and.j.eq.2) then
                  zso(icount)=-facd4/
     #                        (exp(pointi+(l-1))+exp(pointj+(k-1)))**2
                end if
                if (i.eq.5.and.j.eq.6) then
                  zso(icount)=-facd2/
     #                        (exp(pointi+(l-1))+exp(pointj+(k-1)))**2
                end if
                if (i.eq.6.and.j.eq.5) then
                  zso(icount)=facd2/
     #                        (exp(pointi+(l-1))+exp(pointj+(k-1)))**2
                end if
80           continue
70         continue
60       continue
50     continue
       return
      end if
c
c f type functions
c
      if (itype.eq.4) then
        icount=0
        do 90 i=1,ncarti
          do 100 j=1,ncartj
            do 110 k=1,nprimj
              do 120 l=1,nprimi
                icount=icount+1
                if (i.eq.1.and.j.eq.5) then
                  zso(icount)=-f1/
     #                        (exp(pointi+(l-1))+exp(pointj+(k-1)))**3
                end if
                if (i.eq.2.and.j.eq.4) then
                  zso(icount)=f1/
     #                        (exp(pointi+(l-1))+exp(pointj+(k-1)))**3
                end if
                if (i.eq.4.and.j.eq.2) then
                  zso(icount)=-f3/
     #                        (exp(pointi+(l-1))+exp(pointj+(k-1)))**3
                end if
                if (i.eq.4.and.j.eq.5) then
                  zso(icount)=(2.0*f3)/
     #                        (exp(pointi+(l-1))+exp(pointj+(k-1)))**3
                end if
                if (i.eq.5.and.j.eq.1) then
                  zso(icount)=(3.0*f3)/
     #                        (exp(pointi+(l-1))+exp(pointj+(k-1)))**3
                end if
                if (i.eq.5.and.j.eq.4) then
                  zso(icount)=(-2.0*f3)/
     #                        (exp(pointi+(l-1))+exp(pointj+(k-1)))**3
                end if
                if (i.eq.6.and.j.eq.10) then
                  zso(icount)=-f3/
     #                        (exp(pointi+(l-1))+exp(pointj+(k-1)))**3
                end if
                if (i.eq.7.and.j.eq.8) then
                  zso(icount)=-f3/
     #                        (exp(pointi+(l-1))+exp(pointj+(k-1)))**3
                end if
                if (i.eq.8.and.j.eq.7) then
                  zso(icount)=f3/
     #                        (exp(pointi+(l-1))+exp(pointj+(k-1)))**3
                end if
                if (i.eq.9.and.j.eq.10) then
                  zso(icount)=f3/
     #                        (exp(pointi+(l-1))+exp(pointj+(k-1)))**3
                end if
                if (i.eq.10.and.j.eq.6) then
                  zso(icount)=(2.0*f8)/
     #                        (exp(pointi+(l-1))+exp(pointj+(k-1)))**3
                end if
                if (i.eq.10.and.j.eq.9) then
                  zso(icount)=(-2.0*f8)/
     #                        (exp(pointi+(l-1))+exp(pointj+(k-1)))**3
                end if
120          continue
110        continue
100       continue
90     continue
       return
      end if
      end
