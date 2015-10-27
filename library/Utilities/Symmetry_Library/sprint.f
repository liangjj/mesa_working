*deck @(#)sprint.f	5.1  11/6/94
      subroutine sprint(genpt,lambda,lirrep,gen,ngen,nirrep,lengen,
     #     gamma,nop,char,labop,t,bftran,ptbftr,maxmom,
     #     atprmt,natoms,ncart)
c
      implicit integer (a-z)
c
      integer genpt(nirrep),lambda(nirrep)
      integer atprmt(natoms,nop)
      character*(*) lirrep(nirrep),labop(nop)
      real*8 gen(lengen,ngen),gamma(lengen,nop),char(nirrep,nop)
      real*8 t(3,3,nop)
      real*8 bftran(*)
      integer ncart(0:maxmom)
      integer ptbftr(0:maxmom)
      common/io/inp,iout
c
      do 10 n=1,ngen
         write (iout,2) n
    2    format (//,t30,'generator #',i1)
         do 9 irrep=1,nirrep
            pt=genpt(irrep)
            if (lambda(irrep).eq.1) then
               write (iout,3) lirrep(irrep),gen(pt,n)
    3          format (/,t10,a5,t30,f10.6)
            else if (lambda(irrep).eq.2) then
               write (iout,4) lirrep(irrep),
     #              (gen(pt+i,n),i=0,3)
    4          format (/,t10,a5,t25,2f10.6,/,t25,2f10.6)
            else if (lambda(irrep).eq.3) then
               write (iout,5) lirrep(irrep),
     #              (gen(pt+i,n),i=0,8)
    5          format (/,t10,a5,t20,3f10.6,2(/,t20,3f10.6))
            end if
    9    continue
 10   continue
c
      do 20 n=1,nop
         write (iout,12) n
 12      format (//,t30,'operation #',i3)
         do 19 irrep=1,nirrep
            pt=genpt(irrep)
            if (lambda(irrep).eq.1) then
               write (iout,3) lirrep(irrep),gamma(pt,n)
            else if (lambda(irrep).eq.2) then
               write (iout,4) lirrep(irrep),
     #              (gamma(pt+i,n),i=0,3)
            else if (lambda(irrep).eq.3) then
               write (iout,5) lirrep(irrep),
     #              (gamma(pt+i,n),i=0,8)
            end if
 19      continue
 20   continue
c
c     ----- and print the character table -----
c
      mx=0
 30   continue
      mn=mx+1
      mx=min(mx+8,nop)
      write (iout,21) (labop(i),i=mn,mx)
 21   format(/,' character table',/,t15,8a12)
      do 25 irrep=1,nirrep
         write (iout,22) lirrep(irrep),(char(irrep,i),i=mn,mx)
 22      format (t10,a5,8f12.6)
 25   continue
      if (mx.lt.nop) go to 30
c
      do 40 op=1,nop
         write (iout,31) op
 31      format (/,t20,'operator #',i3,
     $          'coordinate transformation matrix',/)
         call matout(t(1,1,op),3,3,3,3,iout)
 40   continue
c
      do 50 op=1,nop
         write (iout,41) op
 41      format (/,t20,'operation #',i3,' s,p,d,... transformations'/)
         do 49 angmom=0,maxmom
            pt=ptbftr(angmom)+ncart(angmom)**2*(op-1)
            call matout(bftran(pt),ncart(angmom),ncart(angmom),
     $           ncart(angmom),ncart(angmom),iout)
 49      continue
 50   continue
c
      write (iout,51)
 51   format (//,t20,'atom permutation table (ops -->, atoms v)',/)
      do 53 atom=1,natoms
         write (iout,52) (atprmt(atom,i),i=1,nop)
 52      format (5x,25i3)
 53   continue
c
c
      return
      end
