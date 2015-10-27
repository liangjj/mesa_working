*deck @(#)project.f	5.1  11/6/94
      subroutine project(fkmtri,fkproj,projop,proj,coord,mass,
     $                   temp1,temp2,natoms,natoms3,nvv,projg)
c***begin prologue     project.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             page, michael(nrl)
c***source             @(#)project.f	5.1   11/6/94
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       project.f
      implicit none
c     --- input variables -----
      integer natoms,natoms3,nvv
      logical projg
c     --- input arrays (unmodified) ---
      real*8 coord(3,natoms),mass(natoms)
c     --- input arrays (scratch) ---
      real*8 temp1(natoms3,natoms3),temp2(natoms3,natoms3)
      real*8 fkmtri(nvv),proj(natoms3,natoms3),projop(natoms3,7)
c     --- output arrays ---
      real*8 fkproj(nvv)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer nproj,i,j,atom,ij,icoord
      logical debug
      real*8 zero,one,two
c
      parameter (zero=0.0d+00,one=1.0d+00,two=2.0d+00)
c
      common /io/ inp,iout
c
      if(projg) then
         nproj=7
      else
         nproj=6
      endif
c
      do 2 i=1,nproj
         do 1 j=1,natoms3
            projop(j,i)=zero
    1    continue
    2 continue
c
      do 4 atom=1,natoms
         do 3 i=1,3
            icoord=(atom-1)*3+i
            projop(icoord,i)=sqrt(mass(atom))
    3    continue
    4 continue
c
      do 5 atom=1,natoms
         projop(3*atom-1,4)=coord(3,atom)*sqrt(mass(atom))
         projop(3*atom  ,4)=-coord(2,atom)*sqrt(mass(atom))
    5 continue
c
      do 6 atom=1,natoms
         projop(3*atom-2,5)=coord(3,atom)*sqrt(mass(atom))
         projop(3*atom  ,5)=-coord(1,atom)*sqrt(mass(atom))
    6 continue
c
      do 7 atom=1,natoms
         projop(3*atom-2,6) =  coord(2,atom)*sqrt(mass(atom))
         projop(3*atom-1,6) = -coord(1,atom)*sqrt(mass(atom))
    7 continue
c
      if(projg) then
         call iosys('read real "cartesian first derivatives" from rwf',
     $               natoms3,projop(1,nproj),0,' ')
         do 8 atom=1,natoms
            projop(3*atom-2,7)=projop(3*atom-2,7)/sqrt(mass(atom))
            projop(3*atom-1,7)=projop(3*atom-1,7)/sqrt(mass(atom))
            projop(3*atom  ,7)=projop(3*atom  ,7)/sqrt(mass(atom))
    8    continue
      endif
c
      if(debug) then
         write(iout,1224)
 1224    format(//,' projection vectors ',/)
         call matout(projop,natoms3,7,natoms3,7,iout)
      endif
c
      call schmidt(projop,temp1,temp2,nproj,natoms3)
      if(debug) then
         write(iout,1294)
 1294    format(//,'schmidted  projection vectors ',/)
         call matout(projop,natoms3,6,natoms3,6,iout)
      endif
c
c
      call embct(proj,projop,projop,natoms3,nproj,natoms3)
      do 13 i=1,natoms3
         proj(i,i)=one+proj(i,i)
   13 continue
      if(debug) then
         write (iout,1225)
 1225    format(//,'   projection matrix',/)
         call matout(proj,natoms3,natoms3,natoms3,natoms3,iout)
      endif
c
      ij=0
      do 10 i=1,natoms3
         do 9 j=1,i
            ij=ij+1
            temp1(i,j)=fkmtri(ij)
            temp1(j,i)=fkmtri(ij)
    9    continue
   10 continue
c
      call ebc(temp2,temp1,proj,natoms3,natoms3,natoms3)
      call ebtc(temp1,proj,temp2,natoms3,natoms3,natoms3)
c
      ij=0
      do 12 i=1,natoms3
         do 11 j=1,i
            ij=ij+1
            fkproj(ij)=(temp1(i,j)+temp1(j,i))/two
   11    continue
   12 continue
c
      if(debug) then
         write (iout,1226)
 1226    format(//,'   projected force-constant matrix',/)
         call matout(temp1,natoms3,natoms3,natoms3,natoms3,iout)
      endif
c
c
      return
      end
