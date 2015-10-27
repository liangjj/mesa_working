*deck @(#)project.f	5.1  11/6/94
      subroutine project(fkmtri,fkproj,projop,proj,coord,mass,
     #                   temp1,temp2,natoms,natoms3,nvv)
c***begin prologue     project
c***date written       850601  
c***revision date      11/6/94      
c
c***keywords           
c***author             page, michael 
c***source             @(#)project.f	5.1   11/6/94
c***purpose            to project translations/rotations from the force constant
c                      matrix. 
c***description
c     
c    
c
c***references
c
c***routines called
c     %m%
c       start here
c
c***end prologue       project
c
      implicit integer (a-z)
c
      logical debug
c
      real*8 fkmtri(nvv),fkproj(nvv),projop(natoms3,6)
      real*8 proj(natoms3,natoms3),coord(3,natoms),mass(natoms)
      real*8 temp1(natoms3,natoms3),temp2(natoms3,natoms3)
c
      common /io/ inp,iout
c
      parameter (debug=.false.)
c
      do 2 i=1,6
         do 1 j=1,natoms3
            projop(j,i)=0.0d+00
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
      nproj=6
c
c     ----- check for linear system -----
c
      do 31 atom=1,natoms
         if (coord(1,atom).ne.0.0d+00.or.coord(2,atom).ne.0.0d+00)
     #   go to 32
   31 continue
      nproj=5
   32 continue
c
      if(debug) then
         write(iout,1224)
 1224    format(//,' projection vectors ',/)
         call matout(projop,natoms3,6,natoms3,6,6)
      end if
c
      call schmidt(projop,temp1,temp2,nproj,natoms3)
c
      if(debug) then
         write(iout,1294)
 1294    format(//,'schmidted  projection vectors ',/)
         call matout(projop,natoms3,6,natoms3,6,6)
      end if
c
c
      call embct(proj,projop,projop,natoms3,nproj,natoms3)
      do 8 i=1,natoms3
         proj(i,i)=1.0d+00+proj(i,i)
    8 continue
c
      if(debug) then
         write (iout,1225)
 1225    format(//,'   projection matrix',/)
         call matout(proj,natoms3,natoms3,natoms3,natoms3,6)
      end if
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
      if(debug) then
         write (iout,1227)
 1227    format(//,'  temp1 before ebc',/)
         call matout(temp1,natoms3,natoms3,natoms3,natoms3,6)
         write (iout,1287)
 1287    format(//,'  proj before ebc',/)
         call matout(proj,natoms3,natoms3,natoms3,natoms3,6)
      end if
c
c
      call ebc(temp2,temp1,proj,natoms3,natoms3,natoms3)
c
      if(debug) then
         write (iout,1297)
 1297    format(//,'  temp2 after ebc',/)
         call matout(temp2,natoms3,natoms3,natoms3,natoms3,6)
      end if
c
      call ebtc(temp1,proj,temp2,natoms3,natoms3,natoms3)
c
      ij=0
      do 12 i=1,natoms3
         do 11 j=1,i
            ij=ij+1
            fkproj(ij)=(temp1(i,j)+temp1(j,i))/2.0d+00
   11    continue
   12 continue
c
      if(debug) then
         write (iout,1226)
 1226    format(//,'   projected force-constant matrix',/)
         call matout(temp1,natoms3,natoms3,natoms3,natoms3,6)
      end if
c
      return
      end
