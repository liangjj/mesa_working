*deck @(#)d1xc.f	5.1 11/6/94
      subroutine d1xc(nat,nbf,nnp,nbtype,iatom,jatom,itype,jtype,
     $                nconti,ncontj,nobf,start,d,dkay,grad)
c***begin prologue     d1xc.f
c***date written       940513  
c***revision date      11/6/94      
c
c***keywords           
c***author             martin, richard(lanl) 
c***source             @(#)d1xc.f	5.1   11/6/94
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       d1xc.f
      implicit none
c     --- input variables -----
      integer nat,nbf,nnp,nbtype
      integer iatom,jatom,itype,jtype,nconti,ncontj
c     --- input arrays (unmodified) ---
      integer start(nat,nbtype),nobf(nbtype)
      real*8 d(nnp),dkay(nbf,nbf,3)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 grad(3,nat)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer numi,numj,istart,jstart,if,jf,ic,jc,i,j,ij,coord
c
c
      do 100 coord=1,3
c
c        --- put this angular momentum block in the right places
         numi=nobf(itype)
         numj=nobf(jtype)
         istart=start(iatom,itype)
         jstart=start(jatom,jtype)
c
         do 4 if=1,numi
            do 3 jf=1,numj
               do 2 jc=1,ncontj
                  j=jstart+(jc-1)*numj+jf
                  do 1 ic=1,nconti
                     i=istart+(ic-1)*numi+if
                     if (i.ge.j) then
                        ij=i*(i-1)/2+j
                     else
                        ij=j*(j-1)/2+i
                     end if
c
                     grad(coord,iatom)=grad(coord,iatom)+
     $                      dkay(i,j,coord)*d(ij)
                     grad(coord,jatom)=grad(coord,jatom)-
     $                      dkay(i,j,coord)*d(ij)
    1             continue
    2          continue
    3       continue
    4    continue
c
  100 continue
c
c
      return
      end
