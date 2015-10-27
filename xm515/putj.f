*deck @(#)putj.f	1.1  4/25/95
      subroutine putj(jmat,nnprim,ndmat,jij,ni,nj,nfi,nfj,is,js,
     $                ijsh)
c***begin prologue     putj.f
c***date written       870702   (yymmdd)  
c***revision date      11/6/94      
c
c   december 8, 1993   rlm at lanl
c     modifying get1dm for use in direct scf code.
c***keywords           coulomb matrices, direct
c***author             saxe, paul (lanl)
c***source             @(#)putj.f	1.1   4/25/95
c***purpose            scatters a shell block of integrals into
c                      the full matrix.
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       putj.f
      implicit none
c     --- input variables -----
      integer nnprim,ndmat
      integer ni,nj
      integer nfi,nfj
      integer is,js
      logical ijsh
c     --- input arrays (unmodified) ---
      real*8 jij(ni,nj,nfi,nfj,ndmat)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 jmat(nnprim,ndmat)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer i,j,if,jf,jpos,ipos
      integer ii,ij,jj,dmat
      integer iftop,itop
      integer inp,iout
      logical debug
c
      parameter (debug=.false.)
      common/io/inp,iout
c
c     -- scatter the local ij portion into the full density ---
c        we are accumulating the local jij into the lower triangle of
c        the full j-matrix. the outer loops assure that ishell.ge.jshell.
c        in the case that ishell=jshell(ijsh.eq.true.), we have the full 
c        diagonal block and want to add only the lower triangle only.
c
      do 50 dmat=1,ndmat
         do 40 jf=1,nfj
            jpos=js+jf-nfj
c           jpos=js+(jf-1)*nfj
            iftop=nfi
            if(ijsh) iftop=jf
            do 30 if=1,iftop
               ipos=is+if-nfi
               do 20 j=1,nj
                  jj=jpos+j*nfj
                  itop=ni
                  if(ijsh.and.if.eq.jf) itop=j
                  do 10 i=1,itop
                     ii=ipos+i*nfi
                     ij=max(ii,jj)*(max(ii,jj)-1)/2+min(ii,jj)
                     if(debug) then
                        if(ijsh) then
                           write(iout,*) 'i,j,if,jf,ij',i,j,if,jf,ij,
     $                                    jij(i,j,if,jf,dmat)
                        endif
                     endif
                     jmat(ij,dmat)=jmat(ij,dmat)+jij(i,j,if,jf,dmat)
 10               continue
 20            continue
 30         continue
 40      continue
 50   continue
c
c
      return
      end
