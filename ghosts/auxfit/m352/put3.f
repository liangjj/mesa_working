*deck  @(#)put3.f	1.2 4/28/92
      subroutine put3(s,conint,start,astart,iatom,jatom,katom,
     $                  itype,jtype,ktype,nconti,ncontj,ncontk,
     #                  nnp,lenblk,nat,nbtype,nobf)
c
c***begin prologue     put3
c***date written       840723  
c***revision date      920417
c      april 17, 1992  rlm at lanl
c         modified to place three-center integrals in appropriate target block.
c***keywords           
c***author             saxe, paul and martin, richard(lanl)
c***source             @(#)put3.f	1.1   4/27/92
c***purpose            transfers an angular-momentum block of three-center
c                      integrals from conint to the array s.
c***description        modified from routine m302/put1el
c     
c    
c
c***references
c
c***routines called    none
c
c***end prologue       put3
c
      implicit integer (a-z)
c
c     ----- input arrays(unmodified) -----
      real*8 conint(nconti,ncontj,ncontk,lenblk)
      integer start(nat,nbtype),astart(nat,nbtype),nobf(nbtype)
c
c     ----- output arrays -----
      real*8 s(nnp,*)
c
      common/io/inp,iout
c
c
      numi=nobf(itype)
      numj=nobf(jtype)
      numk=nobf(ktype)
      istart=start(iatom,itype)
      jstart=start(jatom,jtype)
      kstart=astart(katom,ktype)
c
      intgrl=0
      do 60 kf=1,numk
         do 50 if=1,numi
            do 40 jf=1,numj
               intgrl=intgrl+1
               do 30 jc=1,ncontj
                  j=jstart+(jc-1)*numj+jf
                  do 20 ic=1,nconti
                     i=istart+(ic-1)*numi+if
c    i think we do too much work here
                     if (i.ge.j) then
                        ij=i*(i-1)/2+j
                     else
                        ij=j*(j-1)/2+i
                     end if
                     do 10 kc=1,ncontk
                        k=kstart+(kc-1)*numk+kf
                        s(ij,k)=conint(ic,jc,kc,intgrl)
   10                continue
   20             continue
   30          continue
   40       continue
   50    continue
   60 continue
c
c
      return
      end
