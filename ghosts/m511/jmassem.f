*deck %W% %G%
      subroutine jmassem(conint,nconti,ncontj,ncontk,ncontl,
     $     lenblk,symcen,angmom,start,nat,nbtype,nnp,lenbuf,numbuf,
     $     ncint,numint,nbasis,nobf,nactul,jmat,ncoul,d,ndmat)
c
c***begin prologue     %M%
c***date written       930829
c***revision date      %G%
c
c***keywords           m511, link 511, direct coulomb
c***author             russo, thomas (lanl)
c***source             %W% %G%
c***purpose            assemble the coulomb matrix from integrals in conint
c***description        integrals come in in conint, (ij|kl) where i,j,k,l
c                      are functions on symcen(1) to symcen(4), have angular
c                      momentum type given  by angmom(1) to angmom(4)
c                      We must figure out what function number those correspond
c                      to, multiply by the appropriate density matrix pieces,
c                      form the appropriate elements of the coulomb matrix.
c
c***references        
c                      We don't need no steeeenking references
c
c***routines called
c***end prologue       %M%
c
      implicit none
      real*8 cutoff
      integer nconti,ncontj,ncontk,ncontl,lenblk,nat,nbtype,nnp,lenbuf,
     $     numbuf,ncint,numint,nbasis,nactul,ncoul,ndmat
      integer symcen(4),angmom(4),start(nat,nbtype),nobf(nbtype)

      real*8 conint(numint),jmat(nnp,ncoul),d(nnp,ndmat)
      integer numi,numj,numk,numl,inp,iout,istart,jstart,kstart,lstart
      integer ic,jc,kc,lc,icrt,jcrt,kcrt,lcrt
      integer i,j,k,l,ij,kl,theint,shell
      real*8 mul

      common/toler/cutoff
      common/io/inp,iout

c
c how many functions are on each shell?
c
      numi=nobf(angmom(1))
      numj=nobf(angmom(2))
      numk=nobf(angmom(3))
      numl=nobf(angmom(4))
c
c the *start variables give the basis function number of the first of each block
c
      istart=start(symcen(1),angmom(1))
      jstart=start(symcen(2),angmom(2))
      kstart=start(symcen(3),angmom(3))
      lstart=start(symcen(4),angmom(4))
c
c we now have a buffer of (ij|kl) where i goes from istart+1 to 
c istart+numi*nconti etc. BUT IN A FUNKY ORDER!  Yech.  The way it works is
c that all cartesians for a given contraction are in a block, i.e. if
c we have 2 s contractions and 1 p, integrals (ps|ps) are stored in this
c order:
c (xs1|xs1) (xs1|ys1) (xs2|zs1) (ys1|xs1) (...) (zs1|zs1) (xs1|xs2)...
c
c ok, this rots bad, but maybe it will work.  We can do it more efficiently
c when we're sure the REST of the code works.
c
c pick a j matrix element to work on
c
      theint=1
c this is the broken code:
c$$$      do 10 i=istart+1,istart+numi*nconti
c$$$         do 20 j=jstart+1,jstart+numj*ncontj
c$$$            if (i .ge. j) then
c$$$               ij=((i-1)*i)/2+j
c$$$            else
c$$$               theint=theint+numl*ncontl*numk*ncontk
c$$$               goto 20
c$$$            endif
c$$$            do 30 k=kstart+1,kstart+numk*ncontk
c$$$               do 40 l=lstart+1,lstart+numl*ncontl
c$$$                  if (k .ge. l) then
c$$$                     kl=((k-1)*k)/2+l
c$$$                  else
c$$$                     kl=((l-1)*l)/2+k
c$$$                  endif
c PROBLEM:  Canonical integrals don't include stuff like (ss|ps), only
c (ps|ss), so using (ij| for the J matrix index will not get all the 
c contributions needed: if i and j are S functions then we only get the 
c (ss|ss) contributions.  
      do 10 lc=1,ncontl
         do 20 kc=1,ncontk
            do 30 jc=1,ncontj
               do 40 ic=1,nconti
                  do 50 icrt=1,numi
                     do 60 jcrt=1,numj
                        i=istart+(ic-1)*numi+icrt
                        j=jstart+(jc-1)*numj+jcrt
                        if (i .ge. j) then
                           ij=((i-1)*i)/2+j
                           do 70 kcrt=1,numk
                              do 80 lcrt=1,numl
                                 k=kstart+(kc-1)*numk+kcrt
                                 l=lstart+(lc-1)*numl+lcrt
                                 if (k .ge. l) then
                                    kl=((k-1)*k)/2+l
                                 else
                                    kl=((l-1)*l)/2+k
                                 endif
                                 write(iout,*)'(',i,',',j,'|',k,',',l,
     $                                ')=conint(',theint,')=',
     $                                conint(theint)
                                 do 1000 shell=1,ncoul
                                    jmat(ij,shell)=jmat(ij,shell)+
     $                                   d(kl,shell)*conint(theint)
 1000                            continue 
                                 theint=theint+1
 80                           continue 
 70                        continue 
                        else
                           theint=theint+numl*numk
                        endif
 60                  continue 
 50               continue 
 40            continue 
 30         continue 
 20      continue 
 10   continue  
      return
      end
