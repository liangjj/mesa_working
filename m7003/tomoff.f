*deck @(#)tomoff.f	1.1 9/8/91
c***begin prologue     m7001
c***date written       920525   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m7001, link 7001
c***author             schneider, b. (nsf)
c***source             m7001
c***purpose            transform free functions to basis orthogonal
c***                   to bound and complex functions.
c***
c***description        the free functions and ( H0 - E ) are schmidt
c***                   orthogonalized to the bound and complex orbitals.
c***                   note that the arrays sc and sr are equivalenced
c***                   in the calling program.
c***                 
c***references       
c
c***routines called  
c***end prologue       m7001
      subroutine tomoff(frefn0,frefn1,ddfre0,ddfre1,fns,ddfns,fnsc,
     1                  ddfnsc,s,sc,sc1,sc2,sr,sr1,sr2,mask,
     2                  ecdiag,lf,mf,lb,mb,lc,mc,nlmfre,nchan,n,sze,
     3                  nbfb,nbfc,prnt)
      implicit integer (a-z)
      real *8 fns, ddfns, mask, frefn1, ddfre1, ddfre0
      real *8 ecdiag, sr, sr1, sr2
      complex *16 fnsc, ddfnsc, frefn0, sc, sc1, sc2, s
      logical prnt
      dimension frefn0(n,sze), ddfre0(n,sze), fns(n,nbfb), ddfns(n,nbfb)
      dimension fnsc(n,nbfc), ddfnsc(n,nbfc), frefn1(n,sze)
      dimension mask(*), lb(nbfb), mb(nbfb), lc(nbfc), mc(nbfc)
      dimension lf(sze), mf(sze), ecdiag(sze), sc(n,*), sc1(n,*)
      dimension nlmfre(*), s(nbfb,sze), ddfre1(n,sze), sc2(n,*)
      dimension sr(n,*), sr1(n,*), sr2(n,*)
      common /io/ inp, iout
c**********************************************************************c
c              put the non complex arrays into temporary complex       c
c              arrays. after the orthogonalization they can be         c
c              complex and space was allowed for that in the calling   c
c                             routine                                  c
c**********************************************************************c
c**********************************************************************c
c                schmidt the free to the bound orbitals.               c
c                note that due to the presence of the channel energy   c
c                in the definition of ( H0 - E ) for the free          c
c                functions it is necessary to modify the definition    c
c                of the derivatives before the transformation.         c
c**********************************************************************c
      if (nbfb.ne.0) then
          call mskone(mask,lb,mb,nbfb,lf,mf,sze)
          call ebtcc(s,fns,frefn0,nbfb,n,sze)
          call cvmul(s,s,mask,mask,nbfb*sze,'real')
          call ambcc(frefn0,fns,s,n,nbfb,sze)
          call ebtc(sr,fns,frefn1,nbfb,n,sze)
          call vmul(sr,sr,mask,nbfb*sze)
          call ambc(frefn1,fns,sr,n,nbfb,sze)
          call rtocm(ddfre0,sc1,n*sze) 
          icloc=1 
          do 30 i=1,nchan
             do 40 j=1,nbfb
                do 50 k=1,n
                   ddfns(k,j)=ddfns(k,j)-ecdiag(i)*fns(k,j)
   50           continue
   40        continue
             call ambcc(sc1(1,icloc),ddfns,s(1,icloc),n,nbfb,nlmfre(i))
             call ambc(ddfre1(1,icloc),ddfns,sr(1,icloc),n,nbfb,
     1                 nlmfre(i)) 
             icloc=icloc+nlmfre(i)
   30     continue      
          call cc2opy(sc1,ddfre0,n*sze)
          call rtocm(frefn1,sc1,n*sze)
          call rtocm(ddfre1,sc2,n*sze)
          call cc2opy(sc1,frefn1,n*sze)
          call cc2opy(sc2,ddfre1,n*sze)
c**********************************************************************c
c             at this point all arrays are complex                     c
c                    due to the orthogonalization                      c
c**********************************************************************c
      endif
c**********************************************************************c
c                same as above for complex functions                   c
c**********************************************************************c 
      if (nbfc.ne.0) then
c**********************************************************************c
c                in case first if on nbfb not executed                 c
c**********************************************************************c
          if (nbfb.eq.0) then
              call rtocm(ddfre0,sc1,n*sze)
              call cc2opy(sc1,ddfre0,n*sze)
              call rtocm(frefn1,sc1,n*sze)
              call cc2opy(sc1,frefn1,n*sze)
              call rtocm(ddfre1,sc1,n*sze)
              call cc2opy(sc1,ddfre1,n*sze)
          endif
c**********************************************************************c
c                       proceed                                        c
c**********************************************************************c
          call mskone(mask,lc,mc,nbfc,lf,mf,sze)
          call cebtc(s,fnsc,frefn0,nbfc,n,sze)
          call cvmul(s,s,mask,mask,nbfc*sze,'real')
          call cambc(frefn0,fnsc,s,n,nbfc,sze)
          call cebtc(sc,fnsc,frefn1,nbfc,n,sze)
          call cvmul(sc,sc,mask,mask,nbfc*sze,'real')
          call cambc(frefn1,fnsc,sc,n,nbfc,sze)
          call cc2opy(ddfre0,sc1,n*sze)
          call cc2opy(ddfre1,sc2,n*sze)
          icloc=1 
          do 100 i=1,nchan
             do 110 j=1,nbfb
                do 120 k=1,n
                   ddfnsc(k,j)=ddfnsc(k,j)-ecdiag(i)*fnsc(k,j)
  120           continue
  110        continue
             call cambc(sc1(1,icloc),ddfnsc,s(1,icloc),n,nbfc,nlmfre(i))
             call cambc(sc2(1,icloc),ddfnsc,sc(1,icloc),n,nbfc,
     1                  nlmfre(i))
             icloc=icloc+nlmfre(i)
  100     continue      
          call cc2opy(sc1,ddfre0,n*sze)
          call cc2opy(sc2,ddfre1,n*sze)
      endif
      return
      end









