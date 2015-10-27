      program m950                                                              
c                                                                               
c  **************kohn interface code*****************                           
c                                                                               
c   5/1/90  -- changes for symmetry  bhl                                        
c              nl2 read from route and nmotot = nl2 + nsmall                    
c              written to bndmat                                                
c                                                                               
c   7/19/89 -- changes for multiple scattering energies and for                 
c               getting ci energies from mesa files                             
c                                                                               
c   7/23/89 -- changes to read a print option ... bhl at llnl                   
c                                                                               
c                                                                               
c   makes input files for "quad" suite of kohn scattering codes                 
c                                                                               
c  from prototype interface by byron lengsfield.  modified and                  
c  extended by cwm and tnr 7/5/89                                               
c                                                                               
c   complete set of output files specified as of 7/12/89                        
c                                                                               
c    geobas = geometry and basis set                                            
c    denmat = transition density matrices                                       
c    bndmat  = transformation matrix                                            
c              overlap matrix                                                   
c              direct hamiltonians for channel pairs                            
c              hpp-e                                                            
c              hopt(e)                                                          
c                                                                               
      common // ia(1)                                                           
      dimension z(1)                                                            
      equivalence (ia(1),z(1))                                                  
      common /memory/ ioff                                                      
      integer wpadti                                                            
      call drum                                                                 
      call getscm(need,z,maxcor,'m950',0)                                       
      intoff=wpadti(ioff)                                                       
      call pm950(z(ioff),ia(intoff),maxcor)                                     
      call chainx(0)                                                            
      stop                                                                      
      end                                                                       
      subroutine pm950(z,a,maxcor)                                              
      parameter (ibndmats=62)                                                   
      parameter (iexchmat=63)                                                   
      parameter (igeobas=53)                                                    
      parameter (idenmat=61)                                                    
      parameter (maxchan=10)                                                    
      parameter (maxprims=200)                                                  
      implicit integer(a-z)                                                     
      real*8 z(*),stopj,fpkey                                                   
      real*8 time1,time2,time3,cpu,xio,sys                                      
      real*8 set1,set2,sett,energy                                              
      real*8 eroot, esave(200), echan(maxchan)                                  
      dimension a(*)                                                            
      dimension ihdr(10)                                                        
      character*4096 ops                                                        
      character*32 xform                                                        
      character*128 namkohn                                                     
      character*4 itoc                                                          
      character*16 bflabl(2000)                                                 
      logical logkey                                                            
c                                                                               
      data idebug/1/                                                            
c                                                                               
      common /io/inp,iout                                                       
c                                                                               
      write(iout,*)' '                                                          
      write(iout,*)' m950: '                                                    
c                                                                               
      ncore=maxcor                                                              
c                                                                               
      call iosys('read character options from rwf',-1,0,0,ops)                  
c                                                                               
      call iosys('read integer nwks from rwf',1,mdim,0,' ')                     
c                                                                               
      call iosys('read character namkohn from rwf',                             
     $ -1,0,0,namkohn)                                                          
      call iosys('open kohn as old',0,0,0,namkohn)                              
c                                                                               
      call iosys('read real npvec from kohn',1,npvec,0,0,' ')                   
      nsmall=intkey(ops,'kohn=nsmall',nsmall,' ')                               
      write(iout,*)' values of nsmall and npvec'                                
      write(iout,*) nsmall,npvec                                                
c                                                                               
      call iosys('read integer "number of basis functions"'//                   
     $ ' from rwf',1,nbf,0,' ')                                                 
c                                                                               
c..bhl 5/1/90                                                                   
c                                                                               
      nl2=nbf-nsmall                                                            
      nl2=intkey(ops,'kohn=nl2',nl2,' ')                                        
      nmotot=nl2+nsmall                                                         
      write(iout,*)' values of nl2 and nmotot '                                 
      write(iout,*) nl2,nmotot                                                  
c                                                                               
      if(logkey(ops,'drop',.false.,' ')) then                                   
       idrop=1                                                                  
      call iosys('read integer "old nbf" from kohn',                            
     $ 1,oldnbf,0,' ')                                                          
      else                                                                      
       oldnbf=nbf                                                               
       idrop=0                                                                  
      end if                                                                    
c                                                                               
c                                                                               
c..bhl 7/23/89                                                                  
      if(logkey(ops,'kohn=debug',.false.,' ')) then                             
       iprt=1                                                                   
       ipbas=1                                                                  
       ipmo=1                                                                   
       ipov=1                                                                   
       ip1ao=1                                                                  
       ip1mo=1                                                                  
       ipden=1                                                                  
       ip2ao=1                                                                  
       ip2mo=1                                                                  
       iphpp=1                                                                  
       iphopt=1                                                                 
      else                                                                      
       iprt=0                                                                   
       ipbas=0                                                                  
       ipmo=0                                                                   
       ipov=0                                                                   
       ip1ao=0                                                                  
       ip1mo=0                                                                  
       ipden=0                                                                  
       ip2ao=0                                                                  
       ip2mo=0                                                                  
       iphpp=0                                                                  
       iphopt=0                                                                 
      end if                                                                    
      if(logkey(ops,'m950=print=mo',.false.,' ')) ipmo=1                        
      if(logkey(ops,'m950=print=overlap',.false.,' ')) ipov=1                   
      if(logkey(ops,'m950=print=basis',.false.,' ')) ipbas=1                    
      if(logkey(ops,'m950=print=hopt',.false.,' ')) iphopt=1                    
      if(logkey(ops,'m950=print=hpp',.false.,' ')) iphpp=1                      
      if(logkey(ops,'m950=print=htot',.false.,' ')) iphtot=1                    
      if(logkey(ops,'m950=print=density',.false.,' ')) ipden=1                  
      if(logkey(ops,'m950=print=direct',.false.,' ')) then                      
       ip1mo=1                                                                  
       ip2mo=1                                                                  
      end if                                                                    
      if(logkey(ops,'m950=print=direct-ao',.false.,' ')) then                   
       ip1ao=1                                                                  
       ip2ao=1                                                                  
      end if                                                                    
c..bhl 7/23/89                                                                  
      nnp=nbf*(nbf+1)/2                                                         
c                                                                               
      c=1                                                                       
      scr=c+nbf*nbf                                                             
      s=scr+nbf*nbf                                                             
      vec=s+nbf*nbf                                                             
      hopt=vec+nbf*nbf                                                          
      hpp=hopt+npvec*npvec                                                      
      dir1=hpp+npvec*npvec                                                      
      dir2=dir1+nnp                                                             
      exch=dir2+nnp                                                             
      need=exch+nbf*nbf                                                         
      nroots=intkey(ops,'ci=nroots',1,' ')                                      
c                                                                               
c create binary output files to be read by quad kohn codes                      
c                                                                               
      open(igeobas,file='geobas',form='unformatted',status='unknown')           
      open(22,file='incint',form='unformatted',status='unknown')                
      open(idenmat,file='denmat',form='unformatted',status='unknown')           
      open(ibndmats,file='bndmat',form='unformatted',status='unknown')          
      open(iexchmat,file='exchmat',form='unformatted',status='unknown')         
c                                                                               
      maxlen = nbf*nbf*(maxchan+1)*maxchan/2 + 10*maxchan                       
      nenergy = 1                                                               
      maxbnd = (nroots*(nroots+1)/2 + 1)*nbf*nbf +                              
     x  nenergy*2*npvec*npvec + nnp + nroots + 20                               
c                                                                               
c      call create(igeobas,6hgeobas,1,maxprims*8)                               
c      call create(idenmat,6hdenmat,1,maxlen)                                   
c      call create(ibndmats,6hbndmat ,1,maxbnd)                                 
c                                                                               
c                                                                               
c  write the preliminary information to bndmats                                 
c  this includes the number of channels and their ci energies                   
c  and the number of scattering energies in this run                            
c                                                                               
      do 7001 i=1,nroots                                                        
      call iosys('read real "ci energy '//itoc(i)//'" from kohn',               
     x 1,eroot,0,' ')                                                           
 7001  echan(i) = eroot                                                         
      write(ibndmats) nroots                                                    
      write(iexchmat) nroots                                                    
      write(ibndmats) (echan(i),i=1,nroots)                                     
      write(iout,7003) (echan(i),i=1,nroots)                                    
 7003 format(' channel energies from ci codes :',/,(10x,5(2x,e15.8)))           
c                                                                               
      call iosys('read integer nenergy from kohn',1,numpts,0,' ')               
      call iosys('read real "scattering energies" from kohn',                   
     x numpts,esave,0,' ')                                                      
      write(ibndmats) numpts                                                    
      write(ibndmats) (esave(i),i=1,numpts)                                     
      write(iout,7004) (esave(i),i=1,numpts)                                    
 7004 format(' incident energies from kohn file:',/,5(2x,e15.8))                
c                                                                               
      call iosys('read character "transformation vector" from kohn',            
     $ -1,0,0,xform)                                                            
      call iosys('read real '//xform//' from kohn',nbf*nbf,z(c),0,0)            
      write(iout,*)' '                                                          
      write(iout,*)'  molecular orbitals '                                      
      if(ipmo.ne.0) then                                                        
      call matout(z(c),nbf,nbf,nbf,nbf,iout)                                    
      end if                                                                    
c                                                                               
c  write transformation matrix to bndmats for kohnopt                           
c                                                                               
c..bhl 5/1/90                                                                   
      write(ibndmats) nbf,nmotot                                                
      call wrbinsqr(z(c),nbf,nbf,ibndmats)                                      
c..bhl 5/1/90                                                                   
c                                                                               
      call iosys('read real "overlap integrals" from kohn',                     
     $ nnp,z(scr),0,' ')                                                        
c                                                                               
c                                                                               
c  write overlap matrix to bndmats for kohnopt                                  
      call wrbintri(z(scr),nbf,ibndmats)                                        
c                                                                               
      call trtosq(z(s),z(scr),nbf,nnp)                                          
      write(iout,*)' '                                                          
      write(iout,*)' overlap matrix '                                           
      if(ipov.ne.0) then                                                        
      call matout(z(s),nbf,nbf,nbf,nbf,iout)                                    
      end if                                                                    
c                                                                               
      write(iout,*)' '                                                          
      write(iout,*)' direct ao 1 electron kohn hamiltonian'                     
      call iosys('read real "direct ao 1e kohn hamiltonian"'//                  
     $' from kohn',nnp,z(dir1),0,' ')                                           
      call trtosq(z(scr),z(dir1),nbf,nnp)                                       
      if(ip1ao.ne.0) then                                                       
      call matout(z(scr),nbf,nbf,nbf,nbf,iout)                                  
      end if                                                                    
c                                                                               
      call ebc(z(vec),z(scr),z(c),nbf,nbf,nbf)                                  
      call ebtc(z(scr),z(c),z(vec),nbf,nbf,nbf)                                 
      write(iout,*)' direct mo 1 electron kohn hamiltonian'                     
      if(ip1mo.ne.0) then                                                       
      call matout(z(scr),nbf,nbf,nbf,nbf,iout)                                  
      end if                                                                    
c                                                                               
c                                                                               
      call iosys('rewind "direct ao 2e kohn hamiltonian" on kohn',              
     $ 0,0,0,' ')                                                               
      call iosys('rewind "exchange ao 2e kohn hamiltonian" on kohn',            
     $ 0,0,0,' ')                                                               
c                                                                               
c  write number of channels to file for quad codes                              
c                                                                               
      write(idenmat)nroots,nbf                                                  
      write(ibndmats) nroots                                                    
c                                                                               
      do 1 i=1,nroots                                                           
       do 2 j=1,i                                                               
        call iosys('read real "ao t1pdm:'//itoc(i)//                            
     $  itoc(j)//'" from kohn',nbf*nbf,z(scr),0,' ')                            
        write(iout,*)' '                                                        
        write(iout,*)' ao transition density for states ',i,j                   
      if(ipden.ne.0) then                                                       
        call matout(z(scr),nbf,nbf,nbf,nbf,iout)                                
      end if                                                                    
c                                                                               
c                                                                               
c write density matrix for channel pair i,j to file denmat for quad codes       
c                                                                               
      write(idenmat) i,j                                                        
      write(ibndmats) i,j                                                       
      write(iexchmat) i,j                                                       
      call wrbinsqr(z(scr),nbf,nbf,idenmat)                                     
c                                                                               
      write(iout,*)' '                                                          
      write(iout,*)' direct ao (1+2 electron) kohn hamiltonian'                 
      call iosys('read real "direct ao 2e kohn hamiltonian"'//                  
     $' from kohn without rewinding',nnp,z(dir2),0,' ')                         
      call iosys('read real "exchange ao 2e kohn hamiltonian"'//                
     $ ' from kohn without rewinding',nbf*nbf,z(exch),0,' ')                    
c                                                                               
c  if channel i = channel j add one-electron contribution to                    
c  direct hamiltonian                                                           
c                                                                               
      if(i.eq.j) then                                                           
      do 6 iii=1,nnp                                                            
  6   z(dir2+iii-1) = z(dir2+iii-1) + z(dir1+iii-1)                             
      endif                                                                     
c                                                                               
      call trtosq(z(scr),z(dir2),nbf,nnp)                                       
      if(ip2ao.ne.0) then                                                       
      call matout(z(scr),nbf,nbf,nbf,nbf,iout)                                  
      end if                                                                    
c                                                                               
c  write direct hamiltonian for this channel pair to file for kohnopt           
c                                                                               
      call wrbinsqr(z(scr),nbf,nbf,ibndmats)                                    
c                                                                               
c  write exch hamiltonian for this channel pair to file for kohnopt             
c                                                                               
      call wrbinsqr(z(exch),nbf,nbf,iexchmat)                                   
c                                                                               
      call ebc(z(vec),z(scr),z(c),nbf,nbf,nbf)                                  
      call ebtc(z(scr),z(c),z(vec),nbf,nbf,nbf)                                 
      write(iout,*)' direct mo (1+2 electron) kohn hamiltonian'                 
      if(ip2mo.ne.0) then                                                       
      call matout(z(scr),nbf,nbf,nbf,nbf,iout)                                  
      end if                                                                    
c                                                                               
  2    continue                                                                 
 1    continue                                                                  
c                                                                               
c write the hpp-e and hopt matrices to file for kohnopt                         
c                                                                               
      call iosys('rewind hopt on kohn',0,0,0,' ')                               
      call iosys('rewind hpp on kohn',0,0,0,' ')                                
c                                                                               
c loop on incident energies                                                     
c                                                                               
      do 600 iene=1,numpts                                                      
      write(ibndmats) npvec, nsmall                                             
      ntot=npvec*npvec                                                          
      call iosys('read real hpp from kohn without rewinding',                   
     x ntot,z(hpp),0,' ')                                                       
      write(iout,*)' '                                                          
      write(iout,*)'  hpp-e '                                                   
      if(iphpp.ne.0) then                                                       
      call matout(z(hpp),npvec,npvec,npvec,npvec,iout)                          
      end if                                                                    
c                                                                               
      call iosys('read real hopt from kohn without rewinding',                  
     x ntot,z(hopt),0,' ')                                                      
c                                                                               
c  mesa computes hopt = hpq (hqq - e)**-1 hqp                                   
c  so in has the wrong overall sign.  therefore change its                      
c  sign here.                                                                   
c                                                                               
      npsqr = npvec*npvec                                                       
      do 87 iip = 1, npsqr                                                      
  87  z(hopt+iip-1) = -z(hopt+iip-1)                                            
c                                                                               
      write(iout,*)' '                                                          
      write(iout,*)'  hopt '                                                    
      if(iphopt.ne.0) then                                                      
      call matout(z(hopt),npvec,npvec,npvec,npvec,iout)                         
      end if                                                                    
c                                                                               
      call wrbinsqr(z(hpp),npvec,npvec,ibndmats)                                
      call wrbinsqr(z(hopt),npvec,npvec,ibndmats)                               
c                                                                               
      if(iphtot.ne.0) then                                                      
      write(iout,*)' '                                                          
      write(iout,*)'  hopt + hpp '                                              
      do 88 iip=1,npsqr                                                         
       z(hopt+iip+1)=z(hopt+iip-1)+z(hpp+iip-1)                                 
  88  continue                                                                  
      call matout(z(hopt),npvec,npvec,npvec,npvec,iout)                         
      end if                                                                    
c                                                                               
c close loop on incident energies                                               
  600 continue                                                                  
c                                                                               
c                                                                               
c     get the basis function labels.                                            
c                                                                               
      call iosys('read character "basis function labels" from rwf',             
     $          -1,0,0,bflabl)                                                  
c                                                                               
c     get the lengths of arrays needed for core allocation.                     
c     ntypes is the total number of types used to define the lengths            
c     in l102.  this total is composed of two sets: the first nbtype            
c     types are used to mark the basis function types.  the remainder           
c     refer to the ecp types.                                                   
c                                                                               
      call iosys('read integer "number of atoms" from rwf',1,nat,0,0)           
      call iosys('read integer "number of basis functions" from rwf',           
     $     1,nbasis,0,0)                                                        
      call iosys('length of exponents on rwf',nprim,0,0,' ')                    
      call iosys('length of "contraction coefficients" on rwf',                 
     $     ncont,0,0,' ')                                                       
      call iosys('length of "number of pure functions" on rwf',                 
     $     ntypes,0,0,' ')                                                      
      call iosys('read integer "number basis types" from rwf',                  
     $     1,nbtype,0,0)                                                        
      call iosys('length of "power of x" on rwf',lenxyz,0,0,' ')                
      nnp=(nbasis+1)*nbasis/2                                                   
c                                                                               
c     divide core for basis set information.                                    
c                                                                               
      index=1                                                                   
      zan=index+oldnbf                                                          
      c=zan+nat                                                                 
      cont=c+3*nat                                                              
      ex=cont+ncont                                                             
      ptprim=wpadti(ex+nprim)                                                   
      noprim=ptprim+ntypes*nat                                                  
      nocont=noprim+ntypes*nat                                                  
      ptcont=nocont+ntypes*nat                                                  
      start=ptcont+ntypes*nat                                                   
      nocart=start+ntypes*nat                                                   
      nobf=nocart+ntypes                                                        
      minmom=nobf+ntypes                                                        
      maxmom=minmom+ntypes                                                      
      mintyp=maxmom+ntypes                                                      
      nx=mintyp+ntypes                                                          
      ny=nx+lenxyz                                                              
      nz=ny+lenxyz                                                              
      s=iadtwp(nz+lenxyz)                                                       
c                                                                               
c     determine the amount of scratch space required.                           
c     retrieve information about the most demanding shell block.                
c                                                                               
      call iosys('read integer maxprm from rwf',1,maxprm,0,0)                   
      call iosys('read integer maxcont from rwf',1,mxcont,0,0)                  
      call iosys('read integer maxl from rwf',1,maxl,0,0)                       
      call iosys('read integer maxblk from rwf',1,maxblk,0,0)                   
      call iosys('read integer d1maxblk from rwf',1,dlen,0,0)                   
      call iosys('read integer dolp from rwf',1,dolp,0,0)                       
c                                                                               
c                                                                               
c     read in basis set information from read-write file.                       
c                                                                               
      call iosys('read real exponents from rwf',-1,z(ex),0,0)                   
      call iosys('read real "contraction coefficients" from rwf',               
     $     -1,z(cont),0,0)                                                      
      call iosys('read real "nuclear charges" from rwf',-1,z(zan),0,0)          
      call iosys('read real coordinates from rwf',-1,z(c),0,0)                  
      call iosys('read integer "pointer to primitives" from rwf',               
     $     -1,a(ptprim),0,0)                                                    
      call iosys('read integer "number of primitives" from rwf',                
     $     -1,a(noprim),0,0)                                                    
      call iosys('read integer "pointer to contraction coefficients"'//         
     $     ' from rwf',-1,a(ptcont),0,0)                                        
      call iosys('read integer "number of contraction coefficients" '//         
     $     'from rwf',-1,a(nocont),0,0)                                         
      call iosys('read integer "number of cartesians" from rwf',                
     $     -1,a(nocart),0,0)                                                    
      call iosys('read integer "number of pure functions" from rwf',            
     $     -1,a(nobf),0,0)                                                      
      call iosys('read integer "minimum momentum" from rwf',                    
     $     -1,a(minmom),0,0)                                                    
      call iosys('read integer "maximum momentum" from rwf',                    
     $     -1,a(maxmom),0,0)                                                    
      call iosys('read integer "pointer to cartesians" from rwf',               
     $     -1,a(mintyp),0,0)                                                    
      call iosys('read integer "pointer to first function" from rwf',           
     $     -1,a(start),0,0)                                                     
      call iosys('read integer "power of x" from rwf',-1,a(nx),0,0)             
      call iosys('read integer "power of y" from rwf',-1,a(ny),0,0)             
      call iosys('read integer "power of z" from rwf',-1,a(nz),0,0)             
c                                                                               
      write(iout,*)' '                                                          
      write(iout,*)' geometry '                                                 
      ix=c-1                                                                    
      write(igeobas) nat                                                        
      do 10 i=1,nat                                                             
       write(iout,11) z(ix+1),z(ix+2),z(ix+3)                                   
      write(igeobas) z(ix+1),z(ix+2),z(ix+3),z(zan-1+i)                         
       ix=ix+3                                                                  
  10  continue                                                                  
  11  format(3(2x,f14.9))                                                       
                                                                                
c                                                                               
c     calculate one-electron multipole integrals.                               
c                                                                               
c                                                                               
c..bhl 9/13/89  passing index idrop and oldnbf                                  
c                                                                               
         call basout(z(c),z(ex),z(cont),a(ptprim),                              
     #            a(noprim),a(nocont),a(ptcont),nat,nprim,maxcor-need+1,        
     #            ntypes,nbtype,nnp,ncont,a(start),nbasis,z(zan),               
     #            a(nocart),a(nobf),a(maxmom),a(mintyp),a(nx),a(ny),            
     # a(nz),a(minmom),igeobas,idebug,ipbas,a(index),idrop,oldnbf)              
c                                                                               
      call chainx(0)                                                            
c                                                                               
      stop                                                                      
      end                                                                       
*deck basout                                                                    
      subroutine basout(coords,ex,cont,ptprim,noprim,nocont,ptcont,             
     #           nat,nprim,maxcor,ntypes,nbtype,nnp,ncont,                      
     #           start,nbasis,zan,nocart,nobf,maxmom,mintyp,                    
     #   nx,ny,nz,minmom,igeobas,idebug,iprt,index,idrop,oldnbf)                
c                                                                               
c                                                                               
      implicit integer (a-z)                                                    
c                                                                               
      character*1 itoc                                                          
      character*4 funcnm                                                        
      logical logkey                                                            
      real*8 ex(nprim),cont(ncont)                                              
      real*8 zan(nat)                                                           
      integer ptprim(nat,ntypes),noprim(nat,ntypes),nocont(nat,ntypes)          
      integer ptcont(nat,ntypes),start(nat,ntypes),nocart(ntypes)               
      integer nobf(ntypes),maxmom(ntypes),mintyp(ntypes),minmom(ntypes)         
      integer nx(*),ny(*),nz(*)                                                 
      integer index(oldnbf)                                                     
c                                                                               
c  dimensions for the arrays being built for quad codes                         
c                                                                               
      real*8 coords(3*nat),eta(1000,5),x(3,500),alpha(20,500)                   
      real*8 cc(20,500)                                                         
      integer nstart(1000), nstop(1000),s(500),ncontt(500)                      
      integer lnew(1000),mnew(1000),nnew(1000),ncartt(500)                      
      integer ll(15,500),mm(15,500),nn(15,500)                                  
c                                                                               
c                                                                               
      common/io/inp,iout                                                        
c                                                                               
      nstart(1) = 1                                                             
      ikount = 0                                                                
      iloc = 0                                                                  
c                                                                               
c.. bhl 9/13/89         using index array to drop functions                     
c                                                                               
      if(idrop.ne.0) then                                                       
      write(iout,*)' reading the drop index vector from rwf'                    
      call iosys('read integer "packing index vector" from kohn',               
     $ oldnbf,index,0,' ')                                                      
      else                                                                      
      do 999 i=1,oldnbf                                                         
      index(i)=i                                                                
 999  continue                                                                  
      end if                                                                    
c                                                                               
      kbf=0                                                                     
c                                                                               
c.. bhl 9/13/89             kbf counter                                         
c                                                                               
      ncent=0                                                                   
      do 9 iatom=1,nat                                                          
c        write(iout,*)' '                                                       
c        write(iout,*)' atom number ',iatom                                     
c        write(iout,*)' charge      ',zan(iatom)                                
             do 7 itype=1,nbtype                                                
               if (noprim(iatom,itype).le.0) go to 7                            
c                                                                               
                   mini=mintyp(itype)                                           
                   maxi=mintyp(itype)+nocart(itype)-1                           
                   imax=maxmom(itype)                                           
c        write(iout,*)'  primitive type s=1 p=2 d=3 etc. ',itype                
c        write(iout,*)'  maximum momemtum ',imax                                
c                                                                               
                   nprimi=noprim(iatom,itype)                                   
                   nconti=nocont(iatom,itype)                                   
c        write(iout,*)'    number of primitives functions ',nprimi              
c        write(iout,*)'    number of contracted functions ',nconti              
c                                                                               
c  build arrays for quad codes: nstart, nstop, lnew, mnew, nnew, eta            
c                                                                               
c                                                                               
c  the way the  indexing for quad codes works is that there are                 
c  ncontra contracted functions.  for the ith contracted function nstart(i)     
c  and nstop(i) give the first and last locations of the relevant quantities    
c  in the arrays lnew(i), mnew(i), nnew(i), and eta(i,j=1,5)                    
c   lnew, mnew, and new are the powers of x, y, and z of the primitive          
c   and eta(i,j=1,5) contains the center, exponent and contraction coefficient  
c                                                                               
      do 6 icont=1,nconti                                                       
      ibf=icont+ncent                                                           
      s(ibf)=itype-1                                                            
      x(1,ibf)=coords(3*(iatom-1)+1)                                            
      x(2,ibf)=coords(3*(iatom-1)+2)                                            
      x(3,ibf)=coords(3*(iatom-1)+3)                                            
      ncartt(ibf)=0                                                             
       do 3 m=mini,maxi                                                         
c..bhl 9/13/89                  kbf counter                                     
      kbf=kbf+1                                                                 
      if(index(kbf).ne.0) then                                                  
      ncartt(ibf)=ncartt(ibf)+1                                                 
      ll(ncartt(ibf),ibf)=nx(m)                                                 
      mm(ncartt(ibf),ibf)=ny(m)                                                 
      nn(ncartt(ibf),ibf)=nz(m)                                                 
c..bhl 9/13/89                  if statement                                    
      ikount = ikount + 1                                                       
      iprim=0                                                                   
      do 4 iii=1,nprimi                                                         
      iiirel = iii-1+ (icont-1)*nprimi                                          
      if (abs(cont(ptcont(iatom,itype)+iiirel)).gt.1.e-20) then                 
       iloc = iloc+1                                                            
       iprim=iprim+1                                                            
       eta(iloc,1) = coords(3*(iatom-1)+1)                                      
       eta(iloc,2) = coords(3*(iatom-1)+2)                                      
       eta(iloc,3) = coords(3*(iatom-1)+3)                                      
       eta(iloc,4) = ex(ptprim(iatom,itype)+iii - 1)                            
       eta(iloc,5) = cont(ptcont(iatom,itype)+iiirel)                           
       lnew(iloc) = nx(m)                                                       
       mnew(iloc) = ny(m)                                                       
       nnew(iloc) = nz(m)                                                       
       alpha(iprim,ibf)=eta(iloc,4)                                             
       cc(iprim,ibf)=eta(iloc,5)                                                
      endif                                                                     
   4  continue                                                                  
      nstop(ikount) = iloc                                                      
      nstart(ikount+1) = nstop(ikount)+ 1                                       
      end if                                                                    
c..bhl 9/13/89                  end if                                          
  3    continue                                                                 
      ncontt(ibf)=iprim                                                         
   6  continue                                                                  
                   ncent=ncent+nconti                                           
c                                                                               
c        write(iout,*)' '                                                       
c        write(iout,*)'    primitive exponents '                                
      if(iprt.ne.0) then                                                        
        call matout(ex(ptprim(iatom,itype)),nprimi,1,                           
     $           nprimi,1,iout)                                                 
      end if                                                                    
c        write(iout,*)' '                                                       
c        write(iout,*)'    contraction coefficients '                           
      if(iprt.ne.0) then                                                        
        call matout(cont(ptcont(iatom,itype)),nprimi,nconti,                    
     $           nprimi,nconti,iout)                                            
      end if                                                                    
c        write(iout,*)' '                                                       
c        write(iout,*)'    xyx convention '                                     
        do 5 m=mini,maxi                                                        
        if(iprt.ne.0)write(iout,*)nx(m),ny(m),nz(m)                             
    5   continue                                                                
    7       continue                                                            
    9 continue                                                                  
      write(22)ncent                                                            
      do 20 i1=1,ncent                                                          
      write(22)s(i1),ncontt(i1),(x(j1,i1),j1=1,3)                               
      write(22)(alpha(j1,i1),j1=1,ncontt(i1))                                   
      write(22)(cc(k1,i1),k1=1,ncontt(i1))                                      
      write(22)ncartt(i1)                                                       
      write(22)(ll(j1,i1),mm(j1,i1),nn(j1,i1),j1=1,ncartt(i1))                  
 20   continue                                                                  
 101  format(i5)                                                                
 102  format(i6,10f6.3)                                                         
 103  format(i7,10f10.3)                                                        
 104  format(i8)                                                                
 105  format(5x,3i5)                                                            
 106   format(10e12.4)                                                          
      neta = iloc                                                               
      ncontra = ikount                                                          
c                                                                               
c write the newly indexed basis set information to binary file igeobas          
c                                                                               
c                                                                               
      write(igeobas) ncontra                                                    
      write(igeobas) (nstart(i),nstop(i),i=1,ncontra)                           
      write(iout,*)' '                                                          
      write(iout,*)' geometry and basis written to file ',igeobas               
      if(idebug.ne.0) then                                                      
      write(iout,*)' nbf ',ncontra                                              
      do 12 i=1,ncontra                                                         
      nlow = nstart(i)                                                          
      nhi = nstop(i)                                                            
      do 11 j=nlow,nhi                                                          
      write(igeobas) lnew(j),mnew(j),nnew(j),(eta(j,jj),jj=1,5)                 
      if(iprt.ne.0) then                                                        
      write(iout,*) lnew(j),mnew(j),nnew(j),(eta(j,jj),jj=1,5)                  
      end if                                                                    
  11   continue                                                                 
  12   continue                                                                 
      endif                                                                     
c                                                                               
c                                                                               
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
      subroutine wrbinsqr(amat,nrow,ncol,iunit)                                 
      dimension amat(nrow,ncol)                                                 
      do 1 i=1,ncol                                                             
  1   write(iunit) (amat(j,i),j=1,ncol)                                         
      return                                                                    
      end                                                                       
      subroutine wrbintri(amat,nrow,iunit)                                      
      dimension amat(nrow*(nrow+1)/2)                                           
      nn=nrow*(nrow+1)/2                                                        
      write(iunit) (amat(i),i=1,nn)                                             
      return                                                                    
      end                                                                       
