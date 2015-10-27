*deck @(#)m330.f	5.1  11/6/94
      program m330
c***begin prologue     m330.f
c***date written       840707  
c***revision date      11/6/94      
c
c***keywords           
c***author             saxe, paul(lanl)
c***source             @(#)m330.f	5.1   11/6/94
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       m330.f
      implicit none
c     --- input variables -----
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
      integer a
      real*8 z
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer maxnbf
      parameter (maxnbf=2000)
c
      integer inp,iout
      integer nbf,num,lenunt,nnp,maxcor,iadtwp,ngot,need,idum
      integer wptoin,newnbf,nbasis
      character*4096 ops
      character*128 rints,ints
      character*16 bflabl(maxnbf)
      logical logkey
      logical drop, onel, killr
      common /io/ inp,iout
      pointer (p,z(1)), (p,a(1))
      call drum
c
c     --- recover the options string ---
      call iosys('read character options from rwf',-1,0,0,ops)
      call getmem(0,p,ngot,'first',0)
      call iosys('read integer mxcore from rwf',1,maxcor,0,' ')
c
c     --- see if we are to sort the integrals.
c
c              test for one electron only
      onel=logkey(ops,'m330=only-one-electron',.false.,' ')
      killr=logkey(ops,'m330=destroy-raw-integral-file',.false.,' ')
      drop=logkey(ops,'drop',.false.,' ')
      if (onel) then
          write(iout,*) 'm330:'
          write(iout,*) 'm330:only one electron integrals handled'
      endif
      if(logkey(ops,'noints',.false.,' ')
     $  .or.logkey(ops,'int=reuse',.false.,' ')
     $  .or.logkey(ops,'int=reuse2',.false.,' ')
     $  .or.logkey(ops,'int=reuse-m330',.false.,' ')) then
           write(iout,*)'m330: skip integral sort '
         if (drop) then
            call iosys('read character "integral filename" from rwf',
     $                  0,0,0,ints)
            call iosys('open ints as old',0,0,0,ints)
c           --- restore the portion of the read-write pertaining to 
c               functions dropped
            call iosys('read integer "number of basis functions"'
     $               //' from ints',1,nbf,0,' ')
            call iosys('write integer "number of basis functions"'
     $              //' to rwf',1,nbf,0,' ')
            call iosys('read character "basis function labels"'
     $               //' from ints',len(bflabl(1))*nbf,0,0,bflabl)
            call iosys('write character "basis function labels"'
     $              //' to rwf',len(bflabl(1))*nbf,0,0,bflabl)
            call iosys('read integer "old nbf" from ints',
     $                  1,nbf,0,' ') 
            call iosys('write integer "old nbf" to rwf',
     $                  1,nbf,0,' ') 
         endif
      else
c
c       ---- open the integral units ---
         call iosys('read integer "number of basis functions" from rwf',
     $               1,num,0,' ')
         nnp=(num+1)*num/2
         need=wptoin(nnp)
         call iosys('read character "integral filename" from rwf',
     $               0,0,0,ints)
c
         lenunt=nnp*nnp+nnp+20000
         lenunt=min(lenunt,5000000)
         if(logkey(ops,'unit=ssd=ints',.false.,' '))then
            call iosys('open ints as new on ssd',lenunt,0,0,ints)
         else
            call iosys('open ints as new',lenunt,0,0,ints)
         end if
c
c        --- transfer a copy of the one-electron integrals to ints.
         call getmem(need,p,ngot,'m330',0)
         call iosys('copy real "nuclear repulsion energy"'
     $              //' from rwf to ints',1,z(1),0,' ')
         call iosys('copy real "overlap integrals"'
     $              //' from rwf to ints',nnp,z(1),0,' ')
         call iosys('copy real "kinetic integrals"'
     $              //' from rwf to ints',nnp,z(1),0,' ')
         call iosys('copy real "potential integrals"'
     $              //' from rwf to ints',nnp,z(1),0,' ')
c
c        --- call the routines to make the supermatrices ---
c
         if(onel.and.drop) then
            write(iout,*) 'm330:'
            write(iout,*) 'm330:skipping two electron sort'
            call iosys('read integer "number of basis functions" '//
     $                 'from rwf',1,nbasis,0,' ')
            call iosys('read integer "packing index vector" from rwf',
     $                  nbasis,a(1),0,' ')
            call iosys('read integer "truncated number of basis'
     $                 //' functions" from rwf',1,newnbf,0,' ')
            call iosys ('write integer "number of basis functions" '//
     $                  'to rwf',1,newnbf,0,' ') 
            call iosys ('write integer "number of basis functions"'
     $                  //' to ints',1,newnbf,0,' ') 
            call iosys('read character "new basis function labels"'
     $                 //' from rwf',len(bflabl(1))*newnbf,0,0,bflabl)
            call iosys('write character "basis function labels" to rwf',
     $                  len(bflabl(1))*newnbf,0,0,bflabl)
            call iosys('write character "basis function labels" to '//
     $                 'ints',len(bflabl(1))*newnbf,0,0,bflabl)
            call iosys ('write integer "old nbf" to rwf',1,num,0,' ') 
            call iosys ('write integer "old nbf" to ints',1,num,0,' ') 
            call iosys('write integer "truncated number of basis'
     $                 //' functions" to ints',1,newnbf,0,' ')
            call iosys('write integer "packing index vector"'
     $                 //' to ints',nbasis,a(1),0,' ')
            call getmem(-ngot,p,idum,'m330',idum)
         else
            call iosys('read character "raw integral filename" from '//
     $                  'rwf',0,0,0,rints)
            call iosys('open rints as old',0,0,0,rints)
            call getmem(-ngot,p,idum,'m330',idum)
            call mn330(killr)
         endif
c
         call iosys('close ints',0,0,0,' ')
      endif
c
c     --- and exit with grace ---
      call chainx(0)
c
c
      stop
      end
