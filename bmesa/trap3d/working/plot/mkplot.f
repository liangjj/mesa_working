      subroutine mkplot(z,ia)
c***begin prologue     mkplot
c***date written       920324   (yymmdd)
c***revision date               (yymmdd)
c
c
c***keywords
c***author             
c***source             %W% %G% 
c***purpose
c                       
c***description
c***references
c
c***routines called
c
c***end prologue       mkplot
c
      implicit integer (a-z)
      character*4096 ops
      character*128 filbec, becin
      character*8 cpass
      character*3 ans
      character*240 card, chrkey
      logical logkey, prnt
      real*8 z, natoms, energy, scale, length
      dimension z(*), ia(*)
      common/io/inp,iout
      call iosys ('read character options from rwf',-1,0,0,ops)
      call iosys ('read character "bec filename" from rwf',
     1             -1,0,0,filbec)
      call iosys ('open bec as unknown',0,0,0,filbec)
      call iosys('rewind all on bec',0,0,0,' ')
      call iosys('read integer "size of basis" from bec',1,
     1            n,0,' ')
      call iosys('read integer maxd from bec',1,maxd,0,' ')
      call iosys('read real scale from bec',1,scale,0,' ')
      call posinp('$plot',cpass)
      call cardin(card)
      prnt=logkey(card,'print-only',.false.,' ')
      becin=chrkey(card,'input-plot-filename','"in plot file"',' ')
      call iosys('does '//becin//' exist on bec',0,0,0,ans)
      write(iout,4) becin, ans
      if(ans.eq.'no') then
         call lnkerr('file does not exist on bec. quit')
      endif         
      nfiles=intkey(card,'number-of-files-to-read',100,' ')
      nvals=intkey(card,'number-of-values-to-plot',10,' ')
      if(.not.prnt) then
          psi=1
          vnl=psi+n
          eig=vnl+n
          wt=eig+3*maxd
          pltarr=wpadti(wt+3*maxd)
          atarr=iadtwp(pltarr+nvals)
          call intarr(card,'plot-array',ia(pltarr),nvals,' ')
          write(iout,1) nvals
          call wfn(z(psi),z(vnl),z(eig),z(wt),z(atarr),ia(pltarr),
     1             n,nvals,maxd)
      else
          call iosys('read real '//becin//' from bec',1,length,0,' ')     
          files=length
          count=0
          offset=1
          write(iout,3)
          do while ( count.lt.files)
             count = count + 1
             call iosys('read real '//becin//' from bec',1,
     1                   natoms,offset,' ')
             offset = offset +1
             call iosys('read real '//becin//' from bec',1,energy,
     1                   offset,' ')
             offset = offset + 1
             write(iout,5) natoms, energy 
             offset = offset + n + n
          enddo     
      endif       
      call chainx(0)
c
c
 1    format(/,5x,'number of values to plot = ',i3)
 2    format(/,5x,'the length of the plot file = ',i10,//)
 3    format(15x,'number of atoms',8x,'energy')
 4    format(/,10x,'filename =',a60,/,10x,'answer = ',a3,//)
 5    format(/,13x,e15.8,5x,e15.8)
      stop
      end
