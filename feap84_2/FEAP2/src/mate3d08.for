      subroutine mate3d08(h1,h2,nh,d,md,Eps,Sig,Cmat,nsig,ntyp,plout,
     +             xgp,tgp,dvp,detf,skfy,skfz,ngp,lgp,lay1gp,lay2gp,isw)
c-----------------------------------------------------------------------
c     
c     Purpose: calculate S and C via FE^2 
c              on a micro problem 
c
c     Inputs:
c         h1(nh)    - history array h1
c         d         - local d-array
c         Eps       - strains  
c         nsig      - dimension of Eps,Sig 
c         dvp       - element volume (resp. length/area for beam/shell elements)
c         isw       - solution option from element 
c
c     Input material parameter:                
c
c     Outputs:
c         sig       - stresses
c         cmat      - tangent modulus
c         plout(10) - plot data    
c
c     Allocation of d-array:
c         no terms  - data are defined in micro problem    
c
c     Usage 
c
c     Inputfile micro-problem 
c      .... 
c      nopr
c      end
c
c      batch
c      nopr
c      loop,,1
c       tang,,1
c      next
c      cmat
c      end
c
c-----------------------------------------------------------------------
      USE iofile
      implicit double precision (a-h,o-z)
      dimension h1(*),h2(*),d(*),eps(nsig),sig(nsig),cmat(nsig,nsig),
     +          plout(10)

      if(isw.eq.1) then
c....   input data
        call mati3d08(md,nh)
c
      else 
c....   calculate C,S via FE^2
        call matm3d08(cmat,eps,sig,nsig,ntyp,dvp,skfy,skfz,ngp,lgp,isw)
c
      end if
      
      return
      end
c
      subroutine mati3d08(md,nh)
c-----------------------------------------------------------------------
c
c     Purpose:  mate3d08 FE^2 micro problem 
c       # input material parameter 
c       # set file names for micro problem
c       # up to 10 RVEs possible
c
c      Input 
c        microfile without _i
c        start user defined batchfile 
c        length of restart-file, check in advance
c
c      Output
c        names of all files used for FE^2
c       
c-----------------------------------------------------------------------
c
      USE fe2mat    ! for matfe2=ma
      USE fe2tran   ! for names FE2 
      USE feapprog  ! for feap   
      USE iofile
      USE pdata2    ! for idev
      implicit double precision (a-h,o-z)
      character*229 fmicro,fbatch,fhelp
      dimension td(2)
      integer :: i
      
      character*229 fbatchmode
      fbatchmode = 'feap_console_no_parallel'
      
c.... input parameter
      md = 0
      nh = 0

c.... input microfile without _i
      if(ior.lt.0)write(*  ,1000) 
1000  format(5x,'Input: Name of micro problem with path')
      read(ior,'(a229)') fmicro
      write(iow,1001) fmicro
1001  format(5x,'FE2:',/,
     +          '     Name of micro problem with path',/,5x,a229)

c.... read and start user defined batchfile 
c       - to delete all no longer used files in micro-directory 
c       - to make copies of inputfile fmicro_01 
c     problem: restart files will be deleted always
      if(ior.lt.0)write(*  ,1002) 
1002  format(5x,'Input: Name of batch file with path')
      read(ior,'(a229)') fbatch
      write(iow,1003) fbatch
1003  format('     Name of batch file with path',/,5x,a229)
      call start_process(fbatch)
c.... example
c     cd d:\w\feap\exe\rve
c     copy irve_01 irve_02
c     copy irve_01 irve_03
c     copy irve_01 irve_04
c     copy irve_01 irve_05
c     copy irve_01 irve_06
c     copy irve_01 irve_07
c     copy irve_01 irve_08
c     ...
c     until nproc=no. of threads(par) or nproc=no. of Gauss points(seq)
c     del frve*
c     del brve*
c     del rrve*
c     del orve*
c     del *.lnk
 
      if(ior.lt.0)write(*  ,1004) ! valid only for one material
1004  format(5x,'I: multi/single restart file, record length')
      call dinput(td,2)
      irtyp = td(1)
      irecl = td(2)
      icgp=0
      write(iow,1005) irtyp,irecl
1005  format('     type of restart-file 0=multi 1=single       ',i10,/,
     +       '     record length restart-file (only for single)',i10,//)

      if(idev.eq.3) irecl=irecl/4 ! INTEL  

c.... set file names for micro problem
      ma=matfe2        
      if(ma.gt.10) stop 'material no. 1-10 necessary in FE^2'

      fcis      = ' ' ! FEAP

      fcisi(ma) = ' ' !-Iifile
      fciso(ma) = ' ' !-Iofile
      fcisr(ma) = ' ' !-Irfile 
      fciss(ma) = ' ' !-Isfile
      fcish(ma) = ' ' !-Isfile for hflgu=.false.
      fcisf(ma) = ' ' ! transfer forward
      fcisb(ma) = ' ' ! transfer backward

      fresg     = ' ' ! overal storage file, valid only for one material

      iex=ipos(  fmicro,229)  

c...  program
      if(idev.eq.3) then ! INTEL
        ip =ipos(feapint,229)
        fcis(1:ip) = feapint(1:ip)
      else if(idev.eq.4) then ! SALFORD
        ip =ipos(feapsal,229)
        fcis(1:ip) = feapsal(1:ip)
      end if
#ifdef _BATCH_
        ip =ipos(fbatchmode,229)
        do i=1,229
            fcis(1:i) = ' '
        end do
        fcis(1:ip) = fbatchmode(1:ip)
#endif
c...  -Iifile 
      fhelp = ' '
      fhelp(1:3) = ' -i'               
      ip = 4 
      fhelp(ip:ip+iex) = fmicro(1:iex) 
      fcisi(ma) = fhelp

c...  -Oofile  
      fhelp = ' '
      fhelp = fcisi(ma)
      call dochar2(fhelp,ip)     ! look for i
      call dochar1(fhelp,'o',ip) ! set i -> o
      fhelp(1:3) = ' -o' 
      fciso(ma) = fhelp

c...  -Rrfile (read)
      fhelp = ' '
      fhelp = fcisi(ma)
      call dochar1(fhelp,'r',ip) ! set i -> r
      fhelp(1:3) = ' -r' 
      fcisr(ma) = fhelp

c...  -Srfile (save)
      fhelp = ' '
      fhelp = fcisr(ma)
      fhelp(1:3) = ' -s'      
      fciss(ma) = fhelp

c...  -Srfile (save) for hflgu = false.
      fhelp = ' '
      fhelp = fciss(ma)
      call dochar1(fhelp,'h',ip) ! set r -> h
      fcish = fhelp
 
c...  Transfer files fcisf=E, ficisb=C,S
      fhelp = ' '
      fhelp(1:iex) = fmicro(1:iex) 
      call dochar2(fhelp,ip)     ! look for i
      call dochar1(fhelp,'f',ip) ! set i -> f
      fcisf(ma) = fhelp 

      fhelp = ' '
      fhelp(1:iex) = fmicro(1:iex) 
      call dochar2(fhelp,ip)     ! look for i
      call dochar1(fhelp,'b',ip) ! set i -> b
      fcisb(ma) = fhelp 

c...  FRESG global Restart File, includes restart data at all (ngp,lgp) 
c     valid only for one material
      fresg(1:iex) = fmicro(1:iex) 
      call dochar2(fresg,ip)     ! look for i
      call dochar1(fresg,'r',ip) ! set i -> r

      return
      end 
c
      subroutine matm3d08(cmat,e,sig,nsig,ntyp,dvp,skfy,skfz,ngp,lgp,isw
     &)
c-----------------------------------------------------------------------
c
c     Purpose: mate3d08 FE^2 micro problem  
c
c     Inputs:
c      e(nsig)         - strain vector 
c      nsig            - =6(3D) =8(Shell) =10(Shell ext.) =6(Beam) =15(Beam ext.)
c      ntyp            - 1=3D, 2=shell, 3=shell ext. 4=beam 5=beam ext.
c      dvp             - element length/area/volume
c      skfy            - sqrt of shear correct. factor k_y 
c      skfz            - sqrt of shear correct. factor k_z
c      ngp             - n=element number
c      lgp             - l=GP      number
c      isw             - Solution switch 
c
c     Outputs:
c      cmat(nsig,nsig) - material matrix C
c      sig(nsig)       - stress vector   S
c
c     Comments:
c     up to 10 RVEs possible
c     ifile is input in mati3d08
c
c     micro problem is started via L=GP N=EL
c     feap.exe -Iifile_L -Oofile_L -Rrfile_L   -Srfile_L 
c     feap.exe -Iifile_L -Oofile_L -Rrfile_L_N -Srfile_L_N 
c
c     for hflgu=false no save of history variables
c     feap.exe -Iifile_L -Oofile_L -Rrfile_L   -Shfile_L
c     feap.exe -Iifile_L -Oofile_L -Rrfile_L_N -Shfile_L
c
c     Transferfiles
c_INT ife2f = 100+threadId
c_INT ife2b = 200+threadId
c_SAL ife2f = 100+lgp
c_SAL ife2b = 200+lgp
c
c
c     Transfer ibin=0 with Fresg=rfile and fresl=-rrfile_L_N
c
c     irtyp = 0 all rfiles are separate rfile_L_N
c     irtyp = 1 all rfiles are in File FRESG  
c
c     Transfer in microproblem EPSQ/SIGQ via 33/34
c
c     ibin=0  Transfer von  ngp,lgp via ife2f/b 
c     ibin=1,2: fehlt
c
c
c-----------------------------------------------------------------------
#ifndef _LINUX_
      USE KERNEL32
#endif
      USE fdata
      USE fe2mat    ! for matfe2=ma
      USE fe2tran   ! for names FE2
      USE hdatam    ! hfglu
      USE iofile
      USE pdata2    ! for idev
      USE tdata
      use memorymap
      implicit double precision (a-h,o-z)
      dimension cmat(nsig,nsig),e(nsig),sig(nsig)
      character*229 fcisi2,fciso2,fcisr2,fciss2,fcisf2,fcisb2,fresl
      character*600 fcis1 
      character cgp*2 
      character cel*4
      character ciswb*4
      
#ifdef __INTEL_
      INTEGER :: OMP_GET_THREAD_NUM
      INTEGER*8 :: emap, csmap, nssmap, rsmap, threadId
      CHARACTER(len=2) :: cThreadId
      CHARACTER(len=24) :: rsmapname
#endif

      fcis1  = ' '
      fcisi2 = ' '
      fciso2 = ' '
      fcisr2 = ' '
      fciss2 = ' '
      fcisf2 = ' '
      fcisb2 = ' '
      fresl  = ' '

c...  transfer version
      if(idev.eq.3) ibin=0 ! INTEL  
      if(idev.eq.4) ibin=0 ! SALFORD
      
c...  iswb number
      write(ciswb,'(a2,i2)') '-b',isw
      if(ciswb(3:3) .eq. ' ') ciswb(3:3) = '0' 

c...  GP number
      write(cgp,'(i2.2)') lgp

c.... Element number 0/1 without/with element number, see fe2tran 
cDH.. Do this for irtyp=2 too!
      if(irtyp.eq.0.or.irtyp.eq.2) then
        write(cel,'(i4.4)') ngp
      end if  
      
c.... actual material number
      ma=matfe2

#ifdef __INTEL_
c...  use thread ID for output filenames instead of GP number
c...  to avoid conflicts in parallel version
      threadId = OMP_GET_THREAD_NUM() + 1
      write(cThreadId,'(i2.2)') threadId
#endif
      
c.... -Iifile_l 
      fcisi2 = fcisi(ma)
c.... add GP number
      ia=ipos(fcisi2,229)  
      fcisi2(ia+1:ia+1)='_'  
#ifdef __INTEL_
      fcisi2(ia+2:ia+3)=cThreadId(1:2)
#endif
#ifndef __INTEL_
      fcisi2(ia+2:ia+3)=cgp(1:2)
#endif
      
c...  -Oofile_l
      fciso2 = fciso(ma)
c.... add GP number
      ia=ipos(fciso2,229)  
      fciso2(ia+1:ia+1)='_'  
#ifdef __INTEL_
      fciso2(ia+2:ia+3)=cThreadId(1:2)
#endif
#ifndef __INTEL_
      fciso2(ia+2:ia+3)=cgp(1:2)   
#endif

c...  -Rrfile_l<.n> (read)
      fcisr2 = fcisr(ma)
      
c.... add GP number
      ia=ipos(fcisr2,229)  
      fcisr2(ia+1:ia+1)='_'  
      fcisr2(ia+2:ia+3)=cgp(1:2)   
c...  add element number
      if(irtyp.eq.0) call addext(fcisr2,cel)

#ifdef __INTEL_
c     determine name for restart file mapping
      if(irtyp.eq.2) then
c     add element number here, too
        call addext(fcisr2,cel)
        call feap_shmem_ParseRestartMappingName(fcisr2,rsmapname)
      end if
#endif

      if(hflgu.and.h3flgu) then 
c...    -Ssfile_l<.n> (save)
        fciss2 = fciss(ma)
c....   add GP number
        ia=ipos(fciss2,229)  
        fciss2(ia+1:ia+1)='_'  
        fciss2(ia+2:ia+3)=cgp(1:2)   
c...    add element number
cDH   also for irtyp=2
        if(irtyp.eq.0.or.irtyp.eq.2) call addext(fciss2,cel)
      else  
c...    -Shfile_L for hflgu=.false.
        fciss2 = fcish(ma)
c.... add GP number
        ia=ipos(fciss2,229)  
        fciss2(ia+1:ia+1)='_'  
        fciss2(ia+2:ia+3)=cgp(1:2)   
      end if 
              
c.... -Tffile_l bfile_l
      fcisf2 = fcisf(ma)
      fcisb2 = fcisb(ma)
c.... add GP number
      ia=ipos(fcisf2,229)

#ifdef __INTEL_
      fcisf2(ia+1:ia+1)='_'
      fcisf2(ia+2:ia+3)=cThreadId(1:2)   
      fcisb2(ia+1:ia+1)='_'
      fcisb2(ia+2:ia+3)=cThreadId(1:2)
#endif      
#ifndef __INTEL_
      fcisf2(ia+1:ia+1)='_'  
      fcisf2(ia+2:ia+3)=cgp(1:2)   
      fcisb2(ia+1:ia+1)='_'  
      fcisb2(ia+2:ia+3)=cgp(1:2)
#endif
      
c.... store files into fcis1
      ii=ipos(fcisi2,229)  
      io=ipos(fciso2,229)  
      ir=ipos(fcisr2,229)  
      is=ipos(fciss2,229)  

      ia=1
      ie=ii
      fcis1(ia:ie)=fcisi2(1:ii)  !Ifile
      ia=ie+1
      ie=ia+io 
      fcis1(ia:ie)=fciso2(1:io)  !Ofile
      ia=ie+1
      ie=ia+ir 
      fcis1(ia:ie)=fcisr2(1:ir)  !Rfile
      ia=ie+1
      ie=ia+is 
      fcis1(ia:ie)=fciss2(1:is)  !Sfile
      ia=ie+1
      ie=ia+4 
      fcis1(ia:ie)=ciswb(1:4)     !iswb 
      
c.... check length
      if(ie.gt.600) then 
        write(*,*) 'mate3d08 fcis1 need',ie,'characters' 
        stop 'allowed are 600 characters'
      end if

c.... transfer channels, also ofiles 
c     parallel: thread-number, on one thread is 1 element with all GPs 
c             seq. element/Thread could be seen in write statement below
c     seriell: GP-number, element with seq. GPs seq. 
#ifdef __INTEL_      
      ifadd=threadId
#endif
#ifndef __INTEL_
      ifadd=lgp
#endif
      ife2f=100+ifadd
      ife2b=200+ifadd
c.... control EL-Number,ISW,Material,Thread-Number(INTEL) or GP-Number(FTN95)
c     do below: afterwards  
c      if(pfr.and.lgp.eq.1) then 
c        write(*,'(a16,i6,i6,i6,i6)') ' EL,ISW,MA,Th/GP',ngp,isw,ma,ifadd
c      end if 

c.... send E to micro problem
      if(isw.eq.3.or.isw.eq.4.or.isw.eq.6.or.isw.eq.8)then  
        if(ibin.eq.0) then

c....     extract restart files from fresg, see elmt21s 

c....     write data to transfer file
          open(ife2f,file=fcisf2,form='unformatted',status='unknown')
          rewind(ife2f)
          write(ife2f) dt,ttim,dvp,skfy,skfz,isw,nsig,ntyp
          write(ife2f) e
          write(ife2f) ngp,lgp
          close(ife2f)

        else if (ibin.eq.2) then

#ifdef __INTEL_
c     Create file mappings for E,nss and C,S
          call feap_shmem_CreateMappings(emap, csmap, nssmap, nsig,lgp,
     &   ngp)
                        
c     write E,nss to mapping
          call feap_shmem_writeE(emap,e,dt,ttim,dvp,skfy,skfz,isw,nsig,
     &ntyp)
          call feap_shmem_writeNSS(nssmap,nsig)
          
#endif
        else 
          open(ife2f,file=fcisf2,form='formatted',status='unknown')
          rewind(ife2f)
          write(ife2f,'(5e18.9,3i2)') dt,ttim,dvp,skfy,skfz,isw,nsig,
     &ntyp
          write(ife2f,'(15e18.9)') (e(i),i=1,nsig)
          close(ife2f)

        end if ! ibin
      end if ! isw 
      
#ifdef __INTEL_
      if (irtyp.eq.2) then
c...  If shared memory is used for restart files and the mapping has not
c...  yet been opened, open it now.
#ifdef _LINUX_
          rsmap = shm_open(rsmapname, oread, fperm)
#else
          rsmap = OpenFileMapping( FILE_MAP_ALL_ACCESS, FALSE,
     &                              rsmapname )
#endif
          if(rsmap.eq.0) then
c     mapping was not yet created, do that now
           call feap_shmem_CreateRestartMapping(rsmap,rsmapname,irecl)
          end if
      end if
#endif

c.... start micro problem  
c     feap.exe -Iifile_L -Oofile_L -Rrfile_L -Srfile_L(or -shfile) -bisw
#ifdef __SHMDEBUG_
      write (*,*) "call: ", fcis, fcis1
#endif
      call start_process_fe2(fcis,fcis1)
      
c.... Read C and S from micro problem
      if(isw.eq.3.or.isw.eq.4.or.isw.eq.6.or.isw.eq.8)then  

        if(ibin.eq.0) then
       
c....     read data from transfer file
          open(ife2b,file=fcisb2,form='unformatted',status='unknown')
          rewind(ife2b)
          read(ife2b) sig
          read(ife2b) cmat
          read(ife2b) it1,it2,resife2 ! it-behavior micro-problem  
          close(ife2b)
        

c....     control EL-Number,ISW,Material,Thread-No(INTEL)/GP-No(FTN95),
c                 local its local Resi

          if(pfr.and.lgp.eq.1) then 
            write(*,'(a24,i6,i6,i6,i6,i6,e16.5)') 
     +      ' EL,ISW,MA,Th/GP,IT,|R| ',ngp,isw,ma,ifadd,it1,resife2
          end if 

c....     show it-behavior micro-problem of element at GP/TH
c         if no local convergence (except it2=2)
c         activate, if necessary!   
c          if(pfr.and.lgp.eq.1) then 
            epsresife2=1.e-5 ! set error bound for RVE residual
            if(resife2.gt.epsresife2.and.it2.ne.3) then  ! for 3 local iterat
              write(  *,'(a16,i6,i6,i6,e16.5)') ' EL,it1,it2,Resi',
     +              ngp,it1,it2,resife2
              write(iow,'(a16,i6,i6,i6,e16.5)') ' EL,it1,it2,Resi',
     +              ngp,it1,it2,resife2
            end if
c          end if
          
c....     Save restart files in FRESG, see elmt21s 

        else if(ibin.eq.2) then

#ifdef __INTEL_
        call feap_shmem_readCS(csmap,sig,cmat,nsig)
          
c       close all file mappings
        call feap_shmem_CloseHandle(emap)
        call feap_shmem_CloseHandle(csmap)
        call feap_shmem_CloseHandle(nssmap)
#endif

        else

          open(ife2b,file=fcisb2,form='formatted',status='unknown')
          rewind(ife2b)
          read(ife2b,'(15e18.9)') (sig(i),i=1,nsig) 
          do i = 1,nsig
            read(ife2b,'(15e18.9)') (cmat(i,k),k=1,nsig) 
          end do      
          close(ife2b)

        end if ! ibin
      end if ! isw

      return
      end
c
      subroutine matt3d08(isw,n,numel,mgp)
c----------------------------------------------------------------------
c
c      Purpose: Read:  Copy    file FRESL into in overall file FRESG
c               Write: Extract file FRESL from    overall file FRESG
c
c               1)  isw=1 FRESG>FRESL  
c               2)        Solv Microproblem with FRESL
c               3   isw=2 FRESL>FRESG  
c
c      Inputs:
c         FRESG   - FRESG(numel,mgp,irecl) contains all restart files 
c                   of MICRO-problems
c
c         FRESL   - Name of the restart file of the MICRO-problem at 
c                   element n and gauss-point lgp of the MACRO-problem   
c         n       - actual element number
c         numel   - maximum element number
c         mgp     - maximum gauss point number
c         irecl   - length of data in fresl 
c         icgp    - counter for 1. TANG up to mgp*numel
c
c         iog     - channel FRESG
c         iol     - channel FRESL
c
c      Outputs:
c
c      Comments:       
c
c      Routine has to be implemented into element before(isw=1) and 
c      after (isw=2) call to mate3d08, see ELMT21 
c
c      Restart data are
c        # general values 
c        # displacements
c        # arc length values
c        # history fields  
c
c        # <transient fields v,a>  not included!! 
c        # edge data
c        # <crack values>
c        # <director values>
c
c        RECL in Bytes
c        nnpo,nnlo,nnmo,ndmo,ndfo,nrt,fl(9)  ! 6*int*4+1Log*4  = 28
c        ttim,(b(i),i=1,nneq*iu)             !(1+nneq*3)real*8 = 8+24*nneq 
c        prop,rlnew,c0,cs1,cs2,ds0,r,det0,xn ! 9real*8         = 72
c        ttim,gh1,gh2,gh3                    ! 1real*8+real*8(2*nhmax+nh3max)
c        clfl                                ! 1log*4          = 4
c        RECL= 120+24*nneq+*8(2*nhmax+nh3max)
c
c----------------------------------------------------------------------
      USE fe2tran
      USE pdata2     ! for idev
      implicit double precision (a-h,o-z)
      logical exstg,exstl
      character fresl*229
      character cgp*2 
      dimension mxi(irecl)
      dimension mxs(1)

      irecl1=irecl
      if(idev.eq.3) irecl1=irecl*4
      
      iog = 35
      iol = 36 

      fresl  = ' '

c.... open file FRESG 
      inquire(file=fresg,exist=exstg)
      if(exstg) then ! 
        open(iog,file=fresg,form='unformatted',status='old',err=2,
     +       access='direct',recl=irecl1)
      else ! 
        open(iog,file=fresg,form='unformatted',status='new',err=2,
     +       access='direct',recl=irecl1)
      end if

c...  set name of FRESL
      fresl = fresg
      ia=ipos(fresl,229)  
      fresl(ia+1:ia+1)='_'  

c...  loop over mgp Gauss-Points of actual element
      do 100 lgp =1,mgp    

c...    set name of FRESL
        write(cgp,'(i2)') lgp
        if(cgp(1:1) .eq. ' ') cgp(1:1) = '0' 
        fresl(ia+2:ia+3)=cgp(1:2)   

c....   record to open 
        ir = (n-1)*mgp+lgp

c....   open file FRESL 
        inquire(file=fresl,exist=exstl)
        if(exstl) then 
          open(iol,file=fresl,form='unformatted',status='old',err=1,
     +       access='direct',recl=irecl1)
        else ! on mgp GaussPoints of first element  
c         no action    
        end if

        if(isw.eq.1) then
c....     read general file and copy to local file 
          icgp=icgp+1 
          if(icgp.le.mgp*numel) then 
cww         write(*,*) 'icgp1',icgp
            close(iol,status='delete') 
          else
cww         write(*,*) 'G>L ir',ir
            if(idev.eq.3) then
              read (iog,rec=ir,err=3) mxi 
              write(iol,rec=1, err=4) mxi 
            else if(idev.eq.4) then
              read (iog,rec=ir,err=3) mxs 
              write(iol,rec=1, err=4) mxs 
            end if
          end if
        else if(isw.eq.2) then
c....     read local file and copy to general file 
cww       write(*,*) 'L>G ir',ir        
          if(idev.eq.3) then
            read (iol,rec=1, err=5) mxi 
            write(iog,rec=ir,err=6) mxi 
          else if(idev.eq.4) then
            read (iol,rec=1, err=5) mxs 
            write(iog,rec=ir,err=6) mxs 
          end if
        end if
        close(iol)    
100   continue

      close(iog)
      return
      
1     write(*,1001)n,lgp
      return
1001  format(' Error open file FRESL! N= ',i5,' l= ',i2)

2     write(*,1002)n,lgp
      return
1002  format(' Error open file FRESG! N= ',i5,' l= ',i2)

3     continue
      write(*,1003) n,lgp
      stop
1003  format(' Error Copy FRESG>FRESL on  READ! N= ',i5,' l= ',i2)

4     write(*,1004) n,lgp
      stop
1004  format(' Error Copy FRESG>FRESL on WRITE! N= ',i5,' l= ',i2)

5     write(*,1005) n,lgp
      stop
1005  format(' Error Copy FRESL>FRESG on  READ! N= ',i5,' l= ',i2)

6     write(*,1006) n,lgp
      stop
1006  format(' Error Copy FRESL>FRESG on WRITE! N= ',i5,' l= ',i2)
      end   
      
