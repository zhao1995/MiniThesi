c      program feap
c----------------------------------------------------------------------
c
c.... Finite Element Analysis Program  (FEAP)  for solution of general
c.... problem classes using the finite element method. Problem size
c.... is controlled by the dimension of blank common and value of maxm
c.... as set in module psize in module.for.
c
c.... Precision is controlled by ipr (=1 for single,  =2  for double).
c.... If double precision include: implicit double precision (a-h,o-z)
c.... in each subprogram using double precision quantities.
c
c----------------------------------------------------------------------
c.... Programmed by:
c....            R. L. Taylor
c....            Department of Civil Engineering
c....            University of California
c....            Berkeley, California 94720
c
c----------------------------------------------------------------------
c.... Programmed by:
c....            W. Wagner
c....            Institut fuer Baustatik
c....            University of Karlsruhe
c
c.... Copyright (C) 1984 - 2014
c----------------------------------------------------------------------
c
c.... Input/Output is performed to files during execution of FEAP.
c.... In general the following files are used during executions:
c
c              8 :  finp.tau  1=t! output taurus
c       iop = 11 :  fplt      Used plot outputs
c       ios = 12 :            Used for read/write scratch files
c       ird = 13 :            Used to read results data from disk.
c       iwd = 14 :            Used to write results data to disk.
c       ior = 15 :  finp      Use to read from the input data file
c                             specified when a problem is initiated.
c       iow = 16 :  fout      Use to write ouput result data to file
c                             specified when a problem is initiated.
c             17 :  fplt_x    cgm/ps output under gks
c             18 :  wiss      output under gks
c             19 :  procedure
c             20 :  procedure save
c             22 :  fres.err  output error data
c             23 : -fres.pos  pos output under post macro (data)
c                  -          section
c                  -fres.tec  tec-data
c                  -fres.mxx.tyyyy.vtu  paraview-data
c              
c             25 :  fres.ldd  file for save node defined values 
c             30 :  fplt.pgl  plot output under prin macro
c             31 :  fplt.eps  plot output under prin macro
c     iocfg = 32 :  cfg-file 
c     iom1  = 33 :  ffile     get  E   from macro-problem  for EPSQ/SIGQ
c     iom2  = 34 :  bfile     send C,S to   macro-problem  for EPSQ/SIGQ
c       iog = 35 :  fresg     Rstart global file
c       iol = 36 :  fresl     Rstart local file
c
c       iow = 51 :  foutpar_01   output-processor 1=0     
c       iow = 52 :  foutpar_02   output-processor 2=1     
c       iow = 53 :  foutpar_03   output-processor 3=2     
c       iow = 54 :  foutpar_04   output-processor 4=3     
c       iow = 55 :  foutpar_05   output-processor 5=4    
c       iow = 56 :  foutpar_06   output-processor 6=5    
c       iow = 57 :  foutpar_07   output-processor 7=6    
c       iow = 58 :  foutpar_08   output-processor 8=7    
c
c       ife2f1 = 101 : send E   to   micro-problem for MATE3D08 parallel      
c       ...                     for each GP (e.g.8) 
c       ife2f8 = 108 : send E   to   micro-problem for MATE3D08 parallel      
c       ife2b1 = 201 : get  C,S from micro-problem for MATE3D08 parallel     
c       ...                     for each GP (e.g.8)
c       ife2b1 = 208 : get  C,S from micro-problem for MATE3D08 parallel     
c
c     FILE names
c     finp        Input 
c     finp.lnk    Link  
c     finp.tmp    Adaptivity
c     finp.mrm    Adaptivity
c     
c     fout        Ouput 
c     
c     fres        Rstart
c     fres.name   save on special file
c     fres.ijk    save on numbered files
c     fres.err    Error data
c     fres.ldf    LDF-curves
c     fres.ldd    LDF-curves more results 
c     fres.pos    FEPOST    
c     fres.doc    FEPOST    
c     fres.tec    TECPLOT   
c     fres.mxx.tyyyy.vtu PARVIEW Data   
c     fres.pvd    PARVIEW define all files of a time step     
c     
c     fsav        Rstart(Save)
c     
c     fout        Ouput 
c     
c     fplt        Plot  
c     fplt_x.eps  EPS-file under PLOF/PRIN
c     fplt_x.pgl  PGL-file under PRIN
c     fplt_x.pcx  PCX-file under PRIN
c     fplt_x.cgm  CGM-file under PLOF     
c     
C
c----------------------------------------------------------------------
      USE cdata
cww      USE comfil
      USE iodata
      USE iofile
      USE pdata8
      USE plong
      USE printh
      USE psize
      USE rfeap
      USE soltyp
      USE vdata
      implicit real*8 (a-h,o-z)

c-----set memory capacity-----------------------------------------------
      noff   = 1
      kmax   = 1
c-----------------------------------------------------------------------

c.... set version header 
      versn(2) = 'July, 14,2015'
      versn(3) = 'Baustatik - KIT '
c
c.... set precision
      ipr      = 2

cwwc.... set machne:  1-unix/VMS, 2-CMS
cww      machne   = 1
c
c.... reserve memory size limits size of problems which may be solved
cww      call  memlim (machne,ipr)
c

c.... set default file logical unit numbers
      iop = 11
      ios = 12
      ird = 13
      iwd = 14
      ior = 15
      iow = 16

#ifndef _LINUX_
c.... read Config-file, up to now only valid for windows!
      call feapcfg
#endif

c.... open files and set device
      call filnam(0)

c.... check for nege/gmesh/cylt
      call test_gen(1)

      if(iplot.gt.0) then 
c....   open windows for console and graphics
        call startgr(1)
c....   plot logo  
        call pofeap
      end if

c.... set start time
      call tstart

c.... start execution
10    irfeap = 0 
      call pcontr

c.... restart internal RINP RINP,NEW REME,NEW REME,OLD
      if(irfeap.gt.0) then

c....   Reset variable for IGA/FEA distinction
        call resetNurbs()           

c....   clear pardiso memory            
        if(istyp.eq.4.or.istyp.eq.8) call clear12(neq) 


      if(irfeap.eq.3) then
        irfeap=0
        goto 10  ! REME,NEW REME,OLD  
      end if

c....   'normal'  not REME     
        ior=iabs(ior)
cww     rewind(ior)  ! in pcontr before inte 
cww     close(ior)

c....   reset output-file    
        rewind(iow)
        close(iow) 
        close(25) ! ldd-file ww 

c....   RINP,NEW new file to select 
        if(irfeap.eq.1) call filnam(1)

c....   check input file for NEGE, GMESH is used, not possible for REME   
        call test_gen(2)

c....   reset
        irfeap=0

c....   restart
        go to 10

      end if

c.... calculate used memory (if prin in INPUT) 
      if(prt) call length

c.... close i/o and graphics windows 
      if(iplot.gt.0) call startgr(2)

c.... end program
c...  skip exit window for INTEL
      call fexit()
      stop
      end
c
c----------------------------------------------------------------------
c
      subroutine filnam(irfeap)
c----------------------------------------------------------------------
c
c.....Purpose: general file handling control for FEAP  
c              set type of compiler and graphcial device 
c
c     Inputs: interactive (INTEL/SALFORD) or batch(all) 
c             irfeap=0 not rinp
c             irfeap=1 rinp,new
c             ifreap=2=rinp 
c
c     Outputs: File feapname with 4 entries for inp/outp/rest/plot
c
c
c     W. Wagner BS UKA 04/09
c----------------------------------------------------------------------

c.... Declare variable types
      USE bisw
      USE comfil
      USE errchk
      USE feapprog
      USE inptc
      USE iodata
      USE iofile
      USE pdata2
      USE pdata8
      USE psize
      USE vdata
      logical*4 resultl,initialsettings
      logical pcomp,lfil,linp,lout,lres,lsav,lplt
      character*6  wd(2),y*1
      character    sctyp(4)*13
      character*229 dinp,dout,dres,dplt,blk
      character*229 pinp,pout,pres,psav,pplt,filein
      character*229 iname
      character*229 feapname
      integer iorsav, nargs, ierr, iend,irfeap
      integer iinp, ioup, ires, isav, iplt
      integer idevin

c.... Intrinsics
      intrinsic max, min
c
      data dinp,dout,dres,dplt/'NONE','NONE','NONE','PLOT'/,
     1     blk/' '/, wd/'new','exists'/,iname/'./i*'/
      data sctyp/'1=IBM-PHIGS  ','2=HP-GKS     ','3=WIN-INTEL  ',
     1 '4=WIN-SALFORD'/
      data lfil,linp,lout,lres,lsav,lplt/.false.,.false.,.false.,
     +                                   .false.,.false.,.false./

c     show state of input, do if numel/numnp>1000 or NEGE/GMESH (Windows)
      inptctrl = 0 
      numnpic  = 1000
      numelic  = 1000
c
      iorsav = ior
      nargs  = 0
      iplot  = 0 ! no plot modus for batch

c.... set type of compiler and graphcial device 
      call set_device(idev)

c.... batch executions-read command line and set files.   
      call doargs(finp,fout,fres,fsav,fplt,idev,iswb,nargs)

c     if compiled batch set idev to 1 to bypass window settings
#ifdef _BATCH_
      idev = 1
#endif

c.... windows interactive input via mask
      if((idev.eq.3.or.idev.eq.4) .and. nargs.eq.0) then
         iplot  = 1
c....    set frame and open dialog window (only INTEL) 
         if(irfeap.ne.1) resultl = initialsettings( )
         call pwopn
         
c....    file names
         call filnamw
         return
      end if

c set back to 3 (intel)
#ifdef _BATCH_
      idev = 3
#endif

      if(nargs.gt.0) goto 200
c....   generell interactive input 
        iplot  = 1
        ior    =-iabs(ior)
c....   write head to screen
        write(*,2000) versn,maxm
c....   look to see if any problem has been run
        feapname= './feapname'
        inquire(file=feapname,exist=lfil)
        if(lfil) then
          open(ios,file=feapname,status='unknown')
          read(ios,1003) pinp,pout,pres,psav,pplt
          close(ios)
          finp = pinp
          fout = pout
          fres = pres
          fsav = psav
          fplt = pplt
          go to 200
        else
          pinp = dinp
          pout = dout
          pres = dres
          psav = dres
          pplt = dplt
        end if

c....   name file for input data
100     assign  1 to ierr
        assign 10 to iend
        finp = pinp
1       write(*,2001) pinp
        read (*,1000,err=901,end=902) filein
        if (filein.ne.blk)  finp = filein

c....   check if the input files exists
10      inquire(file=finp,exist=linp)
        if(.not.linp) then
          call prin_inp(iname,filein)
          go to 1
        end if
        pinp = finp

c....   set default files for a filname beginning with 'I'
        call dochar2(pinp,iposc)
        if(pcomp(pinp(iposc:iposc),'I',1)) then
          pout = pinp
          pres = pinp
          psav = pinp
          pplt = pinp
          call dochar1(pout,'O',iposc)
          call dochar1(pres,'R',iposc)
          call dochar1(psav,'R',iposc)
          call dochar1(pplt,'P',iposc)
        end if

c....   name file for output data
        assign  2 to ierr
        assign 20 to iend
        fout = pout
2       write(*,2002) pout
        read (*,1000,err=901,end=902) filein
        if (filein.ne.blk)  fout = filein
20      pout = fout

c....   name file for restart read data
        assign  3 to ierr
        assign 30 to iend
        fres = pres
3       write(*,2003) pres
        read (*,1000,err=901,end=902) filein
        if (filein.ne.blk)  fres = filein
30      pres = fres

c....   name file for restart save data
        assign  4 to ierr
        assign 40 to iend
        fsav = psav
4       write(*,2004) psav
        read (*,1000,err=901,end=902) filein
        if (filein.ne.blk)  fsav = filein
40      psav = fsav

c....   name file for plot data
        assign  5 to ierr
        assign 50 to iend
        fplt = pplt
5       write(*,2005) pplt
        read (*,1000,err=901,end=902) filein
        if (filein.ne.blk)  then
          fplt = filein
        end if
50      pplt = fplt

c....   specify graphic device
        assign  6 to ierr
        assign 60 to iend
6       write(*,2006) sctyp(idev)
        read (*,1002,err=901,end=902) idevin
        if (idevin.ne.0) idev = idevin
60      if (idev.lt.1.or.idev.gt.4)  goto 6 

 
c.... check file status and input if necessary
200   inquire(file=finp,exist=linp)
      if(.not.linp.and.nargs.gt.0) stop 'no input specified'
      if(.not.linp) goto 100

      iinp = 2
      inquire(file=fout,exist=lout)
      ioup = 1
      if(lout) ioup = 2
      inquire(file=fres,exist=lres)
      ires = 1
      if(lres) ires = 2
      inquire(file=fsav,exist=lsav)
      isav = 1
      if(lsav) isav = 2
      inquire(file=fplt,exist=lplt)
      iplt = 1
      if(lplt) iplt = 2

      if(nargs.gt.0) goto 300
        write(*,2007) finp,wd(iinp),fout,wd(ioup),fres,wd(ires)
        write(*,2008) fsav,wd(isav),fplt,wd(iplt),sctyp(idev),idev
        write(*,2009)
        
        assign  7 to ierr
7       read(*,1001,err=901,end=901) y
        if(y.eq.'S' .or. y.eq.'s') then
          call startgr(2)
          stop
        end if
        if(y.ne.'Y' .and. y.ne.'y'.and.y.ne.' ') go to 100
c....   save a copy of the current filenames
        open(ios,file=feapname,status='unknown')
        rewind ios
        write(ios,1003) finp,fout,fres,fsav,fplt
        close(ios)
        write(*,2010)
 
c.... erase the output file if it exists
300   ior = iorsav
      if(lout) then
          open(iow,file=fout,status='unknown')
        close(iow,status='delete')
      end if
      
c.... open the files for input and output
      iorww = iabs(ior)
      open(unit=iorww,file=finp,status='unknown')
      open(unit=  iow,file=fout,status='unknown')
      if(nargs.gt.0) then
        write(iow,2007) finp,wd(iinp),fout,wd(ioup),fres,wd(ires)
        write(iow,2008) fsav,wd(isav),fplt,wd(iplt),sctyp(idev),idev
      end if
      return
c.... error trap
901   write(*,3001)
      goto  ierr
c.... eof encountered
902   call  endclr ('FILNAM',filein)
      goto  iend
c.... format statements
1000  format(a229)
1001  format(a1)
1002  format(i1)
1003  format(5(a229,/))
2000  format(///
     1 '    F I N I T E   E L E M E N T   A N A L Y S I S',
     2        '   P R O G R A M'/
     3        12x,'(C) R.L. Taylor, University of California 2011'
     4              //23x,'VERSION: ',a16/32x,a16/32x,a16/
     5                23x,'STORAGE: ',i12)
2001  format(//1x,'Specify filenames:'//
     2        1x,'Input data   (',1x,a80,') :',$)
2002  format( 1x,'Output data  (',1x,a80,') :',$)
2003  format( 1x,'Restart read (',1x,a80,') :',$)
2004  format( 1x,'Restart save (',1x,a80,') :',$)
2005  format( 1x,'Plot data    (',1x,a80,') :',$)
2006  format( 1x,'Plot device  (default: ',5x,a13,/,
     1           '(1=IBM-PHIGS 2=HP-GKS 3=WIN-INTEL 4=WIN-SALFORD):',$)
2007  format(/1x,'Files are set as:  Filename',43x,'Status'/
     1        1x,'Input   (read ) : ',1x,a80,1x,a6/
     2        1x,'Output  (write) : ',1x,a80,1x,a6/
     3        1x,'Restart (read ) : ',1x,a80,1x,a6)
2008  format( 1x,'Restart (write) : ',1x,a80,1x,a6/
     1        1x,'Plots   (write) : ',1x,a80,1x,a6/ 
     2        1x,'Plot device type: ',5x,a13,/,
     3       1x,'(1=IBM-PHIGS 2=HP-GKS 3=WIN-INTEL 4=WIN-SALFORD): ',i1/
     4        1x,'Caution, existing write files will be overwritten.'/)
2009  format( 1x,'Are filenames correct? ( y or n; s = stop) : ',$)
2010  format(/1x,'R U N N I N G    F E A P    P R O B L E M    N O W')
3001  format(' *** ERROR on read *** reinput')
      end
c
      subroutine doargs(inp,outp,res,sav,plt,idev,iswb,nargs)
c----------------------------------------------------------------------
c....  input files for batch execution  
c                                                                      
c     Batch use:  feap -iIfile -oOfile -rRfile -sSfile -pPfile -dIdev -biswb      
c     or default: feap -iIfile (defaults others)                            
c                                                                     
c     W. Wagner BS UKA 04/09
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      character*229 inp, outp, res, sav, plt,argv,cidev,ciswb
      integer index, nargs, nchars, ipos, idev, iswb, cvtc2i
c
c...  set files to blank
c
      ciswb  = ' ' 
      cidev  = ' '
      inp    = ' '
      outp   = ' '
      plt    = ' '
      res    = ' '
      sav    = ' '
      nargs  = 0 
      argv   = ' ' 
c
c.... count number of arguments on command line ! device dependent
      call doargs1(nargs)

      if (nargs .eq. 0) return
c
      do i = 1,nargs
c....   read one arguments on command line
        call doargs2(argv,nchars,i)  ! device dependent
        nchars = ipos(argv,229)
        if (index(argv,'-') .eq. 1) then
          if (index(argv,'d').eq.2) then
c....       device specification idev=1-4
            if(nchars.ge.3) then 
               cidev = argv(3:nchars)
c              idev  = cvtc2i(cidev)
            end if               
          elseif (index(argv,'b').eq.2) then
c....       batch specification iswb
            if(nchars.ge.3) then 
               ciswb = argv(3:nchars)
               iswb  = cvtc2i(ciswb)
            end if               
          else if (index(argv,'i').eq.2) then
c....       input file specification Ifile
            if(nchars.ge.3) inp = argv(3:nchars)
          else if (index(argv,'o').eq.2) then
c....       output file specification Ofile
            if(nchars.ge.3) outp = argv(3:nchars)
          else if (index(argv,'p').eq.2) then
c....       plot file specification  Pfile
            if(nchars.ge.3) plt = argv(3:nchars)
          else if (index(argv,'r').eq.2) then
c....       restart read file specification  Rfile
            if(nchars.ge.3) res = argv(3:nchars)
          else if (index(argv,'s').eq.2) then
c....       restart save file specification   Sfile
            if(nchars.ge.3) sav = argv(3:nchars)
          else
c....       error on command line
            write( *, 2000) argv(2:nchars)
            stop
          end if
        else
          write( *, 2001)  argv(1:nchars)
          stop
        end if
        end do
c
c.... check that files are correct if nargs > 0
c
      if(inp.ne.' ') then
c...    set further filenames by default
        call dosets(inp,outp,res,sav,plt)
      else
        stop
      end if
      return
2000  format( 'unknown option -> ',a)
2001  format( 'unknown argument -> ',a)
      end
c
c-----------------------------------------------------------------------
c
      subroutine test_gen(isw)
c-----------------------------------------------------------------------
c
c.... Purpose: check input file, if NEGE or GMESH is used
c     Input:
c       isw    - 1 at beginning 
c              - 2 for internal restart with RINP
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      USE comfil
      USE feapprog
      USE inptc
      USE iofile
      USE pdata2
      implicit double precision(a-h,o-z)
      logical pcomp,lnege
      character*229   newf,finp1,finp2
      character*229   fcis
      save finp2

      if(idev.lt.3) return 

c     check if nege/gmesh has been used in case of rinp
      if(isw.eq.1) then
        finp2= ' '
      else
        if1=ipos1(finp,229)  
        if(finp(if1+1:if1+9).eq.'ifeap.tmp') finp = finp2
      end if
c.... look for generation option NEGEF/GMESH/CYLT
c.... ifile is input file for NEGEF/GMESH, NEGEF/GMESH generates ifeap.tmp+nege.log/gmesh.log
      iorww=iabs(ior)
      open(unit=iorww,file=finp,status='unknown')
      read (iorww,'(a80)') newf
      rewind(iorww)
      close (iorww) 

      if(pcomp(newf(1:4),'NEGE',4)) then
c....   delete ifeap.tmp if necessary
        if1=ipos1(finp,229)  
        finp1 = finp
        finp1(if1+1:229) = ' '
        finp1(if1+1:if1+9) = 'ifeap.tmp'  ! path+ifeap.tmp
        inquire(file=finp1,exist=lnege)
        if(lnege) then
          open (iorww,file=finp1,status='unknown')
          close(iorww,status='delete')
        end if

c.....  generate now with NEGEF 
        ip=ipos(fnege,229)
        if=ipos( finp,229)  
        fcis = ' '
        fcis(1:ip) = fnege(1:ip) 
        ip = ip+1
        fcis(ip+1:ip+if) = finp(1:if)
        call start_process(fcis)

c...    new name for calculation 
        finp2 = finp ! save for rinp 
        finp  = finp1  
      end if      

      if(pcomp(newf(1:4),'GMES',4)) then
c....   delete ifeap.tmp if necessary
        if1=ipos1(finp,229)  
        finp1 = finp
        finp1(if1+1:229) = ' '
        finp1(if1+1:if1+9) = 'ifeap.tmp'  ! path+ifeap.tmp
        inquire(file=finp1,exist=lnege)
        if(lnege) then
          open (iorww,file=finp1,status='unknown')
          close(iorww,status='delete')
        endif

c.....  generate now with GMESH 
        ip=ipos(fgmesh,229)
        if=ipos(  finp,229)  
        fcis = ' '
        fcis(1:ip) = fgmesh(1:ip) 
        ip = ip+1
        fcis(ip+1:ip+if) = finp(1:if)
        call start_process(fcis)

c...    new name for calculation 
        finp2 = finp ! save for rinp 
        finp  = finp1  
      endif

      if(pcomp(newf(1:4),'CYLT',4)) then
c....   delete ifeap.tmp if necessary
        if1=ipos1(finp,229)  
        finp1 = finp
        finp1(if1+1:229) = ' '
        finp1(if1+1:if1+9) = 'ifeap.tmp'  ! path+ifeap.tmp
        inquire(file=finp1,exist=lnege)
        if(lnege) then
          open (iorww,file=finp1,status='unknown')
          close(iorww,status='delete')
        endif

c.....  generate now with CYLT -> ini must be in same directory 
        ip=ipos(fcylt,229)
        if=ipos(  finp,229)  
        fcis = ' '
        fcis(1:ip) = fcylt(1:ip) 
        ip = ip+1
        fcis(ip+1:ip+if) = finp(1:if)
        call start_process(fcis) 
        
c
c...    new name for calculation 
        finp2 = finp ! save for rinp 
        finp  = finp1  
      endif 

c.... open the files for input and output
      call test_gen_op(ior,iow,finp,fout)

      return
      end
c
      
