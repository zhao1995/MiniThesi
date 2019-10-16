c$Id:$
      subroutine pnewprob(isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    08/04/2009
c       1. Move some statements to 'pendprob.f' file.       12/05/2009
c       2. Change output file to '8' for sequence of probs. 13/05/2009
c          Add 'comsav.h' for iow_sav and fout_sav.
c       3. Add 'tsplfl' etc. for T-spline options           24/08/2010
c       4. Add call to 'fdate' for isw = 2 option           01/12/2010
c       5. Add 'perflg' for periodic non-homogeneous case   17/03/2011
c       6. Add option to specify control record items       24/01/2011
c          Check that problem has nodes, elements, etc.
c       7. Increase dimension on 'U' to have 4 entries      17/03/2011
c       8. Search for either char 47 or 92 for 'fnamr'      18/12/2011
c          Use findex to set M in fmtl.
c       9. Set 'histpltfl' false and 'hplmax' to zero       05/01/2012
c      10. Remove unused include records, add ddata.h and   19/05/2012
c          set imexfl false
c      11. Add 'bnurfl' to 'qudhsp.h'                       01/08/2012
c      12. Check that iow_sav is open before write          05/05/2012
c      13. Compute 'findex' for material file               06/11/2012
c      14. Add 'hsplfl' etc. for H-spline options           08/11/2012
c      15. Set minimum 'npd' to 300                         22/12/2012
c      16. Set finflg to .false. for hill-mandel projection 07/05/2013
c      17. Set maximum part number to 1                     19/06/2013
c      18. Increase element storage to 11 (for NURBS)       01/07/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Start a new problem

c      Note:    Statements in this routine were removed from pcontr.f
c               to permit better control on starting new problems.

c      Inputs:
c        isw    -  Switch control on actions

c      Outputs:
c        Problem control parameters through common blocks.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'allotd.h'
      include   'allotn.h'
      include   'augdat.h'
      include   'bdata.h'
      include   'cblend.h'
      include   'cdata.h'
      include   'cdat1.h'
      include   'comfil.h'
      include   'compac.h'
      include   'complx.h'
      include   'comsav.h'
      include   'contrl.h'
      include   'cornum.h'
      include   'corset.h'
      include   'c_tanfl.h'
      include   'ddata.h'
      include   'dstars.h'
      include   'edgdat.h'
      include   'eldatp.h'
      include   'elpers.h'
      include   'errchk.h'
      include   'gltran.h'
      include   'idata1.h'
      include   'idptr.h'
      include   'iodata.h'
      include   'iofile.h'
      include   'ioincl.h'
      include   'iosave.h'
      include   'mdata.h'
      include   'mxsiz.h'
      include   'nblend.h'
      include   'part0.h'
      include   'part1.h'
      include   'part3.h'
      include   'part7.h'
      include   'pdata5.h'
      include   'pdata6.h'
      include   'pfeapb.h'
      include   'pglob1.h'
      include   'pointer.h'
      include   'pload1.h'
      include   'print.h'
      include   'pscal.h'
      include   'qudshp.h'
      include   'refng.h'
      include   'sdata.h'
      include   'setups.h'
      include   'umac1.h'
      include   'vdata.h'
      include   'comblk.h'

      character  fileck*128
      logical    errs,oprt,setvar,palloc,tinput,vinput,pcomp,naninfck
      logical    contrfl,lopen
      character  cdate*24, ctext*15
      integer    isw,iii, i,j, l1,l2,l3,l4,l5,l6
      real*8     td(12)

      save

c     Jump to outputs

      if(isw.eq.2) then
        call fdate( cdate )
        go to 200
      endif

c     Close any open multiple problem file before starting new problem

      if(prob_on) then
        call pendprob

c       Remove any existing files from last problem

        call pdelfl()

c       Delete memory use

        do i = ndict,1,-1
          setvar = palloc(dlist(i),dict(i),0,iprec(i))
        end do ! i

      endif

c     Start Problem: Read and print control information

      bflg   = .false.
      gapfl  = .false.
      incf   = .false.
      intr   = .false.
      intx   = .false.
      nurbfl = .false.
      bsplfl = .false.
      bezifl = .false.
      bnurfl = .false.
      tsplfl = .false.
      hsplfl = .false.
      call fdate( cdate )
      ctext   = 'start'
      contrfl = .true.
      do while(.not.pcomp(ctext,'    ',4))
        errck  = tinput(ctext,1,td(2),8)
        if(    pcomp(ctext,'node',4) .or. pcomp(ctext,'numnp',5)) then
          numnp = nint(td(2))
          contrfl = .false.
        elseif(pcomp(ctext,'elem',4) .or. pcomp(ctext,'numel',5)) then
          numel = nint(td(2))
          contrfl = .false.
        elseif(pcomp(ctext,'mate',4) .or. pcomp(ctext,'nummat',5)) then
          nummat = nint(td(2))
          contrfl = .false.
        elseif(pcomp(ctext,'dime',4) .or. pcomp(ctext,'ndm',3)) then
          ndm = nint(td(2))
          contrfl = .false.
        elseif(pcomp(ctext,'dofs',4) .or. pcomp(ctext,'ndf',3)) then
          ndf = nint(td(2))
          contrfl = .false.
        elseif(pcomp(ctext,'elno',4) .or. pcomp(ctext,'nen',3)) then
          nen = nint(td(2))
          contrfl = .false.
        elseif(pcomp(ctext,'add',3)  .or. pcomp(ctext,'nad',3)) then
          nad = nint(td(2))
          contrfl = .false.
        elseif(pcomp(ctext,'prop',4)  .or. pcomp(ctext,'npd',3)) then
          npd = nint(td(2))
          contrfl = .false.
        elseif(pcomp(ctext,'upro',4)  .or. pcomp(ctext,'nud',3)) then
          nud = nint(td(2))
          contrfl = .false.
        elseif(contrfl) then
          errck = vinput(ctext,15,td(1),1)
          if(nint(td(1)).ge.0) then
            numnp  = nint(td(1))
            numel  = nint(td(2))
            nummat = nint(td(3))
            ndm    = nint(td(4))
            ndf    = nint(td(5))
            nen    = nint(td(6))
            nad    = nint(td(7))
            npd    = nint(td(8))
            nud    = nint(td(9))
            go to 101
          endif
        endif
      end do ! while
101   nnn    = 0

c     Adjust storage for material parameters

      npd    = max(npd,300)
      nud    = max(nud,150)
      ndd    = npd + nud + 1

c     Star node/element initialization

      starnd = 0
      starel = 0

c     Blending function initialization

      numsn  = 0
      numsd  = 0
      numbd  = 0

c     NURB blending function initialization

      nursn  = 0
      nursd  = 0
      nurbd  = 0

c     Load table initialization

      ldtot  = 0
      ldnum  = 0
      ldprp  = 0
      sptot  = 0
      spnum  = 0
      ldflg  = .false.
      spflg  = .false.
      perflg = .false.

c     Contact array initialization

      numcels = 0
      ncen    = 0
      optflg  = .false.
      optmsh  = .false.
      opthoit = .false.

c     Serial & parallel solution by unblocked equations

      pfeap_blk  = .false.
      pfeap_glob = .false.

c     Set filenames for multiple problem case

      if(irdef.ne.ior) then

        inquire(unit=ior,name=fnamp,exist=errs)

        prob_on = .false.
        if(errs) then

c         Set multiple problem flag

          prob_on = .true.

c         Save master output file name and unit number

          i = index(flog,' ')
          if(nprob.eq.0) then
            if(isw.gt.0) write(iow,2017) flog(1:i-1)
            iow_sav  = iow
            fout_sav = fout
          endif

c         Extract file name

          i = index(fnamp,' ')
          if(i.eq.0) i = 128
          do j = i,1,-1
            if(pcomp(fnamp(j:j),char(47),1) .or.       ! char(47) = '/'
     &         pcomp(fnamp(j:j),char(92),1)) go to 110 ! char(92) = '\'
          end do ! j
          j = 0
110       fnamr = fnamp(j+1:j+21)

c         Set new plot file name

          fnamr(1:1) = 'P'
          fplt(1:128) = ' '
          fplt(1: 17) = fnamr
          i = index(fplt,'.')
          if(i.gt.0) then
            fplt(i: 21) = ' '
          endif
          i = min(index(fplt,' '), 16)
          if(i.eq.0) then
            i = 16
          endif

c         Increment problem counter or delete output file

          if(keepfl) then
            nprob = nprob + 1
          else
            close(unit = iow, status = 'delete')
            keepfl = .true.
            nprob  = max(1,nprob)
          endif

c         Add problem counter to name

          write(fplt(i:i+2),'(a)') '000'
          if(nprob.lt.10) then
            write(fplt(i+2:i+2),'(i1)') nprob
          elseif(nprob.lt.100) then
            write(fplt(i+1:i+2),'(i2)') nprob
          elseif(nprob.lt.1001) then
            write(fplt(  i:i+2),'(i3)') nprob
          else
            write(*,*) 'Exceeded limit of multiple files (PCONTR)'
          endif

c         Set file names for new problem

          if(isw.gt.0) then
            iow  = 8
            fout = fplt
            fout(1:1) = 'O'
            fres = fplt
            fres(1:1) = 'R'
            fsav = fplt
            fsav(1:1) = 'S'

c           Create clean output file

            inquire(file=fout,exist=initf)
            if(initf) then
              open (unit=iow,file=fout,status='old')
              close(unit=iow,          status='delete')
            endif
            open(unit=iow,file=fout,status='new')
            if(nprob.gt.1) write(ilg,2019)
            write(ilg,2020) nprob,fout
            inquire(unit=iow_sav,opened=lopen)
            if(lopen) write(iow_sav,2021) nprob
          endif

c       Error in file structure

        else
          write(  *,3003)
          write(ilg,3003)
          call plstop()
        endif

c     Single problem solution

      else
        prob_on = .false.
      endif

c     Zero pointer array

      setvar = palloc( 0 ,'START', 0 , 0 )

c     Set initial values for nan and inf checks

      setvar = naninfck(td(1),1,0)

c     Set contact flags for a new problem

      call contact (300)

c     Set rigid body flags for a new problem

      call rigidb (0,0,errs)

c     Zero number of dictionary entries

      ndict     = 0

c     If number of nodes, or elements is zero compute number from data

      if(nocount) then
        oprt        =  prt
        prt         = .false.
        ucount      = .true.
        tsplfl      = .false.
        hsplfl      = .false.
        call pnums()
        irecrd(isf) =  2
        prt         =  oprt
        ucount      = .false.

c       Star node/element re-initialization

        starnd = 0
        starel = 0
      else
        tsplfl      = .false.
        hsplfl      = .false.
      endif

c     Output problem size data

200   write(iow,2000) head,cdate,versn,fincld(isf),
     &               numnp,numel, ndm,ndf,nad,nen, nummat,npd,nud

c     Check that problem has nodes elements, etc.

      if(min(numnp,numel,nummat, ndm,ndf,nen).eq.0) then
        call plstop()
      endif

c     Initialize clock

      call stime()

c     Set parameters for page eject and rotation dof

      o   = '    '
      errck = .false.
      lsave = .false.
      lkflg = .false.
      leflg = .false.
      lcflg = .false.
      initf = .true.
      cxifl = .false.
      eanfl = .false.
      ebcfl = .false.
      ebsfl = .false.
      curfl = .false.
      edifl = .false.
      efcfl = .false.
      eprfl = .false.
      espfl = .false.
      finflg= .false.
      surfl = .false.
      basfl = .false.
      boufl = .false.
      cprfl = .false.
      disfl = .false.
      forfl = .false.
      lfrfl = .false.
      angfl = .false.
      reafl = .false.
      intfl = .false.
      tiefl = .true.
      tief  = .false.
      damfl = .false.
      masfl = .false.
      stifl = .false.

c     Dynamic contact flags

      rattlfl = .false.
      shakefl = .false.
      imexfl  = .false.

c     History plot flag

      histpltfl = .false.
      hplmax    = 0

c     Rotation parameters

      do i = 1,50
        ia(1,i)  = 1
        ia(2,i)  = 2
        ir(1,i)  = 0
        ir(2,i)  = 0
        ea(1,i)  = 1
        ea(2,i)  = 2
        er(1,i)  = 0
        er(2,i)  = 0
        ia3(1,i) = 1
        ia3(2,i) = 2
        ia3(3,i) = 3
        ir3(1,i) = 0
        ir3(2,i) = 0
        ir3(3,i) = 0
        ea3(1,i) = 1
        ea3(2,i) = 2
        ea3(3,i) = 3
        er3(1,i) = 0
        er3(2,i) = 0
        er3(3,i) = 0
        inord(i) = 0
        exord(i) = 0
        do j = 1,30
          ipord(j,i) = 0
          epord(j,i) = 0
        end do ! j
      end do ! i
      nadd   = 0
      nprof  = 0
      nsurf  = 0
      nbouf  = 0
      ndisf  = 0
      nforf  = 0
      nforl  = 0
      nangf  = 0
      nintf  = 0
      neang  = 0
      nebcs  = 0
      nedis  = 0
      nefrc  = 0
      nepro  = 0
      ndamf  = 0
      ncurv  = 0
      nespi  = 0
      nmasf  = 0
      nstif  = 0
      nbasf  = 0
      neule  = 0

c     Zero global parameters

      gtypfl = .false.
      gdeffl = .false.
      gtdofl = .false.
      gomgfl = .false.
      grayfl = .false.
      groufl = .false.
      groupl = .false.
      if(    ndm.le.2) then
        g2type = 2           ! default plane strain
      elseif(ndm.eq.3) then
        g2type = 7           ! default 3-d
      else
        g2type = 9           ! unspecified
      endif
      gdtype = 1
      gtdof  = 0
      gref   = 0
      geqnum = 0
      gpart  = 0
      do i = 1,3
        units(i)  = 1.0d0
        grefx(i)  = 0.0d0
        gtref(i)  = 0.0d0
        gomega(i) = 0.0d0
        gomex(i)  = 0.0d0
        gomev(i)  = 0.0d0
      end do ! i
      do i = 1,2
        gray(i) = 0.0d0
      end do ! i
      do i = 1,14
        gfac(i) = 0.0d0
      end do ! i
      gquadn =  0.0d0        ! Default quadrature: Gauss
      gaugm  =  1.0d0        ! Default augmenting on
      augf   =  1.0d0        ! Augmenting factor multiplier
      augg   =  0.0d0        ! Maximum augmenting gap
      auggfl = .false.
      augmfl = .false.
      lvaug  =  0
      lvsol  =  0

c     Set pointers for allocation of mesh arrays

      nen1      = nen + 9
      nie       = 13
      nst       = max(nen*ndf + nad,1)
      nneq      = ndf*numnp
      do i = 1,ndf
        ndfl(i) = 0
        ndfo(i) = 0
        ndog(i) = 10
      end do ! i

c     Set pointers for allocation of mesh arrays

      if(cplxfl) then
        ipc = 2
      else
        ipc = 1
      end if

c     Allocate size for arrays for mesh and solution vecors

      l1   = ndm*numnp
      l2   = max(ndf*numnp,1)
      l3   = max(nen+1,7*nst,21)
      l4   = numnp*max(ndf,ndm)*ipc
      l5   = ndf*nen
      l6   = max(1,numel)

c     Allocate and zero arrays

      setvar = palloc( 26,'DR   ',l4          ,  2)
      setvar = palloc( 34,'LD   ',l3          ,  1)
      setvar = palloc( 35,'P    ',nst*3       ,  2)
      setvar = palloc( 36,'S    ',nst*nst*2   ,  2)
      setvar = palloc( 39,'TL   ',nen         ,  2)
      setvar = palloc( 41,'UL   ',nst*14      ,  2)
      setvar = palloc( 44,'XL   ',max(4,nen)*3,  2)
      setvar = palloc( 25,'D    ',nummat*ndd  ,  2)
      setvar = palloc( 32,'IE   ',nummat*nie  ,  1)
      setvar = palloc(240,'IEDOF',nummat*l5   ,  1)
      setvar = palloc( 31,'ID   ',l2*2        ,  1)
      setvar = palloc( 33,'IX   ',nen1*l6     ,  1)
      setvar = palloc(190,'NDTYP',numnp       ,  1)
      setvar = palloc(100,'RIXT ',numnp       ,  1)
      setvar = palloc(181,'RBEN ',l6          ,  1)
      setvar = palloc( 43,'X    ',l1          ,  2)
      setvar = palloc( 45,'ANG  ',numnp       ,  2)
      setvar = palloc( 46,'ANGL ',nen         ,  2)
      setvar = palloc( 27,'F    ',2*l2        ,  2)
      setvar = palloc( 28,'F0   ',4*l2        ,  2)
      setvar = palloc( 29,'FPRO ',2*l2        ,  1)
      setvar = palloc( 30,'FTN  ',4*l2        ,  2)
      setvar = palloc( 38,'T    ',numnp       ,  2)
      setvar = palloc( 40,'U    ',4*l2*ipc    ,  2)
      setvar = palloc( 89,'NREN ',numnp*2     ,  1)
      if(ldtot.gt.0) then
        setvar = palloc(265,'LDTAB',ldtot*12  ,  1)
      endif

c     Set ID address pointers

      id31    = np(31)
      idpt(1) = np(31)

c     Set pointers

      npid    = np(31)         ! ID
      npix    = np(33)         ! IX
      npuu    = np(40)         ! U
      npxx    = np(43)         ! X
      nprn    = np(89)         ! NREN
      npty    = np(190)        ! NDTYP

c     Set initial numbering in renumber vector and mark nodes as unused.

      do i = 0,numnp-1
        mr(np( 89)+i      ) = i+1  ! Remap list
        mr(np( 89)+i+numnp) = i+1  ! Reverse list
        mr(np(190)+i      ) = 0
      end do ! i

c     Set pointers for allocation of loads/partition data

      do j = 1,5
        tflp(j)   = .true.
        flp(9,j)  = .false.
        scale(j)  = .false.
        nittyp(j) = -3
        nsolver(j)= solver
        nitolp(j) =  1.d-08
        natolp(j) =  1.d-16
      end do ! j
      do j = 1,ndf
        ndfst(j,1) = 1
        ndfp(j)    = 1
        ndfg(j)    = ndfp(j)
      end do ! j
      nqp(1)    =  1
      npart     =  1
      maxpart   =  1
      nopart    = .true.

c     Open file to store material data

      inquire(unit=iwd,name=fileck, opened=errs)
      if(errs) then
        write(*,*) 'UNIT=',iwd,' FILE=',fileck
      endif

c     Locate offset for setting material file

      fmtl = finp
      do j = len_trim(fmtl),1,-1
        if(pcomp(fmtl(j:j),char(47),1) .or.       ! char(47) = '/'
     &     pcomp(fmtl(j:j),char(92),1)) go to 100 ! char(92) = '\'
      end do ! j
      j = 0
100   findex = j + 1
      fmtl(findex:findex) = 'M'
      open(unit=iwd,file=fmtl,status='unknown')
      rewind iwd

c     Input a mesh from binary file (if it exists)

      if(bflg) then
        call bmesh(mr(np(32)),hr(np(25)),mr(np(31)+nneq),hr(np(43)),
     &             mr(np(33)),hr(np(27)),hr(np(38)),hr(np(45)),
     &             ndd,nie,ndm,ndf,nen,nen1,numnp,numel,nummat)
        close(ios)
        prt   = .true.
        iii   = -2
      else
        iii   =  0
      endif

c     Input mesh data from file

      call pmesh(iii,prt,prth)

c     Perform simple check on mesh to ensure basic data was input

      call meshck(mr(np(32)),mr(np(190)),mr(np(33)),nie,nen,nen1,
     &            numnp,numel,nummat,errs)

      if(errs) then
        call plstop()
      endif

c     Compute boundary nodes (before ties)

      if(tiefl) then
        setvar = palloc( 78,'EXTND',numnp ,1)
        call pextnd()
        tiefl  = .false.
      endif

      tfl = .true.

c     Input/output formats

2000  format(1x,19a4,a3//5x,'Solution date: ',a//14x,a/14x,a/
     &                /5x,'Input Data Filename: ',a/
     &                /5x,'Number of Nodal Points  - - - - - - :',i9
     &                /5x,'Number of Elements  - - - - - - - - :',i9/
     &                /5x,'Spatial Dimension of Mesh - - - - - :',i9
     &                /5x,'Degrees-of-Freedom/Node (Maximum) - :',i9
     &                /5x,'Equations/Element       (Maximum) - :',i9
     &                /5x,'Number Element Nodes    (Maximum) - :',i9/
     &                /5x,'Number of Material Sets - - - - - - :',i9
     &                /5x,'Number Parameters/Set   (Program) - :',i9
     &                /5x,'Number Parameters/Set   (Users  ) - :',i9)

2017  format(/'  Problem definitions are specified by include files.'
     &      //'  Output for each problem is written to separate files.'
     &      //'  Check file ',a,' for problem list and errors.')

2019  format(/'  ',70('-'))

2020  format(/'  --> Problem',i4,': Output in file: ',a)

2021  format(/'  --> End Problem',i4)

3003  format(/' *ERROR* PCONTR: File name error')

      end
