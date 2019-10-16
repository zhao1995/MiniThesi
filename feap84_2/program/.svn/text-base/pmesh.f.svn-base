c$Id:$
      subroutine pmesh(iii,prt,prth)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Change 244 to 256 in palloc of 'NFORC' and
c          change np(244) to np(256) in call genvec         02/11/2006
c       2. Change call pprint('  Enter ...') to write(*,*) 'Enter ...
c          Set xxx(1:13) correctly and remove format 2000.  13/11/2006
c       3. Add 'sweep' option for 2-d to 3-d meshes         17/03/2007
c       4. Set length of lp_file to 15, length for cc/c2    17/04/2007
c       5. Add 'elabel' to 'pelmin' call                    09/07/2007
c       6. Remove reallocation of BNODE for supernodes      15/01/2008
c       7. Remove 'esse' b.c. option and add 'curve' option 13/11/2008
c          Add option for data reads on command line. Use
c          for 'meshrad', 'meshang', 'meshth'
c       8. Add 'load' table input option                    10/01/2009
c       9. Add 'set'/'add' option for 'cbou' command        19/02/2009
c      10. Allow 'rese't of angles, forces, displacements   22/02/2009
c          and boundary codes to zero.
c      11. Add 'spin' command to rotate parts               07/04/2009
c      12. Add error check that 'fpro' does not define      09/11/2009
c          two different proportional load numbers
c      13. Use setext to assign file extenders & set number 20/12/2010
c      14. Increase dimension on spin table to allow edges  04/01/2011
c      15. Add 'peri'odic case                              17/03/2011
c      16. Add options to periodic case; change xc to xd    19/04/2011
c      17. Add 'tx(*)' to umesh and umshlib calls           26/09/2011
c      18. Add 'triad' to input options                     11/02/2013
c      19. Remove periodic prtype = 4                       13/04/2013
c      20. Remove 'ndm' from call to periodic               03/05/2013
c      21. Add 'tbou' to fix all dof's (Taylor b.c.)        15/05/2013
c      22. Set 'prth' to same as 'prt' for prin/nopr        24/09/2013
c      23. Set 'iclink' to 0 to start mesh inputs           28/10/2013
c      24. Extract 'reafi' from 'xxx' (up to 128 chars)     29/11/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Data input routine for mesh description

c      Inputs:
c         iii        - Initialization indicator
c         prt        - Flag, print input data if true
c         prth       - Flag, print title/header if true

c      Outputs:
c         Depends on commands specified
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'bdata.h'
      include  'cblend.h'
      include  'cblktr.h'
      include  'cdata.h'
      include  'cdat1.h'
      include  'cdat2.h'
      include  'chdata.h'
      include  'codat.h'
      include  'comfil.h'
      include  'complx.h'
      include  'corset.h'
      include  'corfil.h'
      include  'cornum.h'
      include  'crotas.h'
      include  'debugs.h'
      include  'dstars.h'
      include  'edgdat.h'
      include  'eldata.h'
      include  'elpers.h'
      include  'eqslv.h'
      include  'hlpdat.h'
      include  'iodata.h'
      include  'iofile.h'
      include  'iosave.h'
      include  'linka.h'
      include  'mdata.h'
      include  'meshtx.h'
      include  'modreg.h'
      include  'pdata3.h'
      include  'pload1.h'
      include  'pointer.h'
      include  'p_point.h'
      include  'p_ptname.h'
      include  'prflag.h'
      include  'prld1.h'
      include  'region.h'
      include  'rigid1.h'
      include  'rigid2.h'
      include  'sdata.h'
      include  'setups.h'
      include  'sldata.h'
      include  'trdata.h'
      include  'umac1.h'
      include  'comblk.h'

      include  'p_int.h'

      integer   list   , numesh
      parameter(list=97, numesh=10)

      logical   setvar,palloc, lp_in, lopen
      logical   prt,error,pcomp,lmesh,errck,pinput,tinput,vinput
      logical   prth, umesh, readfl, savefl, snofl,flgco,eflg
      logical   ldonf,ldfor,lddis
      character lp_file*15, elabel*38
      character wd(list)*4, cc*4,c2*8,fext*8,swd(2)*11
      integer   ed(list)  , i,j,ii,jj, iii, isd, ibn, lmr, side, face
      integer   ll,llo,lp_lun
      integer   nsblk,nblend, nstar,estar, starsvn,starsve
      integer   numn2,nume2,numm2,ndm2,ndf2,nen2,nnmax
      integer   iorold, nn,n1,n2, nseg,ndir,nstyp, is(10)
      real*8    td(6),  xd(3), xx(3,100), xp(3,100), shp(100)

      save

c     List of command names

      data swd /'Cartesian','Cylindrical'/

      data elabel /'E l e m e n t   C o n n e c t i o n s'/

      data wd/'coor','elem','mate','boun','forc','temp','end ','prin',
     1        'nopr','titl','bloc','pola','ebou','angl','sloa','cons',
     2        'sphe','btem','icon','pars','nopa','trib','para','efor',
     3        'eang','cbou','cfor','cang','foll','slav','rese','sblo',
     4        'curv','rota','setn','setr','btra','fpro','cpro','regi',
     5        'tran','damp','mass','stif','csur','ereg','reac','manu',
     6        'body','glob','shif','disp','edis','cdis','debu','side',
     7        'face','snod','blen','move','rigi','moda','flex','base',
     8        'epro','mpro','loop','next','file','cdam','cmas','csti',
     9        'ebas','cbas','eule','ceul','rfor','lfor','load','swee',
     a        'spin','peri','expl','tria','tbou',
     s        '*nod','*ele',
     u        'mes1','mes2','mes3','mes4','mes5','mes6','mes7','mes8',
     u        'mes9','mes0'/

      data ed/    0,     0,     0,     0,     0,     1,     0,     0,
     1            0,     2,     0,     0,     0,     0,     6,     6,
     2            0,     1,     6,     3,     3,     2,     0,     0,
     3            0,     0,     0,     0,     1,     1,     1,     2,
     4            1,     2,     2,     2,     3,     1,     1,     1,
     5            1,     1,     1,     1,     0,     1,     2,     4,
     6            1,     0,     0,     0,     0,     0,     0,     0,
     7            6,     0,     0,     2,     1,     1,     1,     0,
     8            0,     1,     3,     3,     3,     0,     0,     0,
     9            2,     2,     1,     1,     1,     1,     0,     1,
     a            1,     1,     1,     1,     1,
     s            2,     2,
     u            5,     5,     5,     5,     5,     5,     5,     5,
     u            5,     5/

c     Write message to log file

      if(iii.ge.0) write(ilg,4000) finp

c     Initialize arrays and set error detection values

      lp_in = .true.
      error = .false.
      lddis = .false.
      ldfor = .false.
      ldonf = .false.
      lmesh = .false.
      ldflg = .false.
      spflg = .false.
      ll    = 1
      nneq  = ndf*numnp
      nstar = list - numesh - 1
      estar = list - numesh
      if(iii.ge.0) then
        frotas      = .false.
        nffl        = .false.
        nmfl        = .false.
        snofl       = .false.
        anglefl     = .false.
        eulerfl     = .false.
        triadfl     = .false.
        oldfl       = .false.
        iclink      = 0
        ldnum       = 0
        spnum       = 0
        nreg        = 0
        mxreg       = 0
        nblend      = 0
        npstr       = 0
        nrigid      = 0
        numg        = 0
        numsl       = 0
        nim         = 0
        nre         = 0
        nio         = 0
        neo         = 0
        mao         = 0
        ibn         = 20
        isd         = 16
        side        = 0
        face        = 0
        prtype      = 0
        rotopt      = 1

        last_elm    = 0
        ipart       = 0

c       Set angles, boundary code/forced values to zero

        do j = 1,3
          x0(j)    = 0.0d0
          xr(j)    = 0.0d0
          do i = 1,3
            tr(i,j) = 0.0d0
          end do ! i
          tr(j,j) = 1.d0
        end do ! j
        trdet = 1.d0
        if(iii.eq.0) then
          prt = .true.

c         Set node type to undefined

          do j = 0,numnp-1
            mr(np(190)+j) = -1
          end do ! j

          do j = 1,50
            prldv(j) = 1.d0
          end do ! j

c         Insert user functions into list

          do j = 1,10
            i = list-numesh+j
            if(.not.pcomp(umshc(j),wd(i),4)) then
              wd(i) = umshc(j)
              ed(i) = 0
            endif
          end do ! j

        endif

      elseif(iii.eq.-2) then
        nio   = 0
        neo   = 0
        mao   = 0
        nreg  = 0
      endif
100   if(ior.lt.0) then
        write(*,*) ' Enter "help" for list of commands, "end" to exit'
        xxx       = ' '
        xxx(1:13) = '    Mesh    >'
        write(xxx(10:12),'(i3)') ll
        call pprint(xxx)
      endif
      errck = tinput(tx,8,td,0)
      if(errck) go to 100
      utx(1) = tx(1)
      utx(2) = tx(2)
      cc     = tx(1)(1:4)
      c2     = tx(2)(1:4)

c     Look for data on input line

      meshrad = 0.0d0
      meshang = 0.0d0
      meshth  = 0.0d0
      do i = 2,7
        if(pcomp(tx(i),'radi',4)) then
          call setval(tx(i+1),15,meshrad)
        elseif(pcomp(tx(i),'thic',4)) then
          call setval(tx(i+1),15,meshth)
        elseif(pcomp(tx(i),'angl',4)) then
          call setval(tx(i+1),15,meshang)
        endif
      end do ! i

c     Check for options to read data from file

      if( pcomp(cc,'read',4) ) then
        lmesh = readfl(tx(2))
        if(lmesh) then
          llo = ll
        else
          ll  = llo
        endif
        go to 100
      endif
      if(pcomp(cc,'save',4)) then
        lsave = savefl(tx(2))
        go to 100
      endif
      if(ior.lt.0.and.pcomp(cc,'help',4)) then
        call phelp(c2,wd,ed,list,'MESH')
        go to 100
      endif
      go to 120
110   call  errclr ('PMESH ')
      go to 100

c     Transfer to input option based on value of 'cc'

120   do i = 1,list
        if(pcomp(cc,wd(i),4)) go to 130
      end do ! i

c     User mesh commands

      if(.not. pcomp( cc, ' ', 1 ) ) errck = umesh(cc,tx,prt)
      if(.not. errck .and. ior.lt.0) call errclr('PMESH ')
      go to 100
130   ll = ll + 1
      if(iii.ge.0) write(ilg,4001) tx(1),tx(2)

c     [coor]dinates  - Nodal coordinate data input
c     [coor]<all>    - Nodal coordinate data input (no parsing)
c     [coor]<add>    - Accumulate the numbers starting at nio
c     [coor]<bina>ry - Input coordinates in binary mode

      if(i.eq.1) then
        if(pcomp(c2,'all',3)) then
          call pcrdrd(hr(np(43)),mr(np(190)),prt)
        elseif(pcomp(c2,'bina',4)) then
          call pcrbrd(hr(np(43)),mr(np(190)),prt)
        else
          starsvn = starnd
          if(pcomp(c2,'add',3)) then
            starnd  = nio
          endif
          call genvec(ndm,ndm,hr(np(43)),' Coordinates',
     &                prt,prth,error,.true.)
          if(nio.gt.starnd) then
            numnp = max(numnp,nio)
          endif
          starnd = starsvn
        endif

c     [elem]ent      - Data input
c     [elem,all]     - Input all connection data
c     [elem,add]     - Input all connection data
c     [elem<bina>ry  - Input element connections in binary mode
c     [elem,gene,xxx]- Set generation array to xxx

      elseif(i.eq.2) then
        if(pcomp(c2,'all',3)) then
          call pelmrd(mr(np(33)),mr(np(181)),prt)
        elseif(pcomp(c2,'bina',4)) then
          call pelbrd(mr(np(33)),mr(np(181)),prt)
        else
          starsve = starel
          if(pcomp(c2,'add',3)) then
            starel  = neo
          endif
          call pelmin(tx(2),mr(np(34)),mr(np(33)),mr(np(181)),nen1,
     &                prt,prth,error,elabel)
          if(neo.gt.starel) then
            numel = max(numel,neo)
          endif
          starel = starsve
        endif

c     [mate]rial,ma: Data input for material set ma

      elseif(i.eq.3) then
        setvar = palloc(151,'USER1', ndf*(nen+1),1)
        call pmatin(tx,hr(np(25)),hr(np(41)),hr(np(44)),hr(np(39)),
     &                 hr(np(36)),hr(np(35)),mr(np(34)),mr(np(32)),
     &                 mr(np(240)),mr(np(151)),prt,prth)
        setvar = palloc(151,'USER1', 0          ,1)

c     [boun]dary codes,<set,add> - input nodal restraint conditions

      elseif(i.eq.4) then
        if(pcomp(c2,'add',3)) then
          setvar = palloc(151,'USER1', ndf*numnp,1)
          call pzeroi(mr(np(151)),ndf*numnp)
          fp(1) = np(151)
          fp(2) = np(31) + nneq
        else
          fp(1) = np(31) + nneq
        endif
        call pbouin(mr(np(34)),mr(fp(1)),prt,prth)
        if(pcomp(c2,'add',3)) then
          do ii = 0,ndf*numnp-1
            mr(ii+fp(2)) = mr(ii+fp(2)) + mr(ii+fp(1))
          end do ! ii
          setvar = palloc(151,'USER1', 0,1)
        endif

c     [forc]e,<set,add>: Data input for nodal generalized forces

      elseif(i.eq.5) then
        if(pcomp(c2,'add',3)) then
          setvar = palloc(151,'USER1', ndf*numnp,2)
          call pzero(hr(np(151)),ndf*numnp)
          fp(1) = np(151)
          if(ldflg) then
            fp(2) = np(26)
            ldfor = .true.
          else
            fp(2) = np(27)
          endif
        else
          if(ldflg) then
            fp(1) = np(26)
            ldfor = .true.
          else
            fp(1) = np(27)
          endif
        endif
        call genvec(ndf,ndf,hr(fp(1)),' Forces',prt,prth,error,.false.)
        if(pcomp(c2,'add',3)) then
          do ii = 0,ndf*numnp-1
            hr(ii+fp(2)) = hr(ii+fp(2)) + hr(ii+fp(1))
          end do ! ii
          setvar = palloc(151,'USER1', 0,2)
        endif

c     [temp]erature data input

      elseif(i.eq.6) then
        call genvec(1,1,hr(np(38)),' Temperatures',prt,prth,error,
     &             .false.)

c     [end] of mesh data inputs

      elseif(i.eq.7) then
        if(lsave) then
          if(iii.ge.0) write(ilg,3004)
          write(iow,3004)
          if(ior.lt.0) then
            write(iow,3004)
          endif
          error = .true.
        endif

        close(unit=iwd)

c       Complete surface load sets

        call ploadi(prt,prth,error,2)

c       If load table input open, close and save

        if(ldflg) then
          if(prt) write(iow,2013) ldnum
          call pldend(ldfor,lddis,ldflg)
          spflg = .false.
        endif

c       Exit if any errors occurred in mesh data

        if(error) then
          call plstop()

        elseif(iii.ge.0) then

c         Perform delayed mesh generation steps

          if(numbd.gt.0 .and. iii.ge.0) then
            if(snofl) then
              call pblendm(isd,ibn,ndm,nen1,prt,prth,.true.,.true.)
            else
              if(ior.le.0) then
                write(*,3005)
              endif
              write(ilg,3005)
              write(iow,3005)
              call plstop()
            endif
          endif

c         Initialize transformation matrices for all configurations

          call chkrot(mr(np(32)),mr(np(33)),nie,nen1)
          if (frotas) then
            setvar = palloc(111,'TEMP1', 3*numnp, 2)
            call pshsurf(hr(np(43)),hr(np(82)),hr(np(111)),
     &                   mr(np(33)),mr(np(32)),mr(np(81)))
            setvar = palloc(111,'TEMP1', 0, 2)

            call setrot(hr(np(43)),mr(np(81)),hr(np(82)),hr(np(83)),
     &                  numnp,rotopt)

c           Compress rotational data storage

            call cmprot(mr(np(81)),hr(np(82)), numnp, lmr)
            setvar = palloc(82,'MR   ', lmr*54,2)

c           Initialize rotational parameters

            call updrot(hr,ndf,hr(np(82)),mr(np(81)),numnp,0)

          endif

        endif

c       Set load group flag to false and reset residual array size

        ldflg = .false.
        if(ldonf) then
          ldonf = .false.
          setvar = palloc( 26,'DR   ',numnp*max(ndf,ndm)*ipc,2)
        endif
        return

c     [prin]t/[nopr]int of input data

      elseif(i.eq.8 .or. i.eq.9) then
        prt  = i.eq.8
        prth = prt

c     [titl] - set title prints on/off

      elseif(i.eq.10) then
        if(pcomp(c2,'off',3)) then
          prth = .false.
        else
          prth = .true.
        endif

c     [bloc]k - generate block of nodes and elements

      elseif(i.eq.11) then
        if(iii.lt.0) write(iow,3002)
        if(pcomp(c2,'old',3)) then
          oldfl = .true.
        else
          oldfl = .false.
        endif
        call blkgen(ndm,nen1,hr(np(43)),mr(np(33)),mr(np(181)),prt,prth)

c     [pola]r - convert polar to cartesian coordinates

      elseif(i.eq.12) then
        call polar(mr(np(190)),hr(np(43)),ndm,prt,prth)

c     [ebou] - set edge boundary constraints

      elseif(i.eq.13) then
        ebcfl = .true.
        call setext('ebou',nebcs, fext, .true.)
        if(.not.pcomp(c2,'set',3)) then
          c2 = 'add'                 ! default mode 'add' for edge b.c.
        endif
        call plinka(fext,c2,'   ')

c     [angl]e - set boundary angles

      elseif(i.eq.14) then
        call genvec(1,1,hr(np(45)),' Angles',prt,prth,error,.false.)
        anglefl = .true.

c     [sloa]ds - set surface loadings

      elseif(i.eq.15) then
        call ploadi(prt,prth,error,1)

c     [para]meter - set parameter variables

      elseif(i.eq.16 .or. i.eq.23) then
        coflg = .true.
        call pconst(prt)

c     [sphe]re - convert spherical to cartesian coordinates

      elseif(i.eq.17) then
        call sphere(mr(np(190)),hr(np(43)),ndm,prt,prth)

c     [btem] - input block of interpolated temperatures

      elseif(i.eq.18) then
        call blktem(ndm,hr(np(38)),prt,prth)

c     [icon] - input and check contact data

      elseif(i.eq.19) then
c       call incon(hr(np(43)),ndm,prt,error)

c     [pars]ing/[nopa]rsing of statements

      elseif(i.eq.20 .or. i.eq.21) then
        coflg = i.eq.20

c     [trib] - triangular block generator

      elseif(i.eq.22) then
        if(iii.lt.0) write(iow,3002)
        call blktri(ndm,nen1,hr(np(43)),mr(np(33)),prt,prth)

c     [efor] set edge force constraints

      elseif(i.eq.24) then
        efcfl = .true.
        call setext('efor',nefrc, fext, .true.)
        if(.not.pcomp(c2,'set',3)) then
          c2 = 'add'                 ! default mode 'add' for edge force
        endif
        call plinka(fext,c2,'   ')

c     [eang] set edge angle constraints

      elseif(i.eq.25) then
        anglefl = .true.
        eanfl   = .true.
        call setext('eang',neang, fext, .true.)
        if(.not.pcomp(c2,'set',3)) then
          c2 = 'add'                 ! default mode 'add' for edge angle
        endif
        call plinka(fext,c2,'   ')

c     [cbou] set coordinate boundary constraints - based on coordinates

      elseif(i.eq.26) then
        boufl = .true.
        call setext('cbou',nbouf, fext, .true.)
        if(.not.pcomp(c2,'set',3)) then
          c2 = 'add'                 ! default mode 'add' for edge angle
        endif
        call plinka(fext,c2,'   ')

c     [cfor] set coordinate nodal forces - based on coordinates

      elseif(i.eq.27) then
        forfl = .true.
        call setext('cfor',nforf, fext, .true.)
        if(.not.pcomp(c2,'set',3)) then
          c2 = 'add'                 ! default mode 'add' for edge angle
        endif
        call plinka(fext,c2,'   ')

c     [cang] set coordinate angle - based on coordinates

      elseif(i.eq.28) then
        anglefl = .true.
        angfl   = .true.
        call setext('cang',nangf, fext, .true.)
        call plinka(fext,c2,'   ')

c     [foll]ower loads

      elseif(i.eq.29) then
        fext = 'follower'
        call plinka(fext,c2,'   ')
        nfol = iclink

        setvar = palloc(129,'FOLLI ',nfol*2,1)
        setvar = palloc(130,'FOLLR ',nfol,  2)
c       call pfollo(nfol,hr(np(130)),mr(np(129)),prt,prth)

c     [slav]ed nodes to subspace matrix

      elseif(i.eq.30) then
30      errck = pinput(td,1)
        j     = nint(td(1))
        if(j.gt.0 .and. j.le.numnp) then
          errck = palloc(205,'NSLAV',numg+1,1)
          mr(np(205)+numg) = j
          numg             = numg + 1
        else
          call iprint(mr(np(205)),1,numg,1,'Slaved Nodes')
          go to 100
        endif
        go to 30

c     [rese]t  [angl]e: Set all boundary angle values to zeron
c     [rese]t  [disp]l: Set all nodal displace values to zeron
c     [rese]t  [forc]e: Set all nodal force    values to zeron
c     [rese]t  <boun>d: Set all boundary condition codes to zero

      elseif(i.eq.31) then

        if    (pcomp(c2,'angl',4)) then          ! Boundary angles
          call pzero (hr(np(45)     ),numnp)
        elseif(pcomp(c2,'disp',4)) then          ! Nodal displacements
          call pzero (hr(np(27)+nneq),nneq)
        elseif(pcomp(c2,'forc',4)) then          ! Nodal forces
          call pzero (hr(np(27)     ),nneq)
        else                                     ! Boundary codes
          call pzeroi(mr(np(31)+nneq),nneq)
        endif

c     [sblo,nsblk] - input surface blocks to generate 3-d meshes

      elseif(i.eq.32) then
        cc    = yyy(1:4)
        errck = vinput(yyy(16:30),15,td,1)
        nsblk = max(1,nint(td(1)))

c       Allocate arrays for 3-d generation from surface mesh

        setvar = palloc(111,'TEMP1', 9*nsblk,  1)
        setvar = palloc(112,'TEMP2', 9*nsblk,  1)
        setvar = palloc(113,'TEMP3', 5*nsblk,  1)
        setvar = palloc(114,'TEMP4',27*nsblk,  2)
        setvar = palloc(115,'TEMP5', 9*nsblk,  2)
        setvar = palloc(116,'TEMP6',27*nsblk,  2)

        call sblkgn(nsblk,ndm,nen,nen1,hr(np(43)),mr(np(33)),
     &              mr(np(111)),mr(np(112)),mr(np(113)),hr(np(114)),
     &              hr(np(115)),hr(np(116)),prt,prth)

c       Delete arrays used to generate mesh

        setvar = palloc(116,'TEMP6', 0,  2)
        setvar = palloc(115,'TEMP5', 0,  2)
        setvar = palloc(114,'TEMP4', 0,  2)
        setvar = palloc(113,'TEMP3', 0,  1)
        setvar = palloc(112,'TEMP2', 0,  1)
        setvar = palloc(111,'TEMP1', 0,  1)

c     [curv]e,<bound,disp,load>: Specify conditions on curve.

      elseif(i.eq.33) then
        curfl = .true.
        call setext('curv',ncurv, fext, .true.)
        if(.not.pcomp(c2,'set',3)) then
          c2 = 'add'                 ! default mode 'add' for edge b.c.
        endif
        call plinka(fext,c2,'   ')

c     [rota] - set flags for rotational transformations
c     [rota,<opti,rotopt>] - set rotation option

      elseif(i.eq.34) then
        if(pcomp(yyy(16:19),'opti',4)) then
          errck = vinput(yyy(31:45),15,td,1)
          rotopt = max(1,nint(td(1)))
          write(iow,2014) rotopt

        elseif (.not.frotas) then
          frotas = .true.
          setvar = palloc(81,'MO   ',numnp*2 ,1)
          setvar = palloc(83,'MT   ',numnp   ,2)
          setvar = palloc(82,'MR   ',numnp*54,2)
        endif

        call protin(prt,prth)

c     [setn] - define second mesh for director input
c            & set initial normal quantities

      elseif(i.eq.35) then
        call gendir(hr(np(82)),' Director Nodes',prt,prth,error,.true.)

c     [setr] -set initial rotation transformations from nodal directors

      elseif(i.eq.36) then
        call setrot(hr(np(43)),mr(np(81)),hr(np(82)),hr(np(83)),
     &              numnp,rotopt)

c       Read input

        call genvec(9,54,hr(np(82)),' Transforms ',prt,prth,error,
     &             .true.)

c       Set default rotation matrices

        call defrot(hr(np(82)),hr(np(27)),mr(np(31)+nneq),numnp,ndf)

c       Initialize all configurations

        call pmove(hr(np(82)),hr(np(82)+numnp*9 ),numnp*9)
        call pmove(hr(np(82)),hr(np(82)+numnp*18),numnp*9)
        call pmove(hr(np(82)),hr(np(82)+numnp*45),numnp*9)

c     [btra] - transition mesh patch

      elseif(i.eq.37) then
        if(iii.lt.0) write(iow,3002)
        call blktra(ndm,nen1,hr(np(43)),mr(np(33)),prt,prth)

c     [fpro],<set,add> - Proportional load number specification

      elseif(i.eq.38) then
        if(pcomp(c2,'add',3)) then
          setvar = palloc(151,'USER1', ndf*numnp,1)
          call pzeroi(mr(np(151)),ndf*numnp)
          fp(1) = np(151)
          fp(2) = np(29)
        else
          fp(1) = np(29)
        endif

        call genint(ndf,mr(fp(1)),ndf,numnp,'P r o p.  L o a d  N o s.',
     &              '-dof',prt,prth,error,1)

        if(pcomp(c2,'add',3)) then
          do i = 0,ndf-1
            jj = 0
            do ii = i,ndf*numnp-1,ndf
              jj = jj + 1
              if(mr(ii+fp(1)).gt.0) then
                if(mr(ii+fp(2)).gt.0            .and.
     &             mr(ii+fp(2)).ne.mr(ii+fp(1))) then
                   write(iow,3006) jj,i+1,mr(ii+np(fp(2))),
     &                                    mr(ii+np(fp(1)))
                   write(ilg,3006) jj,i+1,mr(ii+np(fp(2))),
     &                                    mr(ii+np(fp(1)))
                else
                  mr(ii+fp(2)) = mr(ii+fp(1))
                endif
              endif
            end do ! ii
          end do ! i
          setvar = palloc(151,'USER1', 0,1)
        endif

c     [cpro]  coordinate proportional load number specification

      elseif(i.eq.39) then
        cprfl = .true.
        call setext('cpro',nprof, fext, .true.)
        call plinka(fext,c2,'   ')

c     [regi,nreg]  set region number: all

      elseif(i.eq.40) then
        cc    = yyy(1:4)
        errck = vinput(yyy(16:30),15,td,1)
        if(errck) go to 110
        nreg  = nint(td(1))
        mxreg = max(mxreg,nreg)
        write(iow,2001) nreg,mxreg
        if(ior.lt.0) then
          write(*,2001) nreg,mxreg
        endif

c     [tran],<inc> - Specify coordinate transformation array

c     Incremental rotation update: tr_new = tinc * tr_old
c                                  xr_new = tinc * xr_old + xr_inc
      elseif(i.eq.41) then
        call ptranf(tx(2),prt)

c     [damp],<set,add> - specify nodal damping quantities
c     [mass],<set,add> - specify nodal mass quantities
c     [stif],<set,add> - specify nodal stiffness quantities

      elseif(i.eq.42 .or. i.eq.43 .or. i.eq.44) then
        if(.not.nmfl) then
          setvar = palloc(88,'NSTI  ',ndf*numnp,2)
          setvar = palloc(87,'NMAS  ',ndf*numnp,2)
          setvar = palloc(86,'NDAM  ',ndf*numnp,2)
          nmfl = .true.
        endif
        if(pcomp(c2,'add',3)) then
          setvar = palloc(151,'USER1', ndf*numnp,2)
          call pzero(hr(np(151)),ndf*numnp)
          fp(1) = np(151)
          fp(2) = np(i+44)
        else
          fp(1) = np(i+44)
        endif
        if(i.eq.42) then
          call genvec(ndf,ndf,hr(fp(1)),' Nodal Dampers',prt,prth,
     &                error,.false.)
        elseif(i.eq.43) then
          call genvec(ndf,ndf,hr(fp(1)),' Nodal Masses' ,prt,prth,
     &                error,.false.)
        elseif(i.eq.44) then
          call genvec(ndf,ndf,hr(fp(1)),' Nodal Springs',prt,prth,
     &                error,.false.)
        endif
        if(pcomp(c2,'add',3)) then
          do ii = 0,ndf*numnp-1
            hr(ii+fp(2)) = hr(ii+fp(2)) + hr(ii+fp(1))
          end do ! ii
          setvar = palloc(151,'USER1', 0,2)
        endif

c     [csur] - surface loading by coordinates

      elseif(i.eq.45) then
        surfl = .true.
        call setext('csur',nsurf, fext, .true.)
        call plinka(fext,c2,'   ')

c     [ereg] - set element regions

      elseif(i.eq.46) then
        call genint(1,mr(np(33)+nen1-2),nen1,numel,'R e g i o n  N o s',
     &              '-regn',prt,prth,error,2)
        do j = nen1-2,numel*nen1,nen1
          mxreg = max(mxreg,mr(np(33)+j))
        end do ! j

c     [reac,file] - input reactions from 'file'

      elseif(i.eq.47) then
        ii = index(xxx,' ')
        if(index(xxx,',').gt.0) then
          ii = min(ii,index(xxx,','))
        elseif(index(xxx,'=').gt.0) then
          ii = min(ii,index(xxx,'='))
        endif
        reafi = xxx(ii+1:ii+127)
        ii = index(reafi,' ')
        if(index(reafi,',').gt.0) then
          ii = min(ii,index(reafi,','))
        elseif(index(reafi,'=').gt.0) then
          ii = min(ii,index(reafi,'='))
        endif
        reafi = reafi(1:ii-1)
        inquire(file=reafi,exist=reafl)
        if(.not.reafl) then
          if(ior.lt.0) then
            write(*,3003) reafi
          else
            if(iii.ge.0) write(ilg,3003) reafi
            write(iow,3003) reafi
            call plstop()
          endif
        else
          if(prt) then
            write(iow,2015) reafi(1:ii-1)
          endif
        endif

c     [manu],hlplev - set Manual help options level

      elseif(i.eq.48) then
        cc    = yyy(1:4)
        errck = vinput(yyy(16:30),15,td,1)
        if(errck) go to 110
        hlplev = max(-1,min(3,nint(td(1))))

c     [body] - compute internal forces, will be stored in nodal loads

      elseif(i.eq.49) then
        intfl = .true.
        call setext('body',nintf, fext, .true.)
        call plinka(fext,c2,'   ')

c     [glob]al - set global parameters

      elseif(i.eq.50) then
        call global()

c     [shif]t:<x0,y0,z0> - origin for polar/spherical conversions

      elseif(i.eq.51) then
51      errck = pinput(x0,3)
        if(errck) go to 51
        if(prt) then
          write(iow,2002) (x0(i),i=1,ndm)
          if(ior.lt.0) then
            write(*,2002) (x0(i),i=1,ndm)
          endif
        endif

c     [disp]l,<set,add> - data input

      elseif(i.eq.52) then
        if(pcomp(c2,'add',3)) then
          setvar = palloc(151,'USER1', ndf*numnp,2)
          call pzero(hr(np(151)),ndf*numnp)
          fp(1) = np(151)
          if(ldflg) then
            fp(2) = np(26) + ndf*numnp
            lddis = .true.
          else
            fp(2) = np(27) + ndf*numnp
          endif
        else
          if(ldflg) then
            fp(1) = np(26) + ndf*numnp
            lddis = .true.
          else
            fp(1) = np(27) + ndf*numnp
          endif
        endif
        call genvec(ndf,ndf,hr(fp(1)),' Displacements',
     &              prt,prth,error,.false.)
        if(pcomp(c2,'add',3)) then
          do ii = 0,ndf*numnp-1
            hr(ii+fp(2)) = hr(ii+fp(2)) + hr(ii+fp(1))
          end do ! ii
          setvar = palloc(151,'USER1', 0,2)
        endif

c     [edis] set edge displacement constraints

      elseif(i.eq.53) then
        edifl = .true.
        call setext('edis',nedis, fext, .true.)
        if(.not.pcomp(c2,'set',3)) then
          c2 = 'add'                 ! default mode 'add' for edge displ
        endif
        call plinka(fext,c2,'   ')

c     [cdis] set coordinate nodal forces - based on coordinates

      elseif(i.eq.54) then
        disfl = .true.
        call setext('cdis',ndisf, fext, .true.)
        if(.not.pcomp(c2,'set',3)) then
          c2 = 'add'                 ! default mode 'add' for edge angle
        endif
        call plinka(fext,c2,'   ')

c     [debu]g,<on,off>  Activate,deactivate debug option

      elseif(i.eq.55) then
        if(pcomp(c2,'off',3)) then
          debug = .false.
        else
          debug = .true.
        endif

c     [side] - for blending interpolations

      elseif(i.eq.56) then
        if(iii.ge.0) then
          numsd  = max(numsd,1)
          setvar = palloc ( 162, 'BSIDE', numsd*isd, 1)
          call psides(mr(np(162)),side,isd,prt,prth,1)
          setvar = palloc ( 162, 'BSIDE', numsd*isd, 1)
        else
          write(*,3001) 'SIDEs'
        endif

c     [face] - for blending interpolations

      elseif(i.eq.57) then
        if(iii.ge.0) then
          setvar = palloc ( 165, 'BFACE', numfc*isd, 1)
          call psides(mr(np(165)),face,isd,prt,prth,2)
        else
          write(*,3001) 'SIDEs'
        endif

c     [snod]e - for blending interpolations

      elseif(i.eq.58) then
        if(iii.ge.0) then
          numsn  = max(numsn,1)
          setvar = palloc ( 161, 'BNODE', numsn*3, 2)
          call pnodes(hr(np(161)),ndm,prt,prth)
          snofl  = .true.
        else
          write(*,3001) 'SNODes'
        endif

c     [blen]ding interpolations (Delayed mesh generation feature)

      elseif(i.eq.59) then
        if(iii.ge.0) then
          nblend = nblend + 1
          numbd  = max(numbd,nblend)
          setvar = palloc ( 163, 'BTRAN', numbd*12          , 2)
          setvar = palloc ( 164, 'BLEND', numbd*ibn         , 1)
          setvar = palloc ( 166, 'BNILR', numbd*max(1,mxilr), 1)
          call pblend(hr(np(163)),mr(np(164)),mr(np(166)),nblend,ibn,
     &                ndm,prt,prth)
        else
          write(*,3001) 'BLENDs'
        endif

c     [move]  relocate coordinates based on some specified locations

      elseif(i.eq.60) then
        call pcormv(hr(np(43)),ndm,numnp,prt,prth)

c     [rigi,nrigid,rbcen(nrigid),rbx0(i,nrigid)]
c     [moda,nrigid,rbcen(nrigid),rbx0(i,nrigid)]

c     Set rigid body number/type
c     rbtype(nrigid)  = 0   :  rigid-flex (default)
c                     = 1   :  rigid-modal
c     rbcen(nrigid)   = 0   :  At center of mass
c                     = >0  :  At rbx0
c     rbx0(ii,nrigid) = x(i):  reference location for dofs

      elseif(i.eq.61 .or. i.eq.62) then
        cc     = yyy(1:4)
        errck  = vinput(yyy(16:105),90,td,6)
        if(errck) go to 110
        nrigid = nint(td(1))

c       Return to non-rigid mode
        if(nrigid.eq.-1) then
          nrigid = 0
          write(iow,2006)
          if(ior.lt.0) write(iow,2006)
          go to 100

c       Set rigid body number if none is given by user

        elseif(nrigid.eq.0) then
          nrigid = nrbody + 1
        end if
        nrbody = max(nrbody,nrigid)

c       Set rigid/modal indicator

        if(i.eq.61) then
          rbtype(nrigid) = 0
        else
          rbtype(nrigid) = 1
        endif

c       Set location for dofs

        rbcen(nrigid)  = nint(td(2))
        if(rbcen(nrigid).gt.0) then
          do ii = 1,ndm
            rbx0(ii,nrigid) = td(ii+2)
          end do ! ii
        endif

c       Set modal body number

        if(rbtype(nrigid).eq.1) then
          nmbody = nmbody + 1
          modbod(nmbody) = nrigid
          write(iow,2003) nrigid
          if(ior.lt.0) then
            write(*,2003) nrigid
          endif
        else
          write(iow,2004) nrigid
          if(ior.lt.0) then
            write(*,2004) nrigid
          endif
        end if
        if(rbcen(nrigid).gt.0) then
          write(iow,2005) (rbx0(ii,nrigid),ii=1,ndm)
          if(ior.lt.0) then
            write(*,2005) (rbx0(ii,nrigid),ii=1,ndm)
          endif
        end if

c     [flex]ible - set element types to deformable (non-rigid)

      elseif(i.eq.63) then
        nrigid = 0
        write(iow,2006)
        if(ior.lt.0) write(iow,2006)

c     [base]  proportional load number specification

      elseif(i.eq.64) then
        if(np(125).eq.0) then
          setvar = palloc(125,'NUBAS',ndf*numnp,1)
        endif
        call genint(ndf,mr(np(125)),ndf,numnp,
     &             'B a s e   P a t t e r n s',
     &            '-dof',prt,prth,setvar,1)

c     [epro] set edge proportional load numbers

      elseif(i.eq.65) then
        eprfl = .true.
        call setext('epro',nepro, fext, .true.)
        if(.not.pcomp(c2,'set',3)) then
          c2 = 'add'                 ! default mode 'add' for edge prop
        endif
        call plinka(fext,c2,'   ')

c     [mpro],<set,add> -  mass proportional load number specification

      elseif(i.eq.66) then
        if(pcomp(c2,'add',3)) then
          setvar = palloc(151,'USER1', ndf*numnp,1)
          call pzeroi(mr(np(151)),ndf*numnp)
          fp(1) = np(151)
          fp(2) = np(29) + nneq
        else
          fp(1) = np(29) + nneq
        endif
        call genint(ndf,mr(fp(1)),ndf,numnp,
     &             'M a s s  P r o p.  L o a d','-dof',
     &              prt,prth,error,1)
        if(pcomp(c2,'add',3)) then
          do ii = 0,ndf*numnp-1
            fp(3) = ii + ii
            mr(ii+fp(2)) = mr(ii+fp(2)) + mr(fp(3))
          end do ! ii
          setvar = palloc(151,'USER1', 0,1)
        endif

c     [loop,#] - Loop start

      elseif(i.eq.67) then
        call ploops(lp_in,tx(2), 1)

c     [next] - Loop end

      elseif(i.eq.68) then
        call ploops(lp_in,tx(2), 2)

c     [file] - File inputs: Open new input file: lp_file

      elseif(i.eq.69) then
        lp_file = tx(2)
        lopen   = .true.
        lp_lun  = icl
        do while(lopen)
          lp_lun       = lp_lun + 1
          inquire(unit = lp_lun, opened = lopen)
        end do ! while
        ior       = lp_lun
        open(unit = ior, file = lp_file, status = 'old')

c     [cdam] set coordinate nodal dampers   - based on coordinates
c     [cmas] set coordinate nodal masses    - based on coordinates
c     [csti] set coordinate nodal stiffness - based on coordinates

      elseif(i.eq.70 .or. i.eq.71 .or. i.eq.72) then
        if(.not.nmfl) then
          setvar = palloc(88,'NSTI  ',ndf*numnp,2)
          setvar = palloc(87,'NMAS  ',ndf*numnp,2)
          setvar = palloc(86,'NDAM  ',ndf*numnp,2)
          nmfl = .true.
        endif

c       [cdam] set coordinate nodal dampers - based on coordinates

        if(i.eq.70) then
          damfl = .true.
          call setext('cdam',ndamf, fext, .true.)
          call plinka(fext,c2,'   ')

c       [cmas] set coordinate nodal masses - based on coordinates

        elseif(i.eq.71) then
          masfl = .true.
          call setext('cmas',nmasf, fext, .true.)
          call plinka(fext,c2,'   ')

c       [csti] set coordinate nodal stiffness - based on coordinates

        elseif(i.eq.72) then
          stifl = .true.
          call setext('csti',nstif, fext, .true.)
          call plinka(fext,c2,'   ')
        endif

c     [ebas] - set edge base proportional loads

      elseif(i.eq.73) then
        if(np(125).eq.0) then
          setvar = palloc(125,'NUBAS',ndf*numnp,1)
        endif
        ebsfl = .true.
        call setext('ebas',nebas, fext, .true.)
        if(.not.pcomp(c2,'set',3)) then
          c2 = 'add'                 ! default mode 'add' for edge b.c.
        endif
        call plinka(fext,c2,'   ')

c     [cbas] set coordinate base proportional load - based on coords

      elseif(i.eq.74) then
        if(np(125).eq.0) then
          setvar = palloc(125,'NUBAS',ndf*numnp,1)
        endif
        basfl = .true.
        call setext('cbas',nbasf, fext, .true.)
        call plinka(fext,c2,'   ')

c     [eule]r  - angles for nodal transformations

      elseif(i.eq.75) then

        if(np(242).eq.0) then
          setvar = palloc(242,'EULER',numnp*3 , 2)
          setvar = palloc(243,'LEULR',nen*3   , 2)
          eulerfl = .true.
        endif
        call genvec(3,3,hr(np(242)),' Euler',prt,prth,error,.false.)

c     [ceul]er  - coordinate angles for nodal transformations

      elseif(i.eq.76) then

        if(np(242).eq.0) then
          setvar = palloc(242,'EULER',numnp*3 , 2)
          setvar = palloc(243,'LEULR',nen*3   , 2)
          eulerfl = .true.
          call setext('ceul',neule, fext, .true.)
          call plinka(fext,c2,'   ')
        endif

c     [rfor]ce -- Specify Radial follower nodal FORces

      elseif(i.eq.77) then
        if(.not.nffl) then
          setvar = palloc(256,'NFORC',ndf*numnp,2)
          nffl = .true.
        endif
        call genvec(ndf,ndf,hr(np(256)),' Nodal Follower Forces',
     &              prt,prth,error,.false.)

c     [lfor]ce -- Specify Line nodal FORces

      elseif(i.eq.78) then

        lfrfl = .true.
        call setext('lfor',nforl, fext, .true.)
        call plinka(fext,c2,'   ')

c     [load,<prop, prop_num> - Load group with proportional load number
c     [load]                 - Load group with total proportional load
c     [load end]             - Load group end
c     [load set]             - Load table set

      elseif(i.eq.79) then

c       Start of load group

        if(pcomp(tx(2),'prop',4) .or. pcomp(tx(2),'    ',4)) then

c         If previous load open, close and save

          if(ldflg) then
            if(prt) write(iow,2013) ldnum
            call pldend(ldfor,lddis,ldflg)
          endif

          call setval(tx(3),15,td(1))
          ldflg = .true.
          spflg = .false.
          ldfor = .false.
          lddis = .false.
          if(pcomp(tx(2),'prop',4)) then
            ldprp = nint(td(1))
          else
            ldprp = 0
          endif
          ldnum = ldnum + 1
          setvar = palloc( 26,'DR   ',nneq*2,2) ! Used to store F/U
          if(prt) write(iow,2012) ldnum,ldprp

c         Check that load table is large enough

          if(ldtot.le.ldnum) then
            ldtot = ldnum
            setvar = palloc(265,'LDTAB',ldtot*12, 1)
          endif

c         Store proportional number load table for force, displ & spin

          call pldtabl(mr(np(265)),ldnum, ldprp, 1, 3) ! Force Prop Ld.
          call pldtabl(mr(np(265)),ldnum, ldprp, 2, 3) ! Displ Prop Ld.

c       End of load group

        elseif(pcomp(tx(2),'end',3)) then

          if(prt) write(iow,2013) ldnum
          call pldend(ldfor,lddis,ldflg)
          spflg = .false.

c       Set load table

        elseif(pcomp(tx(2),'set',3)) then
          ldtot = nint(td(1))
          setvar = palloc(265,'LDTAB',ldtot*12, 1)
        endif

c     [swee]p

      elseif(i.eq.80) then

        errck = tinput(tx(1),1,td,0)

        if(prt) then
          write(iow,2009) head,' 2-d File Name: ',tx(1)
        endif

        inquire(file=tx(1),exist=eflg)

        if(eflg) then

          open(unit=99,file=tx(1))

          iorold = ior
          ior    = 99

          tx(1) = 'start'
          do while( .not.pcomp(tx(1),'end',3))

            errck = tinput(tx(1),1,td,6)

            if(pcomp(tx(1),'feap',4)) then

              errck = tinput(tx(1),0,td,6)
              numn2 = nint(td(1))
              nume2 = nint(td(2))
              numm2 = nint(td(3))
              ndm2  = nint(td(4))
              ndf2  = nint(td(5))
              nen2  = nint(td(6))
              setvar = palloc(151,'USER1 ',numn2*2       ,2)
              setvar = palloc(152,'USER2 ',nume2*(nen2+1),1)

            elseif(pcomp(tx(1),'coor',4)) then

              do nn = 1,numn2
                errck = pinput(td,4)
                n1    = 2*nint(td(1)) - 2
                n2    = n1 + 1
                hr(np(151)+n1)  = td(3)
                hr(np(151)+n2)  = td(4)
              end do ! nn

            elseif(pcomp(tx(1),'elem',4)) then

              do nn = 1,nume2
                errck = pinput(td,nen2+3)
                point = np(152) + (nen2+1)*(nint(td(1)) - 1) - 1
                mr(point+nen2+1) = nint(td(3))
                do jj = 1,nen2
                  mr(jj+point) = nint(td(jj+3))
                end do ! jj
              end do ! nn

            endif
          end do ! while
          close(99)
          ior = iorold
        else
          write(iow,*) 'FILE: ',tx(1),' does not exist.'
          call plstop()
        endif

c       Input Number of segments and axis of rotation

        errck = tinput(tx(1),1,td,3)
        if(pcomp(tx(1),'cart',4)) then
          nstyp = 1
        else
          nstyp = 2
        endif
        nseg = max(1,nint(td(1)))
        ndir = max(1,min(2,nint(td(2))))

c       Input center location

        errck = tinput(tx(1),1,xd,3)

c       Loop input on segment master nodes

        nn    = 0
        td(1) = 1.d0
        if(prt) then
          write(iow,2010) swd(nstyp),nseg,ndir, xd
        endif
        do while(nint(td(1)).gt.0)
           errck = pinput(td,4)
           if(nint(td(1)) .gt.0) then
             nn = nn+ 1
             is(nn)   = nint(td(1))
             xx(1,nn) = td(2)
             xx(2,nn) = td(3)
             xx(3,nn) = td(4)
             if(prt) then
               write(iow,2011) nn,(xx(jj,nn),jj=1,3)
             endif
           endif
        end do ! while

c       Generate the segments

        call interp1(nseg,xx,1,is,nn,3,shp, xp)

c       Generate x-y-z coordinates

        call pgen3dx(hr(np(43)),mr(np(190)),hr(np(151)),
     &               xd,xp, numn2,nstyp,nseg,ndir, nnmax, prt)

c       Generate elements

        call pgen3de(hr(np(43)),mr(np(33)),mr(np(152)),
     &               nen2+1,nume2,numn2,nseg,prt)

        setvar = palloc(152,'USER2',0,1)
        setvar = palloc(151,'USER1',0,2)

c     [spin]  rotate around axis

      elseif(i.eq.81) then

        spnum = spnum + 1
        spflg = .true.
        if(spnum.gt.sptot.or. np(268).eq.0) then
          sptot  = max(sptot,spnum)
          setvar = palloc(268,'SPINS',sptot*12,2)
        endif
        call pspinin(hr(np(268)),spnum,prt)
        if(ldflg) then
          lddis = .true.
        endif

c     [peri]odic case

      elseif(i.eq.82) then

c       Set problem type

        if(pcomp(tx(2),'ther',4)) then        ! Thermal
          prtype = 1
          write(iow,2021)
        elseif(pcomp(tx(2),'cauc',4)) then    ! Cauchy mechanical
          prtype = 2
          write(iow,2022)
        elseif(pcomp(tx(2),'coup',4)) then    ! Coupled thermo-mech
          prtype = 3
          write(iow,2023)
        elseif(pcomp(tx(2),'off' ,3)) then    ! Turn off
          prtype = 0
          write(iow,2025)
        endif

c       Check for file inputs, input periodic state

        if(pcomp(tx(3),'file',4)) then
          filflg   = .true.
          errck    = tinput(tx(8),1,td,0)
          hillfile = xxx(1:128)
          write(iow,2026) hillfile
          if(prtype.gt.0) perflg   = .true.
        else
          filflg = .false.
          if(prtype.gt.0) call periodic(prt)
        endif

c     [expl]icit-implicit element sets

      elseif(i.eq.83) then

        call psetexim(mr(np(33)),nen,nen1,numel)

c     [tria]d - 3-d boundary triads

      elseif(i.eq.84) then

        if(np(274).eq.0) then
          setvar = palloc(274,'TRIAD',numnp*9 , 2)
          setvar = palloc(275,'LTRIA',nen*9   , 2)
          triadfl = .true.
        endif
        call genvec(9,9,hr(np(274)),' Triad',prt,prth,error,.false.)

c       Mark unspecified points

        do j = 1,9*numnp,9
          point = np(274) + j - 1
          eflg  = .true.
          do jj = 0,8
            if(hr(point+jj).ne.0.0d0) then
              eflg = .false.
              exit
            endif
            if(eflg) then
              hr(point) = -100.0d0
            endif
          end do ! jj
        end do !

c     [tbou]nd - Taylor boundary conditions

      elseif(i.eq.85) then

        call ptay_bc(numnp,ndf,mr(np(31)+nneq))

c     [*nod] - Number to add to all input and generated nodes
c     [*ele] - Number to add to all input and generated elements

      elseif(i.eq.nstar .or. i.eq.estar) then
        j = index(record,'=')
        if(j.eq.0) then
          j = index(record,',')
        endif
        xxx(1:75) = ' '
        ii = 0
        do jj = j+1,75
          if(record(jj:jj).ne.' ') then
            ii = ii + 1
            xxx(ii:ii) = record(jj:jj)
          endif
        end do ! jj

c       Set flag to parse expression

        flgco = coflg
        coflg = .true.
        call setval(xxx,75,td(1))
        coflg = flgco

c       Set *node

        if(i.eq.nstar) then
          starnd = nint(td(1))
          write(iow,2007) starnd
          if(ior.lt.0) then
            write(*,2007) starnd
          endif

c       Set *element

        else
          starel = nint(td(1))
          write(iow,2008) starel
          if(ior.lt.0) then
            write(*,2008) starel
          endif
        endif

c     [mesn] -> user defined mesh inputs

      elseif(i.gt.list-numesh) then
        n   = i + numesh - list
        uct = wd(i)
        call umshlib(n,tx,prt)

      endif

      go to 100

c     Formats

2001  format(/' -> Region Number:',i4,' Maximum:',i4)
2002  format(' Coordinate transformation origin set to:'/
     &       '   x0 =',1p,e12.4:,'  y0 =',1p,e12.4:,'  z0 =',1p,e12.4)
2003  format(' -> Rigid body:',i4,' with modal deformations')
2004  format(' -> Rigid body:',i4,' initiated')
2005  format('    Location for Rigid Body Unknowns'/
     &   10x,'x0 =',1p,1e12.4:,'  y0 =',1p,1e12.4:,'  y0 =',1p,1e12.4)
2006  format(' -> Flexible elements')
2007  format(' -> Number to be added to all nodal   inputs =',i8)
2008  format(' -> Number to be added to all element inputs =',i8)

2009  format(/1x,20a4//'   3-D   M e s h   f r o m   2-D   M e s h'//
     &       10x,a,a)

2010  format(/10x,' Spline Definition for 3-d mesh'//
     &        11x,a,' sweep'/
     &        10x,' Generate ',i4,' segments'//
     &        10x,' Rotate on',i4,' axis'//
     &        10x,' Center at:',1p,3e12.4//
     &        10x,' Node     1-Coord     2-Coord     3-Coord')

2011  format(i15,1p,3e12.4)
2012  format(/5x,'Start Load Set',i4/10x,' Proportinal Load =',i4)
2013  format(/5x,'End Load Set',i4)

2014  format(/5x,'Director Rotation Update Option',
     &       /10x,'Type =',i3)

2015  format(/5x,'Reaction forces from file: ',a)

2021  format(/5x,'Thermal periodic boundary condition')
2022  format(/5x,'Cauchy mechanical periodic boundary condition')
2023  format(/5x,'Coupled thermomechanical periodic boundary condition')
2025  format(/5x,'Periodic boundary condition off')
2026  format(/5x,'File input periodic boundary condition'/
     &       10x,' File =',a)

3001  format(' *ERROR* PMESH: Cannot regenerate ',a)
3002  format(' *WARNING* Initial node/element numbers necessary to'
     &      ,' use BLOCk in solution mode.')
3003  format(' *ERROR* PMESH: File:',a,' does not exist.')
3004  format(' *ERROR* PMESH: No SAVE,END statement for input data.')
3005  format(' *ERROR* PMESH: No SNODes were input for blending.')
3006  format(' *ERROR* PMESH: Attempt to define different proportional',
     &        ' load number to:'/5x,'Node =',i8,': DOF =',i3,': Old =',
     &         i4,': New =',i4)

4000  format('  INPUT FILE NAME: ',a//
     &       '  MESH INPUT RECORDS'/'  ------------------')
4001  format(5x,a,' ',a)

      end
