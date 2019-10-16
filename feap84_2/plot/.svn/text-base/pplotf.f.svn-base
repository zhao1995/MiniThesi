c$Id:$
      subroutine pplotf(lci,ct,prop)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add 'strin'g plot caption command                09/02/2007
c       2. Change arguments to call of 'plotelm'            13/02/2007
c       3. Correct print for 'rfor'/'disp'; Change label    10/06/2007
c          107 to 108 (correct entry position)
c       4. Add 'acti' and 'deac' from plot state            11/07/2007
c       5. Set 'evtyp' from 'maci' for eigenvalue labels    28/08/2007
c       6. Plot current mesh after eyes change.             07/09/2007
c          Set hide to false after new 'pers' command.
c       7. Dimention reaction array to max(neq,nneq)        05/10/2007
c       8. Correct setting of 'flux' number for k1=0        05/11/2007
c       9. Correct calling arguments to 'pelink'            17/02/2008
c      10. Add call to ploa1 to set rlnew (and save it)     17/04/2008
c      11. Add print for 'cart' and 'pers' views            23/04/2008
c      12. Output plot prompt for rank 0 processor only     07/02/2009
c      13. Change np(53) to np(43) on call to plcyls        19/02/2009
c      14. Add 'rffl' flag for file reaction outputs        22/02/2009
c          Use also for force plots.
c      15. Set 'pltmfl' to true for reactions               22/03/2009
c      16. Add np(31)+nneq array to setfor call             29/05/2009
c      17. Use 'fpltfl' to control multiple openings of     30/08/2009
c          plot windows.  Set all 'off' options to be 1.0   09/10/2009
c      18. Check for existing velo and acce for contours    30/11/2009
c      19. Add .true. to scalev call                        09/10/2010
c          Add scale by macd for contours
c      20. Change 'nix' to npix in references to ix array   13/10/2010
c      21. Remove call to pload for reaction prints         18/10/2010
c      22. Use 'plfl' to allocate projection arrays         29/10/2010
c      23. Replace k5 by pl(5) in call to prax              06/03/2011
c          Replace all fp(*) by pl(*) and use pl_int.h
c      24. Change k4 to pl(1) for [splo]                    17/04/2011
c      25. Add JPEG option for screen dumps                 03/05/2011
c      26. For 'mate' command use
c          Introduce 'plix' for some uses of 'npix'         16/11/2011
c      27. Set storage for 'HDNP' and 'HSELM' hist plots    05/01/2012
c      28. Remove arg 3 and arg 4 from call to perspz       09/01/2012
c      29. Change color of undeformed mesh to cyan          01/06/2012
c      30. Replace perspz statments by call to phide        20/07/2012
c      31. Pass nne to pndata( )                            09/01/2013
c      32. Remove allocation of FCZM and np(55) use         05/06/2013
c      33. Add plot of shell 'tria'ds                       22/06/2013
c-----[--.----+----.----+----.-----------------------------------------]
c     P l o t   C o n t r o l   R o u t i n e   F o r   F E A P
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Driver for plot commands.

c      Inputs:
c         lci       - Plot command option: if blank interactive inputs
c         ct(3)     - Plot command parameters
c         prop      - Proportional load value

c      Outputs:
c         none      - Plot outputs to screen/file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'comblk.h'
      include  'arclel.h'
      include  'arclei.h'
      include  'arcler.h'
      include  'cblend.h'
      include  'cdata.h'
      include  'cdat1.h'
      include  'chdata.h'
      include  'comfil.h'
      include  'complx.h'
      include  'counts.h'
      include  'ddata.h'
      include  'elacts.h'
      include  'eldata.h'
      include  'eldatp.h'
      include  'elpdat.h'
      include  'endata.h'
      include  'evdata.h'
      include  'fdata.h'
      include  'hdatam.h'
      include  'hlpdat.h'
      include  'idptr.h'
      include  'iodata.h'
      include  'iofile.h'
      include  'mxsiz.h'
      include  'ndata.h'
      include  'part0.h'
      include  'pbody.h'
      include  'pconstant.h'
      include  'pdata1.h'
      include  'pdata2.h'
      include  'pdata3.h'
      include  'pdata4.h'
      include  'pdatap.h'
      include  'pdatas.h'
      include  'pdatay.h'
      include  'pdatps.h'
      include  'pdatxt.h'
      include  'pfeapb.h'
      include  'plcapt.h'
      include  'plclip.h'
      include  'plflag.h'
      include  'pltfac.h'
      include  'pointer.h'
      include  'ppers.h'
      include  'prange.h'
      include  'prflag.h'
      include  'print.h'
      include  'prmptd.h'
      include  'prstrs.h'
      include  'pview.h'
      include  'region.h'
      include  'sdata.h'
      include  'setups.h'
      include  'strnum.h'
      include  'tdata.h'
      include  'umac1.h'
      include  'wdata.h'

      include  'pl_int.h'

      character lci*4,lct*4,button*1,tx(2)*15
      logical   errv,labl,odefalt,outf,pflag,ppr,praxf,symf,timf
      logical   pcomp,pinput,tinput,vinput,palloc, setvar,plmesh, tr,fa
      logical   rffl
      integer   i,icol,icp,ijump,jj,k1,k2,k3,k4,k5,kr,l,list, intvl
      integer   nsinc
      integer   ni,nfaco, iln(2)
      real*8    dumv,es,xm,zcol, prop, swang, rlold
      real*8    ct(3),ss(4),tt(4),td(4)

c     Set number commands

      integer    ncomd
      parameter (ncomd = 124)

      character wd(ncomd)*4
      integer   ed(ncomd),xd(ncomd)

      save

      data      list/ncomd/

c     Plot data command list

      data wd/'fram','wipe','fact','cent','cart','line','symm','cont',
     1        'outl','load','mesh','stre','node','boun','elem','zoom',
     2        'colo','fill','text','size','cvar','eigv','bord','scal',
     3        'axis','pers','disp','show','hide','prin','nopr','defo',
     4        'unde','velo','acce','post','reac','eige','mate','back',
     5        'clip','titl','mark','refr','pick','capt','pbou','pfor',
     6        'pnod','quad','real','imag','eyes','dofs','estr','prof',
     7        'prax','pair','clea','pstr','dplo','splo','manu','prom',
     8        'defa','scre','pdis','pele','proj','labe','nola','snod',
     9        'psno','exno','xypl','wind','logo','time','bplo','rang',
     a        'nora','rect','cyli','sphe','full','nofu','uplo','jint',
     b        'regi','norm','flux','cwir','vwir','awir','swir','pwir',
     c        'ewir','inte','swee','tors','edef','renu','elpl','ndat',
     d        'vplo','aplo','rfor','stri','acti','deac','jpeg','hist',
     e        'stra','tria',
     u        'plt1','plt2','plt3','plt4','plt5','plt6','plt7','plt8',
     u        'plt9','plt0'/

      data ed/    1,     0,     1,     2,     0,     1,     1,     0,
     1            0,     0,     0,     0,     0,     0,     0,     0,
     2            1,     0,     1,     1,     1,     0,     1,     1,
     3            0,     0,     0,     1,     0,     1,     1,     0,
     4            0,     0,     0,     0,     0,     1,     0,     2,
     5            2,     1,     1,     1,     0,     0,     1,     1,
     6            1,     1,     2,     2,     1,     1,     0,     3,
     7            0,     1,     1,     0,     0,     0,     4,     2,
     8            2,     2,     1,     3,     1,     1,     1,     0,
     9            1,     3,     4,     0,     0,     0,     2,     0,
     a            0,     0,     0,     0,     0,     0,     3,     1,
     b            0,     0,     0,     0,     0,     0,     0,     0,
     c            0,     0,     3,     3,     1,     1,     3,     3,
     d            0,     0,     3,     0,     1,     1,     0,     1,
     e            0,     1,
     u            4,     4,     4,     4,     4,     4,     4,     4,
     u            4,     4/

      data xd/    0,     0,     0,     0,     0,     0,     0,     0,
     1            0,     0,     0,     0,     0,     0,     0,     0,
     2            0,     0,     0,     0,     0,     0,     0,     0,
     3            0,     0,     0,     0,     0,     0,     0,     0,
     4            0,     0,     0,     0,     0,     0,     0,     0,
     5            0,     0,     0,     0,     0,     1,     0,     0,
     6            0,     0,     0,     0,     0,     0,     0,     0,
     7            0,     0,     0,     0,     0,     0,     0,     0,
     8            0,     0,     0,     0,     0,     0,     0,     0,
     9            0,     0,     0,     0,     0,     0,     0,     0,
     a            0,     0,     0,     0,     0,     0,     0,     0,
     b            0,     0,     0,     0,     0,     0,     0,     0,
     c            0,     0,     0,     0,     0,     0,     0,     1,
     d            0,     0,     0,     0,     0,     0,     0,     0,
     e            0,     0,
     u            0,     0,     0,     0,     0,     0,     0,     0,
     u            0,     0/

      data ss/.25d0,.75d0,.25d0,.75d0/,tt/.75d0,.75d0,.25d0,.25d0/
      data tr/.true./, fa /.false./

c-----[--.----+----.----+----.-----------------------------------------]

c     Assign storage for deformed position

      setvar = palloc( 53, 'CT   ', 3*numnp, 2)
      button = 'x'

c     Initialize on first call

      if(.not.pfl) then

c       Set pointers

        npid = id31         ! ID
        npix = np(33)       ! IX
        npuu = np(40)       ! U
        npxx = np(43)       ! X
        nper = np(57)       ! NDER
        npnp = np(58)       ! NDNP
        npev = np(77)       ! EVEC
        nprn = np(89)       ! NREN
        npty = np(190)      ! NDTYP

c       Set initial parameters

        plmesh  =  fa
        fwin    =  fa
        fwoff   =  tr
        psoff   =  tr
        cplxpl  =  tr
        pflag   =  tr
        bordfl  =  tr
        clchk   =  fa
        clip    =  tr
        fopn    =  fa
        hdcpy   =  fa
        hdlogo  =  fa
        psfram  =  fa
        blk     =  fa
        ppr     =  tr
        pfl     =  tr
        hide    =  fa
        labl    =  tr
        outf    =  tr
        psmrk   =  fa
        psmmx   =  fa
        rangfl  =  fa
        symf    =  fa
        edgfl   =  fa
        pltact  =  fa
        timf    =  tr
        torsfl  =  fa
        fact    = 1.0d0
        macd    = 0.0d0
        maci    = 0
        evtyp   = 3
        iclear  = 0
        ienter  = 0
        iexit   = 0
        intvl   = 1
        kr      = 0
        nsizt   = 1
        plix    = npix
        nxd     = nen1
        nxn     = nen
        nne     = numel
        ncapt   = 0
        nsym    = 0
        nreg1   = 0
        nreg2   = mxreg

c       Determine 3-d sizes for surface plots

        wmin(1) = 0.0d0
        wmin(2) = 0.0d0
        wmax(1) = 1.4d0
        wmax(2) = 1.0d0

        setvar = palloc(128,'APLOT',numel  ,1)
        call setclp(hr(npxx),ndm,numnp)
        call plfacn(mr(npix),mr(np(128)),nen,
     &              numel,nface,mr(np(32)),nie)
        setvar = palloc( 54,'FCIX ',7*max(nface,1),1)
        nfaco   = max(nface,1)
        call plfacx(mr(npix),mr(np(128)),mr(np(54)),
     &              nen,numel,mr(np(32)),nie)

        setvar = palloc( 61,'OUTL ',max(nface,numnp)+1,1)
        setvar = palloc( 62,'SYMM ',8*max(nface,numel),1)
        call psetip(mr(np(62)),  max(nface,numel))

        setvar = palloc( 66,'VISN ',numnp, 1)

c       Initialize dof list and symmetry table for no refections

        call pzeroi(iquad,8)
        do i = 1,3
          isymm(i,1) = 0
          pdf(i)     = i
        end do ! i
        nfac(1)    = max(nface,numel)
        call ppermu(isymm)
        nfac(1)    = numel

c       Set initial conditions for no perspective

        kpers   = 0
        call pzero (xsyc,3)
        call pzero (eold,3)
        call pzero (vold,3)
        vold(ndm) = 1.0d0
        ifrm    = 0
        iln(1)  = 0
        iln(2)  = 1
        ilno(1) = 0
        ilno(2) = 1
        nzm1 = 0
        nzm2 = 0
        cs   = 0.0d0
        es   = 1.0d0
        call frame(hr(npxx),ndm,numnp,1)

c       Check for user library names

        do i = 1,10
          k4 = ncomd-10+i
          if(.not.pcomp(upltc(i),wd(k4),4)) then
            wd(k4) = upltc(i)
            ed(k4) = 0
          endif
        end do ! i

c     Check for elements currently active

      else
        if(pltact) then
          pltact = fa
          call setclp(hr(npxx),ndm,numnp)
          call plfacn(mr(npix),mr(np(128)),nen,
     &                numel,nface,mr(np(32)),nie)
          if(nface.gt.nfaco) then
            setvar = palloc( 54,'FCIX ',7*max(nface,1),1)
          endif
          call pzeroi(mr(np(54)),7*nface)
          call plfacx(mr(npix),mr(np(128)),mr(np(54)),
     &                nen,numel,mr(np(32)),nie)
          if(max(nface,numel).gt.max(nfaco,numel)) then
            setvar = palloc( 62,'SYMM ',8*max(nface,numel),1)
          endif
          nfaco = max(nfaco,1)
          call pzeroi(mr(np(62)),8*max(nface,numel))
          call psetip(mr(np(62)),  max(nface,numel))
          if(hide) then
            call pdefm(hr(npxx),hr(npuu),cs,ndm,ndf,numnp,
     &                 hr(np(53)))
            call pppcol(1,0)
            dumv = ct(2)
            call phide(dumv,nxd,nxn,nne,nface,iln)
          endif
        endif

      endif

      ijump = 0

c     Start plot or clear screen

60    if(iclear.eq.0 .and. fpltfl) then
        call plopen
        call plclos
        iclear = 1
        fpltfl = .false.
      endif

      if(ior.lt.0.and.pcomp(lci,'    ',4)) then
        ijump  = 1
        if(rank.eq.0) then
          if(pfeap_on) write(*,'(/)')
          if (ppr) then
            write(*,2002)
            write(xxx,2003) ijump
            call pprint(xxx)
          else
            call pprint(' P>')
          endif
        endif
        setvar = tinput(tx,2,td,4)
        lct    = tx(1)(1:4)
        macd   = td(3)
        maci   = nint(td(3))

        if(pcomp(lct,'end ',4) .or. pcomp(lct,'e   ',4)) then
          setvar = palloc( 53, 'CT   ', 0, 2)
          return
        endif

        if(ior.lt.0.and.pcomp(lct,'help',4)) then
          call phelp(tx(2),wd,ed,list,'PLOT')
          go to 200
        endif

      else
        lct   = lci
        tx(1) = '    '
        tx(2) = '    '
      endif

c     Set initial values for parameters

      ipb    = 0
      dtext  = 0.0d0
      do l = 1,list
        if(pcomp(lct,wd(l),4)) go to 190
      end do ! l
      if(ior.gt.0) write(iow,2001) lct
      if(ior.lt.0) write(*,2001) lct
      go to 200

c     Set 'ct' parameters for interactive mode

190   if(ijump.ne.0) then
        if(xd(l).eq.0) then                    ! Case for no extra word
          setvar = vinput(tx(2),15,ct(1),1)
          ct(2)  = td(1)
          ct(3)  = td(2)
          zcol   = td(3)
        else                                   ! Case for extra word
          ct(1)  = td(1)
          ct(2)  = td(2)
          ct(3)  = td(3)
          zcol   = td(4)
        endif
      else                                     ! Batch mode
        zcol   = 0.0d0
      endif

c     Transfer to requested plot operation

      go to( 1, 1, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,

c            f  w  f  i  c  l  s  c  o  l  m  s  n  b  e  z  c  f  t  s
c            r  i  a  s  a  i  y  o  u  o  e  t  o  o  l  o  o  i  e  i
c            a  p  c  o  r  n  m  n  t  a  s  r  d  u  e  o  l  l  x  z
c            m  e  t  m  t  e  m  t  l  d  h  e  e  n  m  m  o  l  t  e

     &      21,22,23,24,25,26,10,28,29,30,30,32,33, 8, 8,36,37,38,39,40,

c            c  e  b  s  a  p  d  s  h  p  n  d  u  v  a  p  r  e  m  b
c            v  i  o  c  x  e  i  h  i  r  o  e  n  e  c  o  e  i  a  a
c            a  g  r  a  i  r  s  o  d  i  p  f  d  l  c  s  a  g  t  c
c            r  v  d  l  x  s  p  w  e  n  r  o  e  o  e  t  c  e  e  k

     &      41,42,43,44,45,46,47,48,49,50,51,52,53,54,12,56,57,58,59,12,

c            c  t  m  r  p  c  p  p  p  q  r  i  e  d  e  p  p  p  c  p
c            l  i  a  e  i  a  b  f  n  u  e  m  y  o  s  r  r  a  l  s
c            i  t  r  f  c  p  o  o  o  a  a  a  e  f  t  o  a  i  e  t
c            p  i  k  r  k  t  u  r  d  d  l  g  s  s  r  f  x  r  a  r

     &      61,12,63,64,65,66,67,68,69,70,70,72,73,74,75,76,77,78,79,80,

c            d  s  m  p  d  s  p  p  p  l  n  s  p  e  x  w  l  t  b  r
c            p  p  a  r  e  c  d  e  r  a  o  n  s  x  y  i  o  i  p  a
c            l  l  n  o  f  r  i  l  o  b  l  o  n  n  p  n  g  m  l  n
c            o  o  u  m  a  e  s  m  j  e  a  d  o  o  l  d  o  e  o  g

     &      81,82,83,84,85,85,87,37,89,90,12, 8, 8, 8,12,12,22,98,99,

c            n  r  c  s  f  n  u  j  r  n  f  c  v  a  s  p  e  i  s
c            o  e  y  p  u  o  p  i  e  o  l  w  w  w  w  w  w  n  w
c            r  c  l  h  l  f  l  n  g  r  u  i  i  i  i  i  i  t  e
c            a  t  i  e  l  u  o  t  i  m  x  r  r  r  r  r  r  e  e

     &      100,101, 13,103, 8, 61, 61, 10,108,109,110,111, 12, 12,

c             t   e   r   e  n   v   a   r   s   a   d   j   h   s
c             o   d   e   l  d   p   p   f   t   c   e   p   i   t
c             r   e   n   p  a   l   l   o   r   t   a   e   s   r
c             s   f   u   l  t   o   o   r   i   i   c   g   t   a

     &      114,

c             t
c             r
c             i
c             a

     &      800,800,800,800,800,800,800,800,800,800),l

c             u   u   u   u   u   u   u   u   u   u
c             p   p   p   p   p   p   p   p   p   p
c             l   l   l   l   l   l   l   l   l   l
c             1   2   3   4   5   6   7   8   9   0
      go to 200

c     New frame
c     [fram,ifrm] - ifrm = 0 for whole screen
c                   ifrm = 1 for upper-left
c                   ifrm = 2 for upper-right
c                   ifrm = 3 for lower-left
c                   ifrm = 4 for lower-right
c                   ifrm = 5 for legend Box

1     ifrm = nint(ct(1))
      ifrm = min(5,ifrm)
1010  if(ifrm.ge.1 .and. ifrm.le.4) then
         scale = 0.5d0*scaleg*fact
         s0(1) = ss(ifrm)
         s0(2) = tt(ifrm)
      elseif(ifrm.eq.0) then
         scale = scaleg*fact
         s0(1) = 0.5d0
         s0(2) = 0.5d0
      endif

c     [wipe,ifrm,noclean] - ifrm as above - clean screen

      if(mod(intvc,intvl).ne.0) go to 200
      if(l.eq.2) then
        xp(1) = real(min(-ct(1),ct(2)))
        if(ifrm.eq.0) iclear = 0
        call plopen
        call pppcol(-1,0)
        if(ifrm.eq.1) call ppbox(.010d0,.485d0,.476d0,.475d0,1)
        if(ifrm.eq.2) call ppbox(.485d0,.485d0,.475d0,.475d0,1)
        if(ifrm.eq.3) call ppbox(.010d0,.010d0,.476d0,.476d0,1)
        if(ifrm.eq.4) call ppbox(.485d0,.010d0,.475d0,.476d0,1)
        if(ifrm.eq.5) call ppbox(0.98d0,0.10d0,0.30d0,0.70d0,1)
        ifrm = mod(ifrm,5)
      endif
      go to 200

c     [fact,value] - Multiply plot area by 'value'

3     if(ct(1).eq.0.0d0) then
        fact = 1.0
      else
        fact = ct(1)
      endif
      go to 1010

c     [cent]er,s0-1,s0-2  - center graphics in window

4     if(abs(ct(1))+abs(ct(2)).ne.0.0d0) then
        s0(1) = ct(1)
        s0(2) = ct(2)
      else
        s0(1) = ss(ifrm)
        s0(2) = tt(ifrm)
      endif
      go to 200

c     [cart] - cartesian view: (l=5)

5     plix  = npix      ! IX
      nxd   = nen1
      nxn   = nen
      nne   = numel
      nfac(1) = nne
      hide  = fa
      edgfl = fa
      kpers = 0
      call psetip(mr(np(62)),nne)
      call frame(hr(npxx),ndm,numnp,1)
      if(ior.lt.0) then
        write(*,*) ' -> Plot in Cartesian view'
      endif
      go to 200

c     [line,value,width] - 'value' is line type
c                          'width' is line width (device dependent)

6     iln(1) = max(0,min(nint(ct(1)),7))
      iln(2) = max(0,min(nint(ct(2)),5))
      if(ior.lt.0) then
        if(iln(1).eq.0) then
          write(*,2021) iln(2)
        elseif(iln(1).eq.1) then
          write(*,2022) iln(2)
        elseif(iln(1).eq.2) then
          write(*,2023) iln(2)
        elseif(iln(1).eq.3) then
          write(*,2024) iln(2)
        elseif(iln(1).eq.4) then
          write(*,2025) iln(2)
        elseif(iln(1).eq.5) then
          write(*,2026) iln(2)
        elseif(iln(1).eq.6) then
          write(*,2027) iln(2)
        elseif(iln(1).eq.7) then
          write(*,2028) iln(2)
        endif
      else
        if(iln(1).eq.0) then
          write(iow,2021) iln(2)
        elseif(iln(1).eq.1) then
          write(iow,2022) iln(2)
        elseif(iln(1).eq.2) then
          write(iow,2023) iln(2)
        elseif(iln(1).eq.3) then
          write(iow,2024) iln(2)
        elseif(iln(1).eq.4) then
          write(iow,2025) iln(2)
        elseif(iln(1).eq.5) then
          write(iow,2026) iln(2)
        elseif(iln(1).eq.6) then
          write(iow,2027) iln(2)
        elseif(iln(1).eq.7) then
          write(iow,2028) iln(2)
        endif
      endif
      call plopen
      call plline(iln)
      go to 200

c     [symm]; symm,x1,x2,x3  Set symmetry conditions for plots
c             xsyc(1),xsyc(2),xsyc(3)
c      N.B.   xsyc(i) =  coordinate defining reflection point

7     do i = 1,ndm
        isymm(i,1) = 0
        if(ct(i).ne.0.0d0) isymm(i,1) = 1
      end do ! i

720   if(defalt) then
        call pzero(xsyc,3)
      else
        if(ior.lt.0) then
          call pprint(' Input: xsyc(i),i=1,ndm ->')
        endif
        setvar = pinput(xsyc,min(ndm,3))
        if(setvar) go to 720
      endif

      nfac(1) = nne
      call ppermu(isymm)
      symf    = tr

c     Do not rescale if zcol is greater than zero!

      if(zcol.le.0.0d0) then
        if(kpers.ne.0) then
          call frame(hr(npxx),ndm,numnp,-1)
        else
          call frame(hr(npxx),ndm,numnp,1)
        endif
      endif
      if(hide) then
        ct(1) = -2.0d0
        ct(2) =  0.0d0
        ct(3) =  0.0d0
        go to 29
      endif
      go to 200

c     [cont,k1,k2,k3] - k1 = component number (1-ndf)
c                       k2 = # lines (fill if <,= 0)
c                       k3 = 0: superpose mesh
c                       k3 > 0: no mesh
c     [cwir,k1,k2,k3] - Wire model contours

8     if(mod(intvc,intvl).ne.0) go to 200
      if(macd.ne.0.0d0) cs = macd
      if(l.lt.92 .and. kpers.eq.1 .and. .not.hide) then
        if(ior.lt.0) then
          write(*,4000)
        else
          write(ilg,4000)
          write(iow,4000)
          call plstop()
        endif
        go to 200
      endif
      k2  = nint(ct(2))
      ipb = nint(ct(3))
      ni  = abs(k2)
      ni  = max(1,ni)
      i   = max(1,abs(int(ct(1))))
      i   = min(i,ndf)
      if(k2.le.0) then
        ni = -i
      endif
      call pdefm(hr(npxx),hr(npuu),cs,ndm,ndf,numnp,hr(np(53)))
      setvar = palloc(111,'TEMP1',nneq,2)

c     Contours of dependent variables

      if(l.eq.8 .or. l.eq.92) then
        call protv(mr(npty),hr(npuu+kr*nneq),ndf,numnp,hr(np(111)))
        if(iview.eq.0) then
          pl(1) = np(111)
          icp   = ndf
        elseif(iview.eq.1) then
          setvar = palloc(112,'TEMP2',numnp,2)
          call plcylc(hr(npxx),hr(np(111)),i,hr(np(112)))
          pl(1) = np(112)
          icp   = 1
          i     = 1
        else
          setvar = palloc(112,'TEMP2',numnp,2)
          write(*,*) ' SPHERICAL NOT IMPLEMENTED'
        endif
        call rprint(mr(plix),nxn,nxd,nne, hr(pl(1)+i-1),icp,k4)
        if(l.eq.8) then
          call pltcon(hr(np(53)),mr(np(32)),mr(plix),mr(np(62)),
     &                hr(pl(1)),nie,3,icp,nxd,nxn,i,ni,i,2,labl)
        elseif(l.eq.92) then
          call pltwir(hr(np(53)),mr(np(32)),mr(plix),mr(np(62)),
     &                hr(pl(1)),nie,3,icp,nxd,nxn,i,ni,2,labl)
        endif
        if(iview.ge.1) then
          setvar = palloc(112,'TEMP2',0,2)
        endif

c     [velo,k1,k2,k3]
c     [vwir,k1,k2,k3] - Wire model velocity contours

      elseif(l.eq.34 .or. l.eq.93) then
        if(np(42).eq.0) then
          if(ior.lt.0) then
            write(  *,4003)
          else
            write(iow,4003)
            call plstop()
          endif
          go to 200
        endif
        call protv(mr(npty),hr(np(42)),ndf,numnp,hr(np(111)))
        if(iview.eq.0) then
          pl(1) = np(111)
          icp   = ndf
        elseif(iview.eq.1) then
          setvar = palloc(112,'TEMP2',numnp,2)
          call plcylc(hr(npxx),hr(np(111)),i,hr(np(112)))
          pl(1) = np(112)
          icp   = 1
          i     = 1
        else
          setvar = palloc(112,'TEMP2',numnp,2)
          write(*,*) ' SPHERICAL NOT IMPLEMENTED'
        endif
        call rprint(mr(plix),nxn,nxd,nne, hr(pl(1)+i-1),icp,k4)
        if(l.eq.34) then
          call pltcon(hr(np(53)),mr(np(32)),mr(plix),mr(np(62)),
     &                hr(pl(1)),nie,3,icp,nxd,nxn,i,ni,i,3,labl)
        elseif(l.eq.93) then
          call pltwir(hr(np(53)),mr(np(32)),mr(plix),mr(np(62)),
     &                hr(pl(1)),nie,3,icp,nxd,nxn,i,ni,3,labl)
        endif
        if(iview.ge.1) then
          setvar = palloc(112,'TEMP2',0,2)
        endif

c     [acce,k1,k2,k3]
c     [awir,k1,k2,k3] - Wire model velocity contours

      elseif(l.eq.35 .or. l.eq.94) then
        if(np(42).eq.0) then
          if(ior.lt.0) then
            write(  *,4003)
          else
            write(iow,4003)
            call plstop()
          endif
          go to 200
        endif
        call protv(mr(npty),hr(np(42)+nneq),ndf,numnp,hr(np(111)))
        if(iview.eq.0) then
          pl(1) = np(111)
          icp   = ndf
        elseif(iview.eq.1) then
          setvar = palloc(112,'TEMP2',numnp,2)
          call plcylc(hr(npxx),hr(np(111)),i,hr(np(112)))
          pl(1) = np(112)
          icp   = 1
          i     = 1
        else
          setvar = palloc(112,'TEMP2',numnp,2)
          write(*,*) ' SPHERICAL NOT IMPLEMENTED'
        endif
        call rprint(mr(plix),nxn,nxd,nne, hr(pl(1)+i-1),icp,k4)
        if(l.eq.35) then
          call pltcon(hr(np(53)),mr(np(32)),mr(plix),mr(np(62)),
     &                hr(pl(1)),nie,3,icp,nxd,nxn,i,ni,i,4,labl)
        elseif(l.eq.94) then
          call pltwir(hr(np(53)),mr(np(32)),mr(plix),mr(np(62)),
     &                hr(pl(1)),nie,3,icp,nxd,nxn,i,ni,4,labl)
        endif
        if(iview.ge.1) then
          setvar = palloc(112,'TEMP2',0,2)
        endif

c     [ndata,xxxx,#] --- xxxx = disp, stre, pstr
c                           # = filenumber

      elseif(l.eq.104) then

        call pndata(tx,ct,nxd,nxn,nne,labl)

      endif
      setvar = palloc(111,'TEMP1',0,2)
      if(timf) call pltime()
      go to 200

c     [outl,k1,k2] - Outline of 2-d meshes
c              k1  - 3-d outline
c              k2 < 0 alters line color

9     if(mod(intvc,intvl).ne.0) go to 200
      call pdefm(hr(npxx),hr(npuu),cs,ndm,ndf,numnp,hr(np(53)))
      call plopen
      call pppcol(5,0)
      if(ndm.eq.3 .and. edgfl) then
        k2 = nint(ct(2))
        call pout3d(mr(np(54)),mr(np(62)),mr(np(64)),hr(np(53)),
     &              hr(np(44)),hr(np(65)),mr(np(63)),k2,ct(3),
     &              nen,3,nface,numnp)
      else
        call pline(hr(np(53)),mr(np(32)),mr(plix),mr(np(61)),
     &             mr(np(62)),numnp,nne,3,nxd,nxn,nie,ct(2),fa)
      endif
      if(cs.ne.0.0d0 .and. timf) call pltime()
      go to 200

c     [load,k1,k2] - Plot load vector: k1 > 0 tip on node,
c                                      k2 > 0 scale factor
c     [rfor,k1,k2] - Plot radial follower force vector

10    if(mod(intvc,intvl).ne.0) go to 200
      if(l.eq.10) then
        rffl  = np(45).ne.0    ! Convert to Cartesian form for 'angl'es.
        pl(1) = np(27)
        if(ior.lt.0) write(*,'(a)') '  Plot nodal forces'
      elseif(l.eq.107) then
        rffl  = .false.
        pl(1) = np(244)
        if(ior.lt.0) write(*,'(a)') '  Plot radial forces'
      endif
      call plopen
      call pppcol(3,0)
      call pdefm(hr(npxx),hr(npuu),cs,ndm,ndf,numnp,hr(np(53)))
      setvar = palloc(111,'TEMP1',numnp,1)
      call pnumna(mr(np(33)),nen1,nen,numel,mr(np(111)))
      k1 = nint(ct(1))
      if(l.eq.10 .or. l.eq.107) then
        if(rlnew.eq.0.0d0) then
          dumv = 1.0d0
        else
          dumv = rlnew*prop
        endif
        setvar = palloc(112,'TEMP2',nneq,2)
        call pzero(hr(np(112)),nneq)
        call ploads(hr(npuu),hr(np(112)),dumv,tr,fa)
        call setfor(hr(pl( 1)),hr(np(28)),mr(np(31)+nneq),
     &              dumv,nneq, hr(np(112)))
c       Convert reactions to cartesian form if necessary
        if(rffl) then
          call pcartre(hr(np(112)),hr(np(45)),ndf,numnp)
        endif
        call pltfor(hr(np(53)),hr(np(112)),mr(npid),
     &              mr(np(111)),3,ndf,numnp,k1,ct(2),1)
        setvar = palloc(112,'TEMP2',0,2)

c     [disp,k1,k2] - Plot displacement vector: k1 > 0 tip on node
c                                              k2 > 0 scale factor

      elseif(l.eq.27) then
        if(ior.lt.0) write(*,'(a)') '  Plot nodal values'
        call pltfor(hr(np(53)),hr(npuu),mr(npid),
     &              mr(np(111)),3,ndf,numnp,k1,ct(2),2)
      endif
      setvar = palloc(111,'TEMP1',    0,1)
      if(timf) call pltime()
      go to 200

c     [mesh,k1] - plot mesh of current material number
c                 k1 < 0 alters mesh color


11    if(mod(intvc,intvl).ne.0) go to 200
      pl(1) = npuu + kr*nneq
      call pdefm(hr(npxx),hr(pl(1)),cs,ndm,ndf,numnp,hr(np(53)))
      call plopen
      call pppcol(4,0)   ! Blue
      call pppcol(6,0)   ! Cyan
      if(cs.gt.0.0d0) call pppcol(1,0)
      call pline(hr(np(53)),mr(np(32)),mr(plix),mr(np(61)),
     &           mr(np(62)),numnp,nne,3,nxd,nxn,nie,ct(1),tr)
      if(cs.ne.0.0d0 .and. timf) call pltime()
      if(button.eq.'l') then
        call plclos
        go to 53
      end if
      go to 200

c     Stress contours
c     [stre,k1,k2,k3] - k1 = component number (1-ndf) - Projected stress
c     [stra,k1,k2,k3] - k1 = component number (1-ndf) - Projected strain
c     [flux,k1,k2,k3] -      component thermal flux 'k1' in elements
c     [hist,k1]       - k1 = history component number

c     [pstr,k1,k2,k3] - k1 = component number (1-7)   - Principal values
c     [estr,k1,k2,k3] - k1 = component number (1-ndf)
c                       k2 = # lines (fill if <,= 0)
c                       k3 = 0: superpose mesh
c                       k3 > 0: no mesh
c     [swir,k1,k2,k3] - Wire model stress contours
c     [pwir,k1,k2,k3] - Wire model principal stress
c     [ewir,k1,k2,k3] - Wire model element stress

12    if(mod(intvc,intvl).ne.0) go to 200
      k1 = nint(abs(ct(1)))
      if(l.eq.91) then
        k1    = min(13,max(1,k1) + 12)
        ncapt = 1
        write(caption,'(a,i2)') '  F L U X   ',k1-6
      elseif(l.eq.112) then
        k1    = min(hplmax,max(1,k1))
        ncapt = 1
        write(caption,'(a,i2)') '  HISTORY   ',k1
      endif
      k2  = nint(ct(2))
      ipb = nint(ct(3))
      if((l.lt.93 .or. l.eq.112 .or. l.eq.113) .and.
     &             kpers.eq.1 .and. .not.hide) then
        if(ior.lt.0) then
          write(*,4000)
        else
          write(iow,4000)
          call plstop()
        endif
        go to 200
      endif

c     Allocate storage for element variable projections

      if(plfl) then
        setvar = palloc( 57,'NDER',numnp*8    ,2)
        setvar = palloc( 58,'NDNP',numnp*npstr,2)
        setvar = palloc( 60,'NDNS',max(nen*npstr,nst*nst),2)
        setvar = palloc(207,'NSCR',numel      ,2)
        nper   = np(57)
        npnp   = np(58)
        plfl   =.false.
        if(histpltfl) then
          setvar = palloc(304,'HSELM',nen*hplmax  ,2)
          setvar = palloc(305,'HDNP ',numnp*hplmax,2)
        endif
      endif
      ner    = nper
      nph    = npnp

      pl(1) = nph + numnp
      pl(2) = ner + numnp
      k5    = 1

c     Compute and project stresses

      if(.not.fl(11)) then
        call pjstrs(trifl)
        praxf = .false.
      endif
      call pdefm(hr(npxx),hr(npuu),cs,ndm,ndf,numnp,hr(np(53)))
      if(l.eq.12 .or. l.eq.91 .or. l.eq.95) then
        k1    = min(npstr-1,max(k1,1))
        pl(2) = pl(1) + numnp*(k1 - 1)
      elseif(l.eq.112) then
        pl(2) = np(305) + numnp*(k1-1)
      elseif(l.eq.113) then
        k1    = min(npstr-1,max(k1,1))
        pl(2) = pl(1) + numnp*(k1 + 5)
      endif
      if(l.eq.12.or.l.eq.91.or.l.eq.95.or.l.eq.112.or.l.eq.113) then
        k5 = 1
        if(l.eq.113) k5 = 10
        if(iview.eq.0) then
        elseif(iview.eq.1) then
          setvar = palloc(112,'TEMP2',numnp,2)
          pl(2)  = np(112)
          call plcyls(hr(npxx),hr(pl(1)),k1,hr(pl(2)))
        else
          setvar = palloc(112,'TEMP2',numnp,2)
          write(*,*) ' SPHERICAL NOT IMPLEMENTED'
        endif
      elseif(l.eq.55) then
        k1    = min(npstr-1,max(k1,1))
        pl(2) = pl(1) + numnp*(k1 - 1)
        k5    = 1
      elseif(l.eq.60 .or. l.eq.96) then
        k1    = min(7,max(k1,1))
        pl(2) = pl(2) + numnp*(k1-1)
        k5    = 5
      else
        go to 61
      endif
      call rprint(mr(plix),nxn,nxd,nne, hr(pl(2)),1,k2)
      if(k2.le.0) then
        k2 = -k1
      endif
      k3 = 1

c     No projection case

      if(l.eq.55) then
        if(.not.hide) then
          setvar = palloc(111,'TEMP1',numel*nen,2)
          call pltelc(hr(np(53)),mr(np(32)),mr(plix),mr(np(62)),
     &                hr(np(111)),nie,3,nxd,k2,k1,k5,labl)
          setvar = palloc(111,'TEMP1',0,2)
          fl(11) = fa
        endif
      else
        if(l.eq.95 .or. l.eq.96) then
          call pltwir(hr(np(53)),mr(np(32)),mr(plix),mr(np(62)),
     &              hr(pl(2)),nie,3,k3,nxd,nxn,k3,k2,k5,labl)
        else
          call pltcon(hr(np(53)),mr(np(32)),mr(plix),mr(np(62)),
     &              hr(pl(2)),nie,3,k3,nxd,nxn,k3,k2,k1,k5,labl)
        endif
        fl(11) = tr
      endif
      if((l.eq.12.or.l.eq.91.or.l.eq.95.or.l.eq.112.or.l.eq.113) .and.
     &    iview.gt.0) then
        setvar = palloc(112,'TEMP2',0,2)
      endif
      if(timf) call pltime()
      go to 200

c     [node,k1,k2] - Plot nodes k1 to k2 & numbers
c                               k1 = 0 show all nodes & numbers
c                               k1 < 0 show all nodes, no numbers

13    if(mod(intvc,intvl).ne.0) go to 200
      call pdefm(hr(npxx),hr(npuu),cs,ndm,ndf,numnp,hr(np(53)))
      setvar = palloc(111,'TEMP1',numnp,1)
      call pnumna(mr(plix),nxd,nxn,nne,mr(np(111)))
      call plopen
      call pppcol(5,1)
      k4 = nint(ct(1))
      if(k4.lt.0) then
        k1 = 0
        k2 = 1
        k3 = numnp
      elseif(k4.eq.0) then
        k1 = 1
        k2 = 1
        k3 = numnp
      else
        k1 = 1
        k2 = max( 1,min(numnp,int(ct(1))))
        k3 = max(k2,min(numnp,int(ct(2))))
        call pppcol(7,0)
      endif
      if(l.eq.102) k1 = -k1
      call pltnod(hr(np(53)),mr(np(111)),  3,numnp,k1,k2,k3)
      setvar = palloc(111,'TEMP1',    0,1)
      go to 200

c     [boun,k1] - Plot/label boundary restraints; k1 > 0

14    if(mod(intvc,intvl).ne.0) go to 200
      call pdefm(hr(npxx),hr(npuu),cs,ndm,ndf,numnp,hr(np(53)))
      setvar = palloc(111,'TEMP1',numnp  ,1)
      setvar = palloc(112,'TEMP2',numnp*3,1)
      call pnumna(mr(np(33)),nen1,nen,numel,mr(np(111)))
      call pnumba(mr(np(33)),mr(np(32)),mr(np(112)))
      call plopen
      call pppcol(2,0)
      k1 = nint(ct(1))
      call pltbou(mr(np(31)),hr(np(53)),hr(np(45)),mr(np(111)),
     &            mr(np(112)),3,ndf,numnp,k1)
      setvar = palloc(112,'TEMP2',    0,1)
      setvar = palloc(111,'TEMP1',    0,1)
      go to 200

c     [elem],n1,n2 - Label elements with numbers n1 to n2 (default=all)

15    if(mod(intvc,intvl).ne.0) go to 200
      call pdefm(hr(npxx),hr(npuu),cs,ndm,ndf,numnp,hr(np(53)))
      call plopen
      call pppcol(3,1)
      if(int(ct(1)).eq.0) then
        k2 = 1
        k3 = nne
      else
        k2 = max( 1,min(nne,int(ct(1))))
        k3 = max(k2,min(nne,int(ct(2))))
      endif
      call pltelm(hr(np(53)),mr(np(32)),mr(plix),scale,nie,3,
     &            nxd,k2,k3)
      go to 200

c     [zoom,k1,k2] Set window for zoom; k1 = 1st node; k2 = 2nd node

16    if(mod(intvc,intvl).ne.0) go to 200
      nzm1 = nint(ct(1))
      nzm2 = nint(ct(2))
      fwin = .false.
161   call pdefm(hr(npxx),hr(npuu),cs,ndm,ndf,numnp,hr(np(53)))

c     Construct perspective projection if necessary

      if(kpers.ne.0) then
        call frame(hr(np(53)),  3,numnp,-1)
        call perspj(hr(np(53)),mr(npty),numnp,errv)
        if(errv) kpers = 0
      else
        call frame(hr(np(53)),  3,numnp,1)
      endif
      go to 200

c     [colo,k1,k2] - set color
c                    k1 <  0  greyscale postscript
c                    k1 >= 0  color postscript
c                    k2 =  0  standard color order
c                    k2 != 0  reversed color order

17    icol = nint(ct(1))
      if(ct(1).lt.0.0d0) then
        if(ior.lt.0) write(*,*) 'PostScript - grayscale'
        pscolr = .false.
        psrevs = .false.
        icol   = 1
      else
        pscolr = .true.
        psrevs = ct(2).ne.0.0d0
        if(ior.lt.0. .and. psrevs) then
          write(*,*) 'PostScript - color order reversed'
        elseif(ior.lt.0) then
          write(*,*) 'PostScript - normal color order'
        endif
      endif
      call plopen
      call pppcol(icol,0)
      go to 200

c     [fill,k1],<k2> - fill current material in color 'k1'
c                    - no 'Time =' displayed if 'k2' non-zero.

18    if(mod(intvc,intvl).ne.0) go to 200
      call pdefm(hr(npxx),hr(npuu),cs,ndm,ndf,numnp,hr(np(53)))
      k1 = nint(ct(1))
      if(ct(3).eq.0.0d0) then
        k3 =  1
      else
        k3 = -1
      endif
      call plopen
      call plot2d(mr(np(32)),mr(plix),mr(np(62)),hr(np(53)),
     &            hr(np(44)),nie,3,nxn,nxd,numnp,nne,k1,k3)
      if(cs.ne.0.0d0 .and. ct(2).eq.0.0d0 .and. timf) call pltime()
      go to 200

c     [text,k1,x,y] - Put text at 'x,y'  (0<x<1.4; 0<y<1) in color 'k1'

19    if(mod(intvc,intvl).ne.0) go to 200
      k1 = max(1,int(ct(1)))
      call plttxt(k1,ct(2))
      go to 200

c     [size,k1] - Select text size 'k1'

20    nsizt = nint(ct(1))
      nsizt = min(3,max(1,nsizt))
      go to 200

c     [cvar,k1,k2,k3] - Plot variable 'k1' from history vector 'k2'
c                       for pair 'k3'
c                       k1 = 0 -> def = 1
c                       k2 = 0 -> def = 2, plot variable from CH2
c                       k3 = 0 -> plot variable for all pairs

21    if(mod(intvc,intvl).ne.0) go to 200
      call contact (308)
      if(timf) call pltime()

      go to 200

c     [eigv,k1,k2,k3] - Plot eigenvector for value 'k1' & material 'k2'
c                       k2 = 0 for all materials. k3-dof for contouring

22    if(mod(intvc,intvl).ne.0) go to 200

c     Check that eigenvectors and values exist

      if(np(77).eq.0) then
        if(ior.lt.0) then
          write(*,4004)
        else
          write(iow,4004)
          call plstop()
        endif
        go to 200
      endif

c     Set mesh description to plot mode

      if(l.lt.97 .and. kpers.eq.1 .and. .not.hide) then
        if(ior.lt.0) then
          write(*,4000)
        else
          write(iow,4000)
          call plstop()
        endif
        go to 200
      endif

c     Set label

      if(maci.gt.0) then
        evtyp = max(1,min(3,maci))
      endif

c     Compute current deformed state of coordinates

      call plopen
      k1 = max(1,min(mq,int(abs(ct(1)))))
      call pdefm(hr(npxx),hr(npuu),cs,ndm,ndf,numnp,hr(np(53)))

c     Add unit scaled eigenvector to deformed coordinates

      setvar = palloc(111,'TEMP1',nneq,2)
      call pmovec(mr(npid),hr(npev+(k1-1)*vneq), hr(np(111)),
     &            nneq)  ! Expand eigenvector into hr(np(111))
      if(np(167).ne.0) then
        call ruplnk(hr(np(111)),hr(np(111)),hr(npxx),mr(np(100)),
     &              mr(np(167)),1,ndm,ndf,numnp,.false.)
      endif
      call scalev(hr(np(111)),pdf,3,ndf,numnp,.true.)
      call pdefm(hr(np(53)),hr(np(111)),es,3,ndf,numnp,hr(np(53)))

      if(kpers.gt.0 .and. hide) then
        dumv = 0.0d0
        call phide(dumv,nxd,nxn,nne,nface,iln)
      endif

      i = nint(abs(ct(3)))
      if(i.eq.0) then
        if(kpers.gt.0 .and. hide) then
          k1 = nint(ct(1))
          k4 = 1
          call plot2d(mr(np(32)),mr(plix),mr(np(62)),hr(np(53)),
     &                hr(np(44)),nie,3,nxn,nxd,numnp,nne,k1,k4)
        else
          call pppcol(1,0)
          call pline(hr(np(53)),mr(np(32)),mr(plix),mr(np(61)),
     &               mr(np(62)),numnp,nne,3,nxd,nxn,nie,ct(2),tr)
        endif

      else
        i  = min(i,ndf)
        k2 = nint(ct(2))
        ni = abs(k2)
        ipb = 0
        if(ct(3).lt.0.0d0) ipb = 1
        if(k2.le.0) then
          ni = -i
        endif
        if(iview.eq.0) then
          pl(1) = np(111)
          icp   = ndf
        elseif(iview.eq.1) then
          setvar = palloc(112,'TEMP2',numnp,2)
          call plcylc(hr(npxx),hr(np(111)),i,hr(np(112)))
          pl(1) = np(112)
          icp   = 1
          i     = 1
        else
          setvar = palloc(112,'TEMP2',numnp,2)
          write(*,*) ' SPHERICAL NOT IMPLEMENTED'
        endif
        call rprint(mr(plix),nxn,nxd,nne, hr(pl(1)+i-1),icp,k4)
        if(l.eq.97) then
          call pltwir(hr(np(53)),mr(np(32)),mr(plix),mr(np(62)),
     &                hr(pl(1)),nie,3,icp,nxd,nxn,i,ni,2,labl)
        else
          call pltcon(hr(np(53)),mr(np(32)),mr(plix),mr(np(62)),
     &                hr(pl(1)),nie,3,icp,nxd,nxn,i,ni,i,2,labl)
        endif
        if(iview.ge.1) then
          setvar = palloc(112,'TEMP2',0,2)
        endif
      endif
      setvar = palloc(111,'TEMP1',0,2)

c     Place value on screen

      if(labl) then
        call pleigt(hr(np(76)+k1-1))
      endif

      go to 200

c     [bord,color] Toggle border

23    if(pcomp(tx(2),'off',3) .or. ct(1).eq.1.0d0) then
        bordfl = .false.
      else
        bordfl = .true.
        icol = nint(ct(1))
        call plopen
        call plbord(icol)
      endif
      go to 200

c     [scal,cs,c2,c3] - Rescale for deformed window

24    nzm1 = 0
      nzm2 = 0
      jj   = 2
      cs   = ct(1)
      if(cs.eq.0.0d0) cs = 1.0
      call pdefm(hr(npxx),hr(npuu),cs,ndm,ndf,numnp,hr(np(53)))
      if(ct(2).eq.0.0d0) call pwind(hr(npxx),hr(np(53)),ndm,3,numnp)
      if(ct(3).eq.0.0d0) jj = 1

c     Construct perspective projection if necessary

      if(kpers.ne.0) then
        call frame(hr(np(53)),  3,numnp,-jj)
        call perspj(hr(np(53)),mr(npty),numnp,errv)
        if(errv) kpers = 0
      else
        call frame(hr(np(53)),  3,numnp,jj)
      endif
      go to 200

c     [axis,x,y] - Plot axes at screen coords x,y

25    if(mod(intvc,intvl).ne.0) go to 200
      xm = 0.1*max(xmax(1)-xmin(1),xmax(2)-xmin(2),xmax(3)-xmin(3))
      call plopen
      call pppcol(1,0)
      if(kpers.eq.0) then
        k1 = ndm
      else
        k1 = 3
      endif
      call pltaxs(ct,k1,xm)
      go to 200

c     [pers] Input parameters for perspective projection
c     [pers,on,zview] - zview = ??

26    if(mod(intvc,intvl).ne.0) go to 200
      if(ior.lt.0) then
        write(*,*) ' -> Plot in perspective view'
      endif
      call plopen
      hide    = fa
      kpers   = 1
      odefalt = defalt
      defalt  = defalt .and. ct(1).eq.0.0d0
      zview   = ct(2)
      call pppcol(5,0)
      call pdefm(hr(npxx),hr(npuu),cs,ndm,ndf,numnp,hr(np(53)))
      call maxcor(hr(np(53)),  3,numnp)
      call perspe(pflag)
      call frame(hr(np(53)),  3,numnp,-1)
      call perspj(hr(np(53)),mr(npty),numnp,errv)
      defalt = odefalt
      if(errv) then
        kpers = 0
        go to 200
      endif
      plix    = npix
      nxd     = nen1
      nxn     = nen
      nne     = numel
      nfac(1) = nne
      call psetip(mr(np(62)),nne)
      if(.not.pflag) go to 29
      go to 200

c     [show] - Display current plot parameters

28    if(ior.lt.0) then
        write(*,2000) fplt,npart,ifrm,iln(1),cs,fact
        if(kpers.eq.0) then
          write(*,2004)
        else
          write(*,2005)
        endif
      endif
      go to 200

c     [hide,nopl,back,color] Plot visible mesh

29    if(l.eq.29 .and. kpers.eq.0) then
        write(iow,4001)
        if(ior.lt.0) write(*,4001)
        go to 200
      endif
      call pdefm(hr(npxx),hr(npuu),cs,ndm,ndf,numnp,hr(np(53)))
      hide  = tr
      pflag = tr
      if(l.eq.29) then
        call plopen
        call pppcol(4,0)
        if(ct(1).lt. 0.0d0) call pppcol(1,0)
        if(ct(1).lt.-1.0d0) call pppcol(0,0)
        k1 = nint(ct(3))
        if(k1.gt.0) call pppcol(k1,0)
      else
        call pppcol(-1,0)
      endif

c     Set plot mesh quantities

      dumv = ct(2)
      call phide(dumv,nxd,nxn,nne,nface,iln)

c     Plot current mesh

      if(plmesh) then
        plmesh = fa
        go to 11
      endif
      go to 200

c     [prin] or [nopr]

30    ppr = l.eq.30
      go to 200

c     [defo,scale,resize,escale]
c           scale  = multiplier on displacements (default=1)
c           resize = non-zero: do not rescale plot region
c           escale = multiplier on eigenvectors  (default=1)

32    cs = ct(1)
      es = ct(3)
      if(max(cs,es).eq.0.0d0) cs = 1.d0

      if(ior.lt.0) write(*,2006) cs,es
      write(iow,2006) cs,es
      if(ct(2).eq.0.0d0) go to 161
      if(ior.lt.0) write(*,2008)

      go to 200

c     [unde,,resize,escale]; do not rescale if resize is non-zero

33    cs = 0.0d0
      if(ior.lt.0) write(*,2007) es

      es = ct(3)

      if(ct(2).eq.0.0d0) go to 161
      if(ior.lt.0) write(*,2008)

      go to 200

c     [post] or [post,1] (0=portrait mode, 1=landscape mode)
c            PostScript -  open or close file

36    if(mod(intvc,intvl).ne.0) go to 200
      if (.not. hdcpy) then

c       Open PostScript file, set hdcpy=.true.

        hdcpy = tr
        psfram = (ct(1) .ne. 0.0d0)
        hdlogo = (ct(2) .eq. 0.0d0)
        if(ct(3).ne.0.0d0) then
          ct(3) = min(1.0d0,max(0.1d0,ct(3)))
        else
          ct(3) = 1.0d0
        endif
        call fppsop(ct(3))
      else

c       Close PostScript file, set hdcpy=.false.

        call fpplcl()
        hdcpy = fa
      endif
      goto 200

c     [reac,k1,k2]  Plot nodal reactions; k1=non-zero gives tip at node
c     [jint,k1,k2]  Plot fracture forces as reactions;
c                                         k2 > 0 scale factor

37    if(mod(intvc,intvl).ne.0) go to 200
      call plopen
      call pppcol(5,0)
      setvar = palloc(111,'TEMP1',numnp,1)
      call pnumna(mr(np(33)),nen1,nen,numel,mr(np(111)))
      k1 = nint(ct(1))
      setvar = palloc(112,'TEMP2',max(neq,nneq),2)
      call pzero(hr(np(112)),max(neq,nneq))
      rlold = rlnew
      call ploa1(ttim,dt)

c     Reaction option

      if(l.eq.37) then
        k2   =  6
        rffl = np(45).ne.0    ! Convert to Cartesian form for 'angl'es.
        pltmfl = .true.

c     Fracture force option

      else
        k2   = 16
        rffl = .false.
      endif
c     Compute element forms
      call formfe(npuu,np(112),np(112),np(112),
     &            fa,tr,fa,tr,k2,1,numel,1)
      pltmfl = .false.
c     Convert reactions to cartesian form if necessary
      if(rffl) then
        call pcartre(hr(np(112)),hr(np(45)),ndf,numnp)
      endif
      call pdefm(hr(npxx),hr(npuu),cs,ndm,ndf,numnp,hr(np(53)))
      call pltfor(hr(np(53)),hr(np(112)),mr(npid),
     &            mr(np(111)),3,ndf,numnp,k1,ct(2),-1)
      setvar = palloc(112,'TEMP2',    0,2)
      setvar = palloc(111,'TEMP1',    0,1)
      if(timf) call pltime()
      rlnew = rlold
      go to 200

c     [eige,k1,c2] - k1 = vector number - plot last element eigenvector

38    if(mod(intvc,intvl).ne.0) go to 200
      if(np(75).ne.0) then

c       Add mode shape to last element coordinates

        if(es.eq.0.0d0) then
          dumv = 1.d0
        else
          dumv = es
        endif
        call pleige(ct,mr(np(32)),mr(npix),hr(npxx),hr(np(53)),
     &              hr(np(75)+nst*nst),hr(np(75)),dumv)
      else
        if(ior.lt.0) then
          write(*,*) ' *ERROR* PPLOTF: Use command EIGE first.'
        endif
      endif
      goto 200

c     [mate,#] - Set material number of plot (reset project flag if new)

39    k1    = maplt
      k2    = msplt
      maplt = max(0,min(nint(ct(1)),nummat))
      msplt = max(0,min(nint(ct(2)),nummat))
      if(maplt.ne.k1 .or. msplt.ne.k2) then
        fl(11) = .false.
        rfl    = .false.
      endif
      if(prnt) then
        if(ior.lt.0) then
          write(  *,2009) maplt,msplt,.not.fl(11)
        else
          write(iow,2009) maplt,msplt,.not.fl(11)
        endif
      endif

c     Set element material set number to all processed

      if(maplt.eq.0) then

c       Set normal material numbers to zero value

        do i = 0,numel-1
          mr(np(128)+i) = 0
        end do ! i

c     Set element material set numbers so only 'maplt' is processed

      else

        k1 = 0
        do i = nen1-1,nen1*numel-1,nen1
          if(mr(np(33)+i).eq.maplt) then
            mr(np(128)+k1) =  0
          else
            mr(np(128)+k1) = -1
          endif
          k1 = k1 + 1
        end do ! i

      endif

c     Set new face list

391   call plfacn(mr(np(33)),mr(np(128)),nen,
     &            numel,nface,mr(np(32)),nie)
      if(nface.gt.nfaco) then
        setvar = palloc( 54,'FCIX ',7*max(nface,1),1)
        nfaco  = nface
      endif
      call pzeroi(mr(np(54)),7*nface)
      call plfacx(mr(np(33)),mr(np(128)),mr(np(54)),
     &            nen,numel,mr(np(32)),nie)
      if(max(nface,numel).gt.max(nfaco,numel)) then
        setvar = palloc( 62,'SYMM ',8*max(nface,numel),1)
      endif
      call psetip(mr(np(62)),  max(nface,numel))

      if(hide) go to 29

      go to 200

c     [back] - Set background for post script - black = white too

40    if(ct(1).ne.0.0d0) then
        blk = .false.
      else
        blk = .true.
      end if
      go to 200

c     [clip,x-dir,x-min,x-max] - clip in x-dir between x-min anad x-max
c     [clip]                   - restore full plot

41    k1 = nint(ct(1))
      k1 = min(k1,ndm)

c     Set a clip plane min/max

      if(k1.gt.0) then
        cmin(k1) = min(ct(2),ct(3))
        cmax(k1) = max(ct(2),ct(3))
        if(ior.lt.0) write(  *,2010) k1,cmin(k1),cmax(k1)
        if(ior.gt.0) write(iow,2010) k1,cmin(k1),cmax(k1)
      else

c       Reset all clip plane min/max

        call setclp(hr(npxx),ndm,numnp)
        if(ior.lt.0) write(  *,*) 'Reset all clip planes to full mesh'
        if(ior.gt.0) write(iow,*) 'Reset all clip planes to full mesh'
      endif
      if( (k1.gt.0) .and. (cmin(k1).eq.cmax(k1)) ) then
        if(ior.lt.0) write(  *,*) 'Zero thickness clip defined'
        if(ior.gt.0) write(iow,*) 'Zero thickness clip defined'
      else
        call plfacn(mr(npix),mr(np(128)),nen,
     &              numel,nface,mr(np(32)),nie)
        if(nface.gt.nfaco) then
          setvar = palloc( 54,'FCIX ',7*max(nface,1),1)
          nfaco  = nface
        endif
        call pzeroi(mr(np(54)),7*nface)
        call plfacx(mr(npix),mr(np(128)),mr(np(54)),
     &              nen,numel,mr(np(32)),nie)
        if(max(nface,numel).gt.max(nfaco,numel)) then
          setvar = palloc( 62,'SYMM ',8*max(nface,numel),1)
        endif
        call psetip(mr(np(62)),  max(nface,numel))

c       Initialize symmetry table for no refections

        isymm(1,1) = 0
        isymm(2,1) = 0
        isymm(3,1) = 0
        nfac(1)    = max(nface,numel)
        call ppermu(isymm)
        nfac(1)    = numel
        ct(1)      = 0.0d0
        ct(2)      = 0.0d0
        if(ndm.eq.3) go to 29
      endif
      go to 200

c     [titl]e - Place problem title on plot

42    if(mod(intvc,intvl).ne.0) go to 200
      k1 = max(1,nint(ct(1)))
      call pltitl(k1)
      go to 200

c     [mark,n1] - Place marks on contour plots to locate max/min
c                 any nonzero 'n1' turns off marks.

43    psmrk = ct(1).eq.0.0d0
      psmmx = ct(2).ne.0.0d0
      if(prnt) then
        if(ior.lt.0) then
          write(*,2030) psmrk,psmmx
        endif
        write(iow,2030) psmrk,psmmx
      endif
      go to 200

c     [refr]esh - Refresh an X-window

44    call plopen
      go to 200

c     [pick] - Pick a point

45    call plopen
      call pdefm(hr(npxx),hr(npuu),cs,ndm,ndf,numnp,hr(np(53)))
      call pick()
      go to 200

c     [capt]ion - Define Caption for next contour plot

46    if(mod(intvc,intvl).ne.0) go to 200
      if(ior.lt.0) then
        ncapt   = 1
        caption = tx(2)
      else
        write(iow,3003)
      endif
      go to 200

c     [pbou,idir] - Pick boundary conditions with a mouse

47    call plopen
      call pdefm(hr(npxx),hr(npuu),cs,ndm,ndf,numnp,hr(np(53)))
      k1 = nint(ct(1))
      k1 = max(1,min(k1,ndf))
      call psboun(mr(np(31)+nneq),mr(npty),hr(np(53)),k1,ndf,3,numnp)

      do i = 0,nneq-1
        mr(npid+i) = mr(np(31)+i+nneq)
      end do ! i

      call profil(mr(np(21)),mr(np(34)),mr(npid),mr(npix),1,
     &            .false.)
      if(lkflg) then
        call plink(mr(npid),mr(npty),ndf,numnp,neq,.false.)
      endif
      if(leflg) then
        setvar = palloc(120,'TEMP0',numnp,1)
        call pelink(mr(npid),hr(npxx),ndm,ndf,numnp,neq,.false.,
     &              mr(np(257)),mr(np(120)))
        setvar = palloc(120,'TEMP0',    0,1)
      endif
      call profil(mr(np(21)),mr(np(34)),mr(npid),mr(npix),2,
     &            .false.)

      go to 200

c     [pfor,idir,val] - Pick nodal forces with a mouse

48    call plopen
      call pdefm(hr(npxx),hr(npuu),cs,ndm,ndf,numnp,hr(np(53)))
      k1 = nint(ct(1))
      k1 = max(1,min(k1,ndf))
      if(ct(2).ne.0.0d0) then
        call psforc(hr(np(27)),mr(npty),hr(np(53)),k1,ct(2),
     &              ndf,3,numnp)
      else
        write(*,*) ' ** WARNING: Specify forced condition is zero.'
      endif

      go to 200

c     [pnod] - Pick a node

49    call plopen
      call pdefm(hr(npxx),hr(npuu),cs,ndm,ndf,numnp,hr(np(53)))
      call psnode(hr(np(53)),  3,numnp, k1)
      go to 200

c     [quad] - Set active plot quadrant

50    if(symf) then
        call pltqud(iquad,ndm)
        call ppermu(isymm)
      else
        write(iow,3002)
        if(ior.lt.0) write(*,3002)
      endif
      go to 200

c     [real] - Set active plot to real part

51    write(iow,2013)
      if(ior.lt.0) write(*,2013)
      kr     = 0
      cplxpl = tr
      go to 200

c     [imag] - Set active plot to imaginary part

52    write(iow,2014)
      if(ior.lt.0) write(*,2014)
      kr     = 3
      cplxpl = fa
      go to 200

c     [eyes] - Pick 3-d view using eyes

53    if(kpers.eq.0) then
        write(iow,4001)
        if(ior.lt.0) write(*,4001)
        go to 200
      else
        write(iow,2011)
        if(ior.lt.0) write(*,2011)
        call eyepik(button)
        if(button.ne.'l') then
          pflag  = tr
          plmesh = fa
          go to 200
        endif
        iclear = 0
        pflag  = fa
        plmesh = tr
        go to 26
      endif

c     [dofs,k1,k2,k3] - Reset degree-of-freedoms for deformed plots

54    do i = 1,3
        pdf(i) = nint(ct(i))
      end do ! i

      write(iow,2012) pdf
      if(ior.lt.0) write(*,2012) pdf

      go to 200

c     [prof]ile,<ct.eq.0> - Plot layout of upper profile
c     [prof]ile,<ct.ne.0> - Plot layout of total profile

56    call plopen
      call pltprf(mr(np(21)),neq,ct(1).ne.0.0d0)
      go to 200

c     [prax,k1,k2,k3] - Compute and plot principal stress directions
c     k1 -- principal stress vector to be plotted, 0 = all
c     k2 -- < 0 negative vectors only, >0 positive vectors only
c           = 0 all.
c     k3 -- >= 0 optional color shift
c           <  0 scaling factor for vectors

57    if(mod(intvc,intvl).ne.0) go to 200
      k1 = nint(ct(1))
      k2 = nint(ct(2))
      if(ct(3) .lt. 0.0) then
       dumv = 0.3d0*abs(ct(3))
       k3 = 1
      else
       dumv = 0.3d0
       k3 = max(1,min(int(ct(3)),7))
      endif
      call plopen

c     Allocate storage for element variable projections

      if(plfl) then
        setvar = palloc( 57,'NDER',numnp*8    ,2)
        setvar = palloc( 58,'NDNP',numnp*npstr,2)
        setvar = palloc( 60,'NDNS',max(nen*npstr,nst*nst),2)
        setvar = palloc(207,'NSCR',numel      ,2)
        nper   = np(57)
        npnp   = np(58)
        plfl   =.false.
        if(histpltfl) then
          setvar = palloc(304,'HSELM',nen*hplmax  ,2)
          setvar = palloc(305,'HDNP ',numnp*hplmax,2)
        endif
      endif
      ner    = nper
      nph    = npnp

c     Compute and project stresses

      if(.not.fl(11)) then
        call pjstrs(.false.)
        praxf = .false.
      endif

c     Get storage and calculate principal directions

      call pdefm(hr(npxx),hr(npuu),cs,ndm,ndf,numnp,hr(np(53)))

      setvar = palloc( 59,'NPAX', numnp*(ndm-1)*6 , 2)
      if(.not.praxf) then
        pl(5) = nph + numnp
        call prax(mr(npty),hr(np(59)),hr(pl(5)),numnp,ndm)
        praxf = .true.
      endif

      xm = dumv*max(xmax(1)-xmin(1),xmax(2)-xmin(2),xmax(3)-xmin(3))
     &         /sqrt(dble(numel))

      call praxp(hr(np(59)),hr(np(53)),numnp,ndm,xm,k1,k2,k3)

      go to 200

c     [pair,k1,k2,k3] - Plot for pairs from 'k1' to 'k2' surface 'k3'
c                       k3 = 0 -> plot slave and master surf
c                       k3 = 1 -> plot slave surf
c                       k3 = 2 -> plot master surf

58    if(mod(intvc,intvl).ne.0) go to 200
      call pdefm(hr(npxx),hr(npuu),cs,ndm,ndf,numnp,hr(np(53)))
      call maxcor(hr(np(53)),  3,numnp)
      call contact (305)
      if(timf) call pltime()

      go to 200

c     [clea]r screen area

59    if(mod(intvc,intvl).ne.0) go to 200
      call plopen
      call pppcol(-1,0)
      call ppbox(.01d0,.02d0,.95d0,.95d0,1)
      go to 200

c     [dplo,k1,k2,k3] set up node list for displ. plots along a line
c     [vplo,k1,k2,k3] set up node list for veloc. plots along a line
c     [aplo,k1,k2,k3] set up node list for accel. plots along a line
c     [splo,k1,k2,k3] set up node list for stress plots along a line

61    k1 = nint(abs(ct(1)))
      k2 = nint(abs(ct(2)))
      k2 = min(k2,12)
      k3 = nint(ct(3))
      k4 = max(1,min(k1,npstr-1))
      call pdefm(hr(npxx),hr(npuu),cs,ndm,ndf,numnp,hr(np(53)))
      call plopen
      setvar = palloc(115,'TEMP5',nneq  ,1)
      setvar = palloc(116,'TEMP6',nneq*2,2)

c     [dplo] plot curve of displacements

      if(l.eq.61) then
       k4 = min(k4,ndf)
       setvar = palloc(111,'TEMP1',nneq,2)
       call protv(mr(npty),hr(npuu+kr*nneq),ndf,numnp,hr(np(111)))
       call plotn(hr(np(111)+k4-1),hr(np(53)),mr(npix),k4,k2,k3,
     &            ndf,nen,nen1,numnp,numel,
     &            mr(np(115)),hr(np(116)),'d')
       setvar = palloc(111,'TEMP1',0,2)

c     [vplo] plot curve of velocities

      elseif(l.eq.105) then
       if(np(42).eq.0) then
         if(ior.lt.0) then
           write(  *,4002) 'velocity'
         else
           write(iow,4002) 'velocity'
         endif
       else
         k4 = min(k4,ndf)
         setvar = palloc(111,'TEMP1',nneq,2)
         call protv(mr(npty),hr(np(42)),ndf,numnp,hr(np(111)))
         call plotn(hr(np(111)+k4-1),hr(np(53)),mr(npix),k4,k2,k3,
     &              ndf,nen,nen1,numnp,numel,
     &              mr(np(115)),hr(np(116)),'v')
         setvar = palloc(111,'TEMP1',0,2)
       endif

c     [aplo] plot curve of accelerations

      elseif(l.eq.106) then
       if(np(42).eq.0) then
         if(ior.lt.0) then
           write(  *,4002) 'acceleration'
         else
           write(iow,4002) 'acceleration'
         endif
       else
         k4 = min(k4,ndf)
         setvar = palloc(111,'TEMP1',nneq,2)
         call protv(mr(npty),hr(np(42)+nneq),ndf,numnp,hr(np(111)))
         call plotn(hr(np(111)+k4-1),hr(np(53)),mr(npix),k4,k2,k3,
     &              ndf,nen,nen1,numnp,numel,
     &              mr(np(115)),hr(np(116)),'a')
         setvar = palloc(111,'TEMP1',0,2)
       endif

c     [splo] plot curve of stresses

      elseif(l.eq.62) then
       k1    = max(k1,1)
       pl(1) = nph + k4*numnp
       call plotn(hr(pl(1)),hr(np(53)),mr(npix),k1,k2,k3,1,
     &            nen,nen1,numnp,numel,
     &            mr(np(115)),hr(np(116)),'s')
      endif
      setvar = palloc(116,'TEMP6',   0,2)
      setvar = palloc(115,'TEMP5',   0,1)
      go to 200

c     [manu]al,hlplev - Set manual help display level

63    hlplev = max(-1,min(3,int(ct(1))))
      go to 200

c     [prom]pt,<on,off> - Set graphics prompt level

64    if(pcomp(tx(2),'off',3) .or. ct(1).eq.1.d0) then
        prompt = .false.
        if(ior.lt.0) write(*,*) ' Prompt: OFF'
      else
        prompt = .true.
        if(ior.lt.0) write(*,*) ' Prompt: ON'
      endif
      go to 200

c     [defa]ult,<on,off> - Set graphics defaults

65    if(pcomp(tx(2),'off',3) .or. ct(1).eq.1.d0) then
        defalt = .false.
        if(ior.lt.0) write(*,*) ' Defaults: OFF'
      else
        defalt = .true.
        if(ior.lt.0) write(*,*) ' Defaults: ON'
      endif
      go to 200

c     [scre]en,<on,off> or <0,1>  Turn on/off screen for PostScript only

66    if(pcomp(tx(2),'off',3) .or. ct(1).eq.1.d0) then
        screfl = .false.
        if(ior.lt.0) write(*,*) ' Screen: OFF during plots'
      else
        screfl = .true.
        if(ior.lt.0) write(*,*) ' Screen: ON during plots'
      endif
      go to 200

c     [pdis,idir,val] - Pick nodal displacements with a mouse

67    call plopen
      call pdefm(hr(npxx),hr(npuu),cs,ndm,ndf,numnp,hr(np(53)))
      k1 = nint(ct(1))
      k1 = max(1,min(k1,ndf))
      call psforc(hr(np(27)+nneq),mr(npty),hr(np(53)),k1,ct(2),
     &            ndf,3,numnp)
      go to 200

c     [pele,val(i),i=1,3] - Plot element variables

68    if(mod(intvc,intvl).ne.0) go to 200
      call plopen
      call pppcol(1,0)
      elplt(1) = ct(1)
      elplt(2) = ct(2)
      elplt(3) = ct(3)
      call formfe(npuu,np(26),np(26),np(26),
     &            fa,fa,fa,fa,20,1,numel,1)
      go to 200

c     [proj] - Set flag to reproject stresses

69    fl(11) = .false.
      if(ior.lt.0) write(*,*) '  Project element values for next plot'
      go to 200

c     [labe]l     - Set flag to add  plot label/scales
c     [labe]l off - Set flag to omit plot label/scales
c     [nola]bel   - Set flag to omit plot label/scales

70    labl = l.eq.70
      if(pcomp(tx(2),'off',3) .or. ct(1).eq.1.d0) then
        labl = .false.
      endif
      if(ior.lt.0) then
        if(labl) then
          write(*,*) '  Add label/scale to contour plots'
        else
          write(*,*) '  No label/scale for contour plots'
        endif
      endif
      go to 200

c     [snod,k1,k2] - Plot supernodes nodes k1 to k2 & numbers
c                               k1 < 0 show all nodes
c                               k1 = 0 show all nodes & numbers

72    if(mod(intvc,intvl).ne.0) go to 200
      k4 = nint(ct(1))
      if(k4.lt.0) then
        k1 = 0
        k2 = 1
        k3 = numsn
      elseif(k4.eq.0) then
        k1 = 1
        k2 = 1
        k3 = numsn
      else
        k1 = 1
        k2 = max( 1,min(numsn,int(ct(1))))
        k3 = max(k2,min(numsn,int(ct(2))))
      endif
      call plopen
      call pppcol(5,0)
      setvar = palloc(111,'TEMP1',numsn,1)
      do k4 = 0,numsn-1
        mr(np(111)+k4) = 1
      end do ! k4
      call pltnod(hr(np(161)),mr(np(111)),3,numsn,k1,k2,k3)
      setvar = palloc(111,'TEMP1',    0,1)
      go to 200

c     [psno]de - Pick a super node

73    if(numsn.gt.0) then
        call plopen
        call psnode(hr(np(161)),  3,numsn, k1)
        if(k1.gt.0) then
          pl(1) = np(161) + 3*k1 - 4
          write(*,2015) (hr(pl(1)+k3),k3=1,ndm)
          call pprint('  Specify new coordinates (y or n): ')
          setvar = tinput(tx,1,td,0)
          if(pcomp(tx(1),'y',1)) then
            call pprint('  Specify new coordinates:')
            setvar = pinput(td,ndm)
            do k3 = 1,ndm
              hr(pl(1)+k3) = td(k3)
            end do ! k3
          endif

c         Perform delayed mesh generation steps

          if(numbd.gt.0) then
            call pblendm(15,20,ndm,nen1,prt,prth,.false.,.false.)
          endif

        endif

      else
        if(ior.lt.0) write(*,*) '   No super nodes in mesh'
      endif

      go to 200

c     [exno] - k1 = 0: Plot exterior nodes numbers
c              k1 < 0: Plot exterior nodes, no numbers

74    if(mod(intvc,intvl).ne.0) go to 200
      call pdefm(hr(npxx),hr(npuu),cs,ndm,ndf,numnp,hr(np(53)))
      k4 = nint(ct(1))
      if(k4.lt.0) then
        k1 = 0
      else
        k1 = 1
      endif
      call plopen
      call pppcol(5,0)
      pl(1) = npty - 1
      do k2 = 1,numnp
        if(mr(pl(1)+k2).ge.0. and. mr(np(78)-1+k2).gt.0) then
          call pltnod(hr(np(53)),mr(np(78)), 3,numnp,k1,k2,k2)
        end if
      end do ! k2
      go to 200

c     [xypl] - ct: x-y plots

75    write(*,*) ' Option not active yet'
c75    call plotxy(ct)
      go to 200

c     [wind]ow,<number> - Set window number

76    iwindow = max(1,min(nint(ct(1)),4))
      call plopen
      go to 200

c     [logo] - Place program logo on window

77    if(mod(intvc,intvl).ne.0) go to 200
      call plopen
      call pfeap(0.983d0,0.017d0,0.250d0,3,3)
      call pfeap(0.982d0,0.018d0,0.250d0,2,1)
      call pfeap(0.981d0,0.019d0,0.250d0,2,1)
      call pfeap(0.980d0,0.020d0,0.250d0,2,1)
      go to 200

c     [time],<on,off> - Set graphics time flag

78    if(pcomp(tx(2),'off',3) .or. ct(1).eq.1.d0) then
        timf = .false.
        if(ior.lt.0) write(*,*) '  No TIME added to plots'
      else
        timf = .true.
        if(ior.lt.0) write(*,*) '  TIME added to plots'
      endif
      go to 200

c     [bplo] - ct: beam surface plots

79    if(mod(intvc,intvl).ne.0) go to 200
      call plopen
      call pppcol(1,1)
      call bplot(ct)
      go to 200

c     [rang] - ct: range for fill plots

80    if(pcomp(tx(2),'off',3)) then
        rangfl = .false.
        if(ior.lt.0) write(*,2016)
      else
        rangfl = .true.
        rangmn = min(ct(1),ct(2))
        rangmx = max(ct(1),ct(2))
        if(ior.lt.0) write(*,2017) rangmn,rangmx
      endif
      go to 200

c     [nora]nge - off (batch mode)

81    rangfl = .false.
      if(ior.lt.0) write(*,2016)
      go to 200

c     [rect]angular - stress plots

82    iview = 0
      zview0(1) = 0.0d0
      zview0(2) = 0.0d0
      zview0(3) = 0.0d0
      if(ior.lt.0) write(*,2018)
      go to 200

c     [cyli]ndrical - stress plots

83    iview = 1
      zview0(1) = ct(1)
      zview0(2) = ct(2)
      zview0(3) = ct(3)
      if(ior.lt.0) write(*,2019) (i,ct(i),i=1,ndm)
      go to 200

c     [sphe]rical - stress plots

84    iview = 2
      zview0(1) = ct(1)
      zview0(2) = ct(2)
      zview0(3) = ct(3)
      if(ior.lt.0) write(*,2020) (i,ct(i),i=1,ndm)
      go to 200

c     [full] - Set to full screen

85    k1 = l - 85
      call plopen
      call pppcol(1,1)
      call pfullscr(k1)
      iclear = 0
      go to 200

c     [uplo] - ct: user plots

87    call plopen
      call pppcol(1,0)
      call uplot(ct)
      if(timf) call pltime()
      go to 200

c     [regi,n1,n2] - Set plots for regions n1 to n2.
c     [regi]       - Plot all regions

89    nreg1 = min(nint(ct(1)),nint(ct(2)))
      nreg2 = max(nint(ct(1)),nint(ct(2)))
      if(nreg2.eq.0) then
        nreg1 = 0
        nreg2 = mxreg
      elseif(nreg1.eq.0) then
        nreg1 = nreg2
      endif
      if(prnt) then
        if(ior.lt.0) then
          write(*,2029) nreg1,nreg2
        endif
        write(iow,2029) nreg1,nreg2
      endif
      go to 391

c     [norm,,k2] - Plot frame normals; k2 > 0 for scaling

90    if(mod(intvc,intvl).ne.0) go to 200
      call plopen
      call pppcol(3,0)
      call pdefm(hr(npxx),hr(npuu),cs,ndm,ndf,numnp,hr(np(53)))
      setvar = palloc(111,'TEMP1',numnp,1)
      call pnumna(mr(plix),nxd,nxn,nne,mr(np(111)))
      k1 = nint(ct(1))
      call pltfor(hr(np(53)),hr(np(206)),mr(npid),
     &            mr(np(111)),3, 3,numnp,k1,ct(2),2)
c     if(ndm.eq.3) then
c       call formfe(npuu,np(26),np(26),np(26),
c    &              fa,fa,fa,fa,30,1,numel,1)
c     else
c       write(*,*) ' Normal for 3-d elements only'
c     endif
      go to 200

c     [inte]rval,<intvl> - Plot interval for contours

98    intvl = max(1,nint(ct(1)))
      intvc = 0
      if(prnt) then
        write(iow,2031) intvl
        if(ior.lt.0) then
          write(*,2031) intvl
        endif
      endif
      go to 200

c     [swee]p,<angle> - Plot 3-d mesh for plane problem

99    if(mod(intvc,intvl).ne.0) go to 200
      call pdefm(hr(npxx),hr(npuu),cs,ndm,ndf,numnp,hr(np(53)))
      swang = ct(1)
      nsinc = nint(ct(2))
      if(swang.le.0.0d0) swang = 45.0
      if(nsinc.eq.0    ) nsinc = 5
      if(prnt) then
        write(iow,2032) swang,nsinc
        if(ior.lt.0) then
          write(*,2032) swang,nsinc
        endif
      endif
      call plopen
      call psweep(swang,nsinc,nxd,nxn)
      go to 200

c     [tors]ional projection for deformed

100   if(pcomp(tx(2),'off',3) .or. ct(1).eq.1.0d0) then
        torsfl = .false.
        if(ior.lt.0) then
          write(  *,2033)
        else
          write(iow,2033)
        endif
      else
        torsfl = .true.
        if(ior.lt.0) then
          write(  *,2034)
        else
          write(iow,2034)
        endif
      endif
      go to 200

c     [edef]orm eigenvector

101   call evscal(hr(npxx),ct(1), es)
      if(prnt) then
        write(iow,2035) es
        if(ior.lt.0) then
          write(*,2035) es
        endif
      endif
      go to 200

c     [elpl]ot,ne - plot of element ne

103   k1 = max(1,min(numel,nint(ct(1)))) - 1
      k2 = nint(ct(2))
      k3 = max(0,min(nen,nint(ct(3))))
      if(prnt) then
        write(iow,2036) k1+1
        if(ior.lt.0) then
          write(*,2036) k1+1
        endif
      endif
      call plopen
      call pppcol(1,0)
      call plotelm(mr(np(32)), mr(plix+nen1*k1),
     &             hr(npxx), k2.ne.0, k3)
      go to 200

c     [string,<text>,k1  - Put string on plot

108   if(mod(intvc,intvl).ne.0) go to 200
      k1     = max(1,nint(ct(1)))
      if(ior.lt.0) then
        write(*,2037)
      endif
      setvar = tinput(tx,2,td,0)
      call plstring(k1)
      go to 200

c     [acti,k1,k2,k3  - Activate   region k1 to k2 in increments k3

109   call pactivate(tx(1),ct, 1)
      go to 200

c     [deac,k1,k2,k3  - Deactivate region k1 to k2 in increments k3

110   call pactivate(tx(1),ct,-1)
      go to 200

c     [jpeg] - Dump screen to JPEG file

111   call pjpeg()
      go to 200

c     [tria]ds for shells

114   if(mod(intvc,intvl).ne.0) go to 200
      k1 = nint(ct(1))
      k2 = nint(ct(2))
      if(ct(3) .lt. 0.0) then
       dumv = 0.4d0*abs(ct(3))
       k3 = 1
      else
       dumv = 0.4d0
       k3 = max(1,min(int(ct(3)),7))
      endif
      call plopen

c     Set deformed configuration

      call pdefm(hr(npxx),hr(npuu),cs,ndm,ndf,numnp,hr(np(53)))

c     Allocate storage for director vectors

      setvar = palloc( 59,'NPAX', numnp*(ndm-1)*6 , 2)

c     Move director vectors into array for plots

      call pldir(mr(npty),hr(np(59)),hr(np(82)),numnp)

      xm = dumv*max(xmax(1)-xmin(1),xmax(2)-xmin(2),xmax(3)-xmin(3))
     &         /sqrt(dble(numel))

      call pdirp(hr(np(59)),hr(np(53)),numnp,ndm,xm,k3)

c     setvar = palloc( 59,'NPAX', 0 , 2)

      go to 200

c     User plot option library

800   k4 = l + 10 - ncomd
      call plopen
      call pppcol(1,0)
      uct = upltc(k4)
      call upltlib(k4,ct(1))
      go to 200

c     Close plot

200   call plclos
      if(ijump.ne.0) go to 60

      setvar = palloc( 53, 'CT   ', 0, 2)
      return

c     Formats

2000  format(3x,'P L O T    P A R A M E T E R S'//
     &       5x,'Plot file   = ',a36/5x,'Partition   =',i3/
     &       5x,'Frame No.   =',i3/5x,'Line Type   =',i3/
     &       5x,'Disp. Scale =',e12.5/5x,'Fact  Scale =',e12.5)

2001  format(' ** WARNING: no match on a plot,',a4,' request.')

2002  format(' Input PLOT instructions, use "help" for a list,',
     &       ' exit with "end".')
2003  format('     Plot ',i3,'> ')

2004  format(5x,'Cartesian view')

2005  format(5x,'Perspective view')

2006  format('  -> Plots on deformed mesh with scale:',1p,e12.5/
     &       '                     Eigenvector scale:',1p,e12.5)

2007  format('  -> Plots on undeformed mesh:',
     &       ' Eigenvector scale:',1p,e12.5)

2008  format('     No rescale for plot.')

2009  format('     Plot results for material set identifier',i3/
     &       '                      material set number    ',i3/
     &       '     Project = ',l2,'  (0 - for all sets)')

2010  format('  Set Clip plane: direction =',i3,' Min =',1p,e11.4,
     &       ' Max =',1p,e11.4)

2011  format('  Pick 3-d view from location of eyes, exit by pressing',
     &       ' RIGHT button or ESC')

2012  format('  Plot degree-of-freedoms reassigned to:',
     &       ' 1 =',i2,': 2 =',i2,': 3 =',i2/1x)

2013  format('  Plot real part of complex solution')

2014  format('  Plot imaginary part of complex solution')

2015  format('  Coordinates of node are:',1p,3e12.4)

2016  format('  Range: Automatic (off)')

2017  format('  Range: Min =',1p,1e12.4,', Max =',1p,1e12.4)

2018  format('  Plots in rectangular cartesian components')

2019  format('  Plots in cylindrical coordinate components'/
     &       '  Origin at: x0_',i1,' = ',1p,1e12.4:/
     &               (13x,'x0_',i1,' = ',1p,1e12.4:))

2020  format('  Plots in spherical coordinate components'/
     &       '  Origin at: x0_',i1,' = ',1p,1e12.4:/
     &               (13x,'x0_',i1,' = ',1p,1e12.4:))

2021  format('  Line type 0: Solid          ; Width =',i3)
2022  format('  Line type 1: Dotted         ; Width =',i3)
2023  format('  Line type 2: Dash-dot       ; Width =',i3)
2024  format('  Line type 3: Short-dash     ; Width =',i3)
2025  format('  Line type 4: Long-dash      ; Width =',i3)
2026  format('  Line type 5: Dot-dot-dash   ; Width =',i3)
2027  format('  Line type 6: Short+long dash; Width =',i3)
2028  format('  Line type 7: wide dash      ; Width =',i3)

2029  format('     Plot results for regions'/
     &       '        Minimum = ',i4/
     &       '        Maximum = ',i4/)
2030  format(/'     Plot label markers'
     &       /'        Min/Max  = ',l2
     &       /'        Location = ',l2/)

2031  format(/'  Plot interval = ',i5/)

2032  format(/'  Axisymmetric sweep angle (degrees) = ',1p,1e12.4/
     &        '  Number of increments =',i4/)

2033  format(/'  Cartesian view in deformed plots')

2034  format(/'  Axisymmetric view with torsion in deformed plots')

2035  format(/'  Eigenvalue scale factor = ',1p,1e12.4/)

2036  format(/'  Plot of individual element = ',i10/)

2037  format( '  Specify string: ',$)

3002  format(' *WARNING* Specify SYMMetry before using QUADrant'/)

3003  format(' *WARNING* Issue CAPTion in command mode for BATCh'/)

4000  format(' *ERROR* PPLOTF: Must specify HIDE first'/)

4001  format(' *ERROR* PPLOTF: Must specify PERSpective view first'/)

4002  format(' *ERROR* PPLOTF: Cannot plot ',a,' for static problem'/)

4003  format(' *ERROR* PPLOTF: Must specify transient option first'/)

4004  format(' *ERROR* PPLOTF: Must compute eigen-solution first'/)

      end
