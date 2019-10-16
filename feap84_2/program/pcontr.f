c$Id:$
      subroutine pcontr()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add 'gaugm' parameter                            14/03/2007
c       2. Set 'pfeap_glob' to false                        13/04/2007
c       3. Add 'lagr'ange multiplier nodal indicators       15/06/2007
c       4. Increase 'nie' from 11 to 12 to store element    20/07/2007
c          multiplier partion number
c       5. Add end of mesh and end of contact messages      02/08/2007
c       6. Relocate 2009/2010 outputs                       21/08/2007
c       7. Add .true. to faceset call                       04/09/2007
c       8. Add call to pidset after ties and coord inputs   21/09/2007
c       9. Add 'np(257)' array use for ELINk                26/12/2007
c      10. Correct set of umacr number to be 1 to 0         09/01/2008
c      11. Correct setting of g2type material model type    05/03/2008
c      12. Add initialization of 'geqnum'.                  27/03/2008
c          Increase 'nie' from 12 to 13 to store global
c          equation flag.
c      13. Increase nen1 to nen + 9: nen+6 = type elmt      12/11/2008
c      14. Set ncurv = 0 and curfl to false                 13/11/2008
c      15. Remove 'ndf' from call to meshck                 13/12/2008
c      16. Adjust call of acheck to 256                     21/12/2008
c      17. Initialize 'prt' and set false for pnums use     08/01/2009
c      18. Add load table option and set parameters         10/01/2009
c      19. Remove 'tief' set to .false. in 'tie off'        02/03/2009
c          Permits new ties later.
c      20. Set pltmfl to .false. for default calls to ormfe 05/03/2009
c      21. Add read of 'ldnum,ldprp,ldflg'; 'pload1.h'      09/03/2009
c          Set 'sptot','spnum' & 'spflg', write spnum,spflg
c      22. Add new subroutine 'pnewprob' and calls at       08/04/2009
c          statements 100 and binary file inputs
c          Set format 1006 to 3i8 not 3i5
c      23. Add 'ufea'p option                               10/04/2009
c      24. Remove 'geof'eap option (use 'feap')             21/04/2009
c      25. Exit include file when 'stop' commaned found     02/05/2009
c      26. Add call to 'pendprob' for multiple runs.        13/05/2009
c          Add common 'comsav.h' for 'prob_on' and set.
c      27. Check 'initf' to set initial partition           08/06/2009
c      28. Remove 'prt' argument from usetlib               09/07/2009
c      29. Add 'file' input option, needed for loops        15/07/2009
c      30. Add 'fpltfl = .true.' for first call to pplotf   30/08/2009
c      31. Add 'keepfl' options to keep or delete outputs   16/11/2009
c      32. Add 'nreg = -1' to force history initializing    23/11/2009
c      33. Add initialization for hsize(1:2) & hksize(1:2)  17/03/2010
c      34. Add and set 'cplxpr' to true                     09/08/2010
c      35. Set length of file extension to 8                25/01/2011
c          Change 'lnk' to 'link' and 'eln' to 'elin'
c      36. Set 'tfl' true for 'link', 'elin', 'clin'        06/06/2011
c          Add 'setext' for 'link' and 'elin'
c      37. Add call to pdxmsh to compute size of coords     29/06/2011
c      38. Add 'tx' to argument on call to umshlib          26/09/2011
c      39. Check 20 user macro routines                     25/01/2012
c      40. Change usub*15, remove .false. form upltlib      01/05/2012
c      41. Make option for rve's 'fe2feap'                  07/05/2012
c      42. Increase dof count on 'orde','part','lagr'       30/05/2012
c      43. Remove unused common block includes              19/05/2012
c      44. Add call to contact(316); move location to after 25/01/2013
c          contact is active. Activates contact nodes for
c          parallel solutions.
c      45. Save 'maxpart' for checking purposes             19/06/2013
c      46. Add 'cprt' - set to .false. by default           10/11/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Control program for FEAP problem input and solution.

c      Inputs:
c        none

c      Outputs:
c        none
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'bdata.h'
      include   'cdata.h'
      include   'cdat1.h'
      include   'cdat2.h'
      include   'chdata.h'
      include   'codat.h'
      include   'comfil.h'
      include   'compac.h'
      include   'complx.h'
      include   'comsav.h'
      include   'contrl.h'
      include   'conval.h'
      include   'corset.h'
      include   'crotas.h'
      include   'debugs.h'
      include   'eltran.h'
      include   'errchk.h'
      include   'gltran.h'
      include   'hdatam.h'
      include   'hlpdat.h'
      include   'idptr.h'
      include   'iodata.h'
      include   'iofile.h'
      include   'ioincl.h'
      include   'iosave.h'
      include   'mxsiz.h'
      include   'part0.h'
      include   'part1.h'
      include   'part3.h'
      include   'part7.h'
      include   'pdatps.h'
      include   'pfeapb.h'
      include   'pointer.h'
      include   'plflag.h'
      include   'pload1.h'
      include   'prflag.h'
      include   'print.h'
      include   'p_point.h'
      include   'qudshp.h'
      include   'region.h'
      include   'rigid1.h'
      include   'rjoint.h'
      include   'sdata.h'
      include   'trdata.h'
      include   'umac1.h'
      include   'comblk.h'

      logical    setvar,palloc,pinput,tinput,vinput,pcomp,lp_in
      logical    errs,evint,cprt,mchflg,oprt,oprth,newprob,usetfl(10)
      logical    lopen,partfl,endwh
      character  titl*81,dnam*15, fext*8, type*4, usub*15
      character  uset(10)*4, wfac(2)*8, tx(8)*15
      integer    pln,pre, iorsv
      integer    i, iif, j,jj, l1,l2,l3,l4,l5, lp_lun, ndum
      integer    usetno(10)
      real*8     gaptol,td(16)

      save

c     Default names for manipulation sets

      data      uset / 'man1', 'man2', 'man3', 'man4', 'man5',
     &                 'man6', 'man7', 'man8', 'man9', 'man0' /

      data      wfac / 'Material','Region' /

c     Destroy old output and log files if they exist

      inquire(file=fout,exist=initf)
      if(initf) then
        open (unit=iow,file=fout,status='old')
        close(unit=iow,          status='delete')
      endif

c     Set file for log

      call fileset(fout,flog,'O', 'L')

c     Open files for input, output and log

      open(unit=ior,file=finp,status='old')
      open(unit=iow,file=fout,status='new')
      open(unit=ilg,file=flog,status='new')
      call uscreen(1,iow)
      write(ilg,2018)

c     Initial values for include options

      prob_on = .false.
      chflg   = .false.
      cprt    = .true.
      debug   = .false.
      eofile  = .false.
      eralloc = .false.
      everon  = .false.
      evint   = .false.
      fpltfl  = .true.
      hdcpy   = .false.
      incf    = .false.
      intr    = .false.
      intx    = .false.
      keepfl  = .true.
      lagrfl  = .false.
      lmate   = .false.
      lp_in   = .true.
      lread   = .false.
      lsave   = .false.
      mchflg  = .false.
      newprob = .false.
      nocount = .true.
      pltmfl  = .false.
      prt     = .true.
      conprt  = .false.
      ucount  = .false.
      lfile   = ios
      icf     = icl
      isf     = 1
      iif     = 1
      irdef   = ior
      fincld(1) = finp
      irecrd(1) = 0

c     Set program constants

      call pvalues()

c     Set problem counter and initialize clock

      nprob    = 0
      call stime()

c     Initial value for variable storage

      cplxfl = .false.     ! Problem default is real
      cplxpr = .true.      ! Print   default is real

c     Set default to print headers

      prth   = .true.

c     Flags for user manipulation commands

      do j = 1,10
        usetfl(j) = .false.
        usetno(j) = 0
      end do ! j

c     Install user functions

c     Set user mesh command names

      do j = 1,10
        if(j.lt.10) then
          write(usub,'(a3,i1)') 'mes',j
        else
          write(usub,'(a4)') 'mes0'
        endif
        uct = usub(1:4)
        call umshlib(j,tx,prt)
        umshc(j) = uct
      end do ! j

c     Set user solution command names

      do j = 1,20
        if(j.lt.10) then
          write(usub,'(a3,i1)') 'mac',j
        else
          write(usub,'(a2,i2)') 'ma',j
        endif
        uct   = usub(1:4)
        fnamp = ' '
        call umaclib(j,fnamp,td,.false.)
        umacc(j) = uct
      end do ! j

c     Set user material model names

      do j = 1,10
        if(j.lt.10) then
          write(usub,'(a3,i1)') 'mat',j
        else
          write(usub,'(a4)') 'mat0'
        endif
        uct = 'mate'
        call uconst(usub,td,td,td,l1,l2,l3)
      end do ! j

c     Set user plot command names

      do j = 1,10
        if(j.lt.10) then
          write(usub,'(a3,i1)') 'plt',j
        else
          write(usub,'(a4)') 'plt0'
        endif
        uct = usub(1:4)
        call upltlib(j,td)
        upltc(j) = uct
      end do ! j

c     Set user manipulation command names

      do j = 1,10
        uct     = uset(j)
        call usetlib(j)
        uset(j) = uct
      end do ! j

c     Input with interactive interactive statements

1     if(intx) then
        if(cprt) then
          call pprint(
     &        ' Continue with interactive input options for control?')
          call pprint('<y or n> :')
        endif
        errck = tinput(dnam,1,td,0)

c       Read command interactively

        if(pcomp(dnam,'y',1)) then
          call pprint(
     &        ' Specify command (INTEractive, INCLude, STOP, etc.) >')
          read (*,1001,err=900,end=910) yyy
          cprt = .true.

c       Read command from current file and turn off intx flag

        else
          cprt = .false.
          intx = .false.
          ior  =  abs(ior)
          read(ior,1001,err=900,end=910) yyy
        endif

c     Input from current file

      else
        ior = abs(ior)
        read(ior,1001,err=900,end=910) yyy
      endif

c     Compare with command list

      call pstrip(xxx,yyy,1)
      l1   = min(80,len(xxx))
      titl = xxx(1:l1)

c     Start solution of new problem

      if(pcomp(titl(1:4),'feap',4)) then
        go to 100

c     Set count/nocount or input mode to binary

      elseif(pcomp(titl(1:4),'noco',4)) then
        nocount = .false.

      elseif(pcomp(titl(1:4),'coun',4)) then
        nocount = .true.

c     Set keep/nokeep flags for output file retension

      elseif(pcomp(titl(1:4),'keep',4)) then
        keepfl = .true.

      elseif(pcomp(titl(1:4),'noke',4)) then
        keepfl = .false.

c     Set flag for constructing interface mesh

      elseif(pcomp(titl(1:4),'matc',4) .or.
     &       pcomp(titl(1:4),'clus',4)) then

c       Increase storage for local solution arrays

        iif    = 2
        setvar = palloc( 34,'LD   ',max(nen+1,nst*4,21),1)
        setvar = palloc( 35,'P    ',nst*6*ipc          ,2)
        setvar = palloc( 36,'S    ',nst*nst*4*ipc      ,2)
        setvar = palloc( 39,'TL   ',nen*2              ,2)
        setvar = palloc( 41,'UL   ',nst*14*ipc         ,2)
        setvar = palloc( 46,'ANGL ',nen*6              ,2)
        setvar = palloc( 44,'XL   ',max(4,nen)*6       ,2)

        if(pcomp(titl(1:4),'matc',4)) then
          mchflg = .true.
        else
          setvar = palloc(210,'MATCH', 8*numel, 1)
          call faclset(hr(np(25)),mr(np(32)),mr(np(33)),
     &                 mr(np(210)),l2, l4)
          setvar = palloc(210,'MATCH', 8*l2   , 1)

c         Allocate storage for interface history

          setvar = palloc(214,'HINTE',max(1,l4),2)
          setvar = palloc(215,'HINT1',max(1,hnimax),2)
          setvar = palloc(216,'HINT2',max(1,hnimax),2)
          setvar = palloc(217,'HINT3',max(1,hni3max),2)
        endif

c     Set input for a binary mode

      elseif(pcomp(titl(1:4),'bina',4)) then
        newprob = .true.
        bflg    = .true.
        call acheck(titl,yyy,18,80,80)
        fnamr   = ' '
        fnamr   = yyy(19:36)
        fext    = ' '
        fext    = titl(1:4)
        call opnfil(fext,fnamr, 3,ios,lread)
        rewind ios
        read(ios) head,numnp,numel,nummat,ndm,ndf,nen,ndd,nud
        npd = ndd - nud - 1
        nad = 0
        call pnewprob(2)
        go to 1

c     User command sets

      elseif(pcomp(titl(1:4),uset(1),4)) then
        usetno(1) = usetno(1) + 1
        usetfl(1) = .true.
        fext      = 'u1a'
        go to 300
      elseif(pcomp(titl(1:4),uset(2),4)) then
        usetno(2) = usetno(2) + 1
        usetfl(2) = .true.
        fext      = 'u2a'
        go to 300
      elseif(pcomp(titl(1:4),uset(3),4)) then
        usetno(3) = usetno(3) + 1
        usetfl(3) = .true.
        fext      = 'u3a'
        go to 300
      elseif(pcomp(titl(1:4),uset(4),4)) then
        usetno(4) = usetno(4) + 1
        usetfl(4) = .true.
        fext      = 'u4a'
        go to 300
      elseif(pcomp(titl(1:4),uset(5),4)) then
        usetno(5) = usetno(5) + 1
        usetfl(5) = .true.
        fext      = 'u5a'
        go to 300
      elseif(pcomp(titl(1:4),uset(6),4)) then
        usetno(6) = usetno(6) + 1
        usetfl(6) = .true.
        fext      = 'u6a'
        go to 300
      elseif(pcomp(titl(1:4),uset(7),4)) then
        usetno(7) = usetno(7) + 1
        usetfl(7) = .true.
        fext      = 'u7a'
        go to 300
      elseif(pcomp(titl(1:4),uset(8),4)) then
        usetno(8) = usetno(8) + 1
        usetfl(8) = .true.
        fext      = 'u8a'
        go to 300
      elseif(pcomp(titl(1:4),uset(9),4)) then
        usetno(9) = usetno(9) + 1
        usetfl(9) = .true.
        fext      = 'u9a'
        go to 300
      elseif(pcomp(titl(1:4),uset(10),4)) then
        usetno(10) = usetno(10) + 1
        usetfl(10) = .true.
        fext       = 'u0a'
        go to 300

c     User problem selection

      elseif(pcomp(titl(1:5),'ufeap',5)      .or.
     &       pcomp(titl(1:7),'fe2feap',7) ) then

        newprob = .true.
        do i = 1,20
          l2 = 4*i+ 1
          l1 = l2 - 3
          head(i) = titl(l1:l2)
        end do ! i
        call uprob(titl,prt)

c     Print/noprint

      elseif(pcomp(titl(1:4),'prin',4)) then

        prt = .true.

      elseif(pcomp(titl(1:4),'nopr',4)) then

        prt = .false.

c     Perform inputs from an include file

      elseif(pcomp(titl(1:4),'incl',4)) then
        call acheck(titl,yyy,15,80,80)
        read(yyy,1002,err=900,end=900) titl(1:4),dnam
        if(pcomp(dnam,'end',3)) then
          call pincld(dnam)
          if(evint) then
            write(*,2005) fnamr
          endif
          write(iow,2005) fnamr
        else
          fnamr =  dnam
          inquire(file = fnamr, exist = errs)
          if(.not.errs) then
            write(ilg,3003) fnamr(1:15)
            write(iow,3003) fnamr(1:15)
            write(  *,3003) fnamr(1:15)
            call plstop()
          endif
          call pincld(dnam)
        endif
        incf = .true.
        cprt = .false.

c     Perform inputs for initial conditions

      elseif(pcomp(titl(1:4),'init',4)) then
        call acheck(titl,yyy,15,80,80)
        titl( 1: 4) = yyy( 1: 4)
        titl(16:19) = yyy(16:19)
        setvar = vinput(yyy(31:75),45,td(1),3)
        call pinitl(dnam,td,errs)
        if(errs) call plstop()

c     Set title prints on/off

      elseif(pcomp(titl(1:4),'titl',4)) then
        call acheck(titl,yyy,15,80,80)
        read(yyy,1002,err=900,end=911) titl(1:4),titl(16:19)
        if(pcomp(titl(16:19),'off',3)) then
          prth = .false.
        else
          prth = .true.
        endif

c     Solution mode

      elseif(pcomp(titl(1:4),'inte',4)) then

        ior   = -abs(ior)
        evint = .true.
        intr  = .true.
        intx  = .true.
        cprt  = .true.
        call pltcur()
        go to 400

      elseif(pcomp(titl(1:4),'batc',4) .or.
     &       pcomp(titl(1:4),'macr',4)) then
        intr = .false.
        go to 400

c     Manual level set: 0 = basic; 1 = advanced; 2 = expert

      elseif(pcomp(titl(1:4),'manu',4)) then
        call acheck(titl,yyy,15,80,80)
        read(yyy,1003,err=900,end=911) titl(1:4),hlplev
        hlplev = max(-1,min(3,hlplev))

c     Mesh manipulations: Link and tie

c     Reset id list to link dof's on different nodes - set by node #

      elseif(pcomp(titl(1:4),'link',4)) then
        ndum = 0
        call setext('link', ndum, fext, .false.)
        call plinka(fext,'set','   ')
        lkflg = .true.
        tfl   = .true.

c     Reset id list to link dof's on different nodes - set by coord.

      elseif(pcomp(titl(1:4),'elin',4)) then
        ndum = 0
        call setext('elin', ndum, fext, .false.)
        call plinka(fext,'set','   ')
        leflg = .true.
        tfl   = .true.
        if(np(257).eq.0) then
          setvar = palloc(257,'ELINK',numnp*ndm,  1)
          do i = 1,ndm*numnp
            mr(np(257)+i-1) = 0
          end do ! i
          do i = 1,3
            dxlnk(i) = 0.0d0
          end do ! i
        endif

      elseif(pcomp(titl(1:4),'clin',4)) then
        call acheck(titl,yyy,15,80,256)
        read(yyy,1003,err=900,end=911) titl(1:4),(clnk(i),i=1,ndf)
        if(prt) then
          write(iow,2016) (i,clnk(i),i=1,ndf)
        endif
        lcflg = .true.
        tfl   = .true.
        if(np(79).eq.0) then
          setvar = palloc( 79,'IPOS ',numnp,  1)
          call pseqn(mr(np(79)),numnp)
        endif
        setvar = palloc(111,'TEMP1',numnp, 1)
        setvar = palloc(112,'TEMP2',numnp, 1)

        l1 = 1
        l2 = numnp
        l3 = 0
        l4 = 0
        l5 = 0
        td(2) = 0.0d0
        if(.not.gapfl) then
          gaptol = 1.d-3/sqrt(dble(max(1,numnp)))
        endif
        call tienod(mr(np(33)),hr(np(43)),mr(np(79)),mr(np(111)),
     &              mr(np(112)),mr(np(78)),ndm,nen,nen1,
     &              numnp,numel,l1,l2,l3,l4,l5,gaptol,td(2))

        setvar = palloc(112,'TEMP2',0, 1)
        setvar = palloc(111,'TEMP1',0, 1)

      elseif(pcomp(titl(1:4),'tie' ,3)) then
        go to 500

c     Parameter sets

      elseif(pcomp(titl(1:4),'para',4) .or.
     &       pcomp(titl(1:4),'cons',4)) then
        coflg = .true.
        call pconst(prt)

c     Rigid body and master/slave sets

c     Create rigid bodies

      elseif(pcomp(titl(1:4),'rigi',4)) then
        call rigidb(1,1,errs)
        if(errs) call plstop()
        tfl = .true.

c     Create list of joints

      elseif(pcomp(titl(1:4),'join',4)) then
        call rigidb(1,2,errs)
        if(errs) call plstop()

c     Rigid body loads

      elseif(pcomp(titl(1:4),'rloa',4)) then
        call rigidb(1,3,errs)
        if(errs) call plstop()

c     Rigid body boundary conditions

      elseif(pcomp(titl(1:4),'rbou',4)) then
        call rigidb(1,4,errs)
        if(errs) call plstop()

c     Rigid body displacements

      elseif(pcomp(titl(1:4),'rdis',4)) then
        call rigidb(1,5,errs)
        if(errs) call plstop()

c     Master/slave

      elseif(pcomp(titl(1:4),'mast',4)) then
        setvar = palloc( 167, 'RLINK', ndf*numnp, 1)
        call pmastr(mr(np(100)),mr(np(167)),mr(np(190)),hr(np(43)),prt)

c     Optimize profile

      elseif(pcomp(titl(1:4),'opti',4)) then
        call acheck(titl,yyy,15,80,80)
        read(yyy,1002,err=900,end=911) titl(1:4),dnam
        if(pcomp(dnam,'hoit',4)) then
          opthoit = .true.
        else
          opthoit = .false.
        endif
        call optid()
        tfl = .true.

c     Dictionary search for name of array

      elseif(pcomp(titl(1:4),'dict',4)) then
        call acheck(titl,yyy,15,80,80)
        read(yyy,1002,err=900,end=911) titl(1:4),dnam
        call pgetd(dnam,point,pln,pre,errs)
        write(*,*) point,pln,pre,errs

c     Loop start

      elseif(pcomp(titl(1:4),'loop',4)) then
        call acheck(titl,yyy,15,80,80)
        read(yyy,1002,err=900,end=911) titl(1:4),dnam
        call ploops(lp_in,dnam,1)

c     Loop end

      elseif(pcomp(titl(1:4),'next',4)) then
        call acheck(titl,yyy,15,80,80)
        read(yyy,1002,err=900,end=911) titl(1:4),dnam
        call ploops(lp_in,dnam,2)

c     File inputs

      elseif(pcomp(titl(1:4),'file',4)) then
        call acheck(titl,yyy,15,80,80)
        read(yyy,1002,err=900,end=911) titl(1:4),dnam
        lopen   = .true.
        lp_lun  = icl
        do while(lopen)
          lp_lun       = lp_lun + 1
          inquire(unit = lp_lun, opened = lopen)
        end do ! while
        ior       = lp_lun
        open(unit = ior, file = dnam, status = 'old')

c     Mesh generator outputs

      elseif(pcomp(titl(1:4),'tri2',4)) then
        call tridat(hr(np(43)),mr(np(33)),numnp,numel,nen,ndm,nen1,
     &              evint)

c     Unit multipliers for length, force and time

      elseif(pcomp(titl(1:4),'unit',4)) then
        call punits()

c     Data storage modes: Real/complex

      elseif(pcomp(titl(1:4),'*rea',4)) then
        cplxfl = .false.
        write(iow,2006)

      elseif(pcomp(titl(1:4),'*com',4)) then
        cplxfl = .true.
        write(iow,2007)

c     Contact data inputs

      elseif(pcomp(titl(1:4),'cont',4)) then
        backspace (ior)
        call plinka('cxi ','set','end')
        cxifl = .true.

c     Partition sets

      elseif(pcomp(titl(1:4),'part',4)) then

        npart = 0
        endwh = .false.
        do while (.not.endwh)

          npart  = npart + 1
          partfl = .false.

c         Read partition data
          do i = 1,ndf,16
            errck = pinput(td(1),min(16,ndf-i+1))
            jj = 0
            do j = i,min(i+15,ndf)
              jj           = jj + 1
              if(nint(td(jj)).gt.0) then
                ndfst(j,npart) = npart
                partfl  = .true.
              else
                ndfst(j,npart) = 0
              endif
            end do ! j
          end do ! i

          if(.not.partfl) then
            npart = npart - 1
            endwh = .true.
          endif
        end do ! while

c       Save maximum partition number

        maxpart = npart

c       Output stagger data

        if(prt) then
          write(iow,2002) (i,i=1,ndf)
          do j = 1,npart
            write(iow,2027) j,(ndfst(i,j),i=1,ndf)
          end do ! j
        endif

c       Allocate storage for stagger equation numbers

        setvar = palloc( 31,'ID   ',ndf*numnp*(npart+1), 1)
        idpt(1) = np(31)
        id31    = np(31)

        do j = 1,npart
          if(j.gt.1) then
            idpt(j) = np(31) + nneq*j
          endif
          call pzeroi(nqp , 5)
          do i = 1,ndf
            ndfp(i)      = ndfst(i,j)
            ndfg(i)      = ndfp(i)
            nqp(ndfp(i)) = ndfp(i)
          end do ! i
          do i = 1,4
            npart = nqp(i)
            if(npart.ge.1 .and. npart.le.4) then
              call partpt(npart,tflp(npart),.false.)
            elseif(npart.ne.0) then
              write(  *,3000) npart
              write(ilg,3000) npart
            endif
            ndlp(i) = 0
          end do ! i
        end do ! j

c       Set initial allocation for monolithic case

        do i = 1,ndf
          ndfst(i,5) = 5
        end do ! i
        npart   = 5
        idpt(5) = np(31)
        call partpt(npart,tflp(npart),.false.)
        ndlp(5) = 0

c     ODE order for each degree of freedom

      elseif(pcomp(titl(1:4),'orde',4)) then
        do i = 1,ndf,16
          errck = pinput(td(1),min(16,ndf-i+1))
          jj = 0
          do j = i,min(i+15,ndf)
            jj           = jj + 1
            ndfo(j) = nint(td(jj))
            ndog(j) = ndfo(jj)
          end do ! j
        end do ! i
        write(iow,2003) (i,ndfo(i),i=1,ndf)

c     Set list of Lagrange multiplier nodal variables

      elseif(pcomp(titl(1:4),'lagr',4)) then
        do i = 1,ndf,16
          errck = pinput(td(1),min(16,ndf-i+1))
          jj = 0
          do j = i,min(i+15,ndf)
            jj           = jj + 1
            ndfl(i) = nint(td(jj))
          end do ! j
        end do ! i
        write(iow,2024) (i,ndfl(i),i=1,ndf)
        lagrfl = .true.

c     Remarks to output file

      elseif(pcomp(titl(1:4),'rema',4)) then
        write(*,2008) titl(1:78)

c     Debug set

      elseif(pcomp(titl(1:4),'debu',4)) then
        debug = .true.

c     Domain inputs for parallel solutions

      elseif(pcomp(titl(1:4),'doma',4)) then

       errck = pinput(td,3)

       numpn  = nint(td(1))
       numtn  = nint(td(2))
       numteq = nint(td(3))

       write(iow,2023) numpn, numtn, numteq

       call pdomain(prt)

c     Stop execution

      elseif(pcomp(titl(1:4),'stop',4)) then
        if(evint) write(*,2004) fout
        if(abs(ior).eq.irdef) then
          if(prob_on) then   ! Multiple problem is active
            call pendprob
          endif
          call plstop()
        elseif(ior.eq.icf) then
          call pincld('end')
          incf = .false.
          intx = .false.
          cprt = evint
        endif

c     Return without a call to plstop

      elseif(pcomp(titl(1:4),'cntn',4)) then
        return

      endif

c     Read again

      go to 1

c     Start Problem: Read and print control information

100   newprob = .true.
      do i = 1,20
        l2 = 4*i
        l1 = l2 - 3
        head(i) = titl(l1:l2)
      end do ! i
      call pnewprob(1)
      go to 1

c     [mani] - Perform user manipulation commands

300   errs  = .true.
      j     = 0
      do while(errs .and. j.lt.26)
        j     = j + 1
        fnamr = fsav
        write(fext(3:3),'(a1)') char(96+j)
        call addext(fnamr,fext,128,8)
        inquire(file = fnamr, exist = errs)
      end do !
      call plinka(fext,'set','   ')
      go to 1

c     [batc]h and [inte]ractive execution
c      Establish profile of resulting equations for stiffness, mass, etc

400   if(.not.newprob) then
        write(  *,3001)
        write(ilg,3001)
        call plstop()
      elseif(intx .and. .not.intr .and. .not.incf) then
        write(  *,3002)
        write(ilg,3002)
        go to 1
      endif

      if(tfl) then

c       Check if nadd greater than nad

        if(nadd.gt.nad) then

          nst    = max(ndf*nen + nadd,1)
          setvar = palloc( 34,'LD   ',max(nen+1,7*nst*iif,21), 1)
          setvar = palloc( 35,'P    ',nst*3*iif*ipc          , 2)
          setvar = palloc( 36,'S    ',nst*nst*2*iif*ipc      , 2)
          setvar = palloc( 41,'UL   ',nst*14*ipc             , 2)

        endif

c       Merge boundary conditions after ties: forces & contact

        if(tief) then
          call pmerbc(nopart,prt)
c         tief = .false.
        endif

c       Compute active boundary loading after ties

        call pecmes()

c       Compute boundary nodes (after ties)

        call pextnd()

c       Rotational 5/6 degree-of-freedom checks

        if(frotas) then
          call rotred(mr(np(81)),mr(np(31)+nneq),mr(np(190)),ndf,numnp)
        endif

c       Allocate memory to store all possible equations
c       (include: rigid body dof x rigid body)
c       (include: rigid joints   x 6         )

        do j = 1,4
          neq = 0
          do i = 1,ndf
            if(ndfp(i).eq.j) then
              neq = neq + numnp
            endif
          end do ! i
          if(nrbprt.eq.j) then
            neq = neq + nrbdof*nrbody + numjts*6
          elseif(numjts.gt.0) then
            neq = neq + numjts*6
          endif
          if(neq.gt.0) then
            write(dnam,'(a2,i1)') 'JP',j
            setvar = palloc( 20+j, dnam, neq, 1)
          endif
        end do ! j

c       Set user commands

        do j = 1,10
          fext = 'u1a'
          if(j.lt.10) then
            write(fext(2:2),'(i1)') j
          else
            write(fext(2:2),'(i1)') 0
          endif
          if(usetfl(j)) then
            do l3 = 1,26
              write(fext(3:3),'(a1)') char(96+l3)
              fnamr =  fsav
              call addext(fnamr,fext,128,8)
              inquire(file = fnamr, exist = errs)
              if(errs) then
                call opnfil(fext,fnamr,-1,ios,prt)

c               Read data from file

                iorsv = ior
                ior   = ios

                do l1 = 0,36
                  do l2 = 1,26
                    vvsave(l2,l1) = vvv(l2,l1)
                  end do ! l2
                end do ! l1
                do l1 = 1,3
                  x0sav(l1) = x0(l1)
                end do ! l1
                oprt  = prt
                oprth = prth

                read(ior,1004) type,fincld(isf),irecrd(isf),prt,prth
                read(ior,1005) vvv
                read(ior,1005) tr,xr,trdet,x0
                read(ior,1006) ldnum,ldprp,spnum,ldflg,spflg

                call usetlib(j)

                close(ior,status='delete')
                ior   = iorsv

                do l1 = 0,36
                  do l2 = 1,26
                    vvv(l2,l1) = vvsave(l2,l1)
                  end do ! l2
                end do ! l1
                do l1 = 1,3
                  x0(l1) = x0sav(l1)
                end do ! l1
                prt  = oprt
                prth = oprth

              endif
            end do ! l3
          endif
        end do ! j

c       Set default partition data

        if(nopart) then
          nopart = .false.
          npart  = 1
          call partpt(npart,tflp(npart),.false.)
        endif

c       Construct interface mesh

        if(mchflg) then

c         Compute sparse structure for matrix

          setvar = palloc(111,'TEMP1', numnp+1, 1)         ! IC
          call optic(numnp,numel,numcels,nen,nen1,ncen,ncen1,
     &               mr(np(33)),mr(np(168)),mr(np(111)), l1)

          setvar = palloc(112,'TEMP2', l1, 1)              ! IP
          call pzeroi(mr(np(112)),l1)

          setvar = palloc(113,'TEMP3', numel+numcels, 1)   ! NNEL
          call opcon(numel,numcels,nen,nen1,ncen,ncen1,mr(np(33)),
     &                mr(np(168)),mr(np(111)),mr(np(112)),mr(np(113)))
          setvar = palloc(113,'TEMP3', 0, 1)

c         Set an upper bound for number of faces

          setvar = palloc(210,'MATCH', 50*numel, 1)

          call faceset(hr(np(25)),mr(np(32)),mr(np(33)),mr(np(111)),
     &                 mr(np(112)),mr(np(210)),l2, l4, .true.)

          setvar = palloc(210,'MATCH', 8*l2, 1)

          setvar = palloc(112,'TEMP2', 0, 1)
          setvar = palloc(111,'TEMP1', 0, 1)

c         Allocate storage for interface history

          setvar = palloc(214,'HINTE',max(1,l4),2)
          setvar = palloc(215,'HINT1',max(1,hnimax),2)
          setvar = palloc(216,'HINT2',max(1,hnimax),2)
          setvar = palloc(217,'HINT3',max(1,hni3max),2)

        endif
      endif ! tfl = .true.

c     Check if contact has been specified

      if(cxifl) then

c       Generate the contact surfaces, pairs and material sets

        call pextnd()
        fext = 'cxi'
        call pinpfl('PCONTR',fext, type, 1)
        call contact (1)
        call pinpfl('PCONTR',fext, type, 2)

c       Check contact surfaces for eliminated tied nodes

        if(tief) then
          call contact (312)
          tief = .false.
        endif

        call contact (313)

        cxifl = .false.
        tfl   = .true.

c       End of contact message

        write(iow,2025)
        write(ilg,2026)

      endif

c     Set initial history in elements and contact

      if(tfl) then

c       Perform check on mesh again to set final boundary codes

        setvar = palloc(111,'TEMP1',max(numnp*(ndf+1),1), 1)

c       Check contact elements

        call contact (316)

c       Check remaining elements

        call pidset(mr(np(111)),mr(np(32)),mr(np(240)),mr(np(31)+nneq),
     &              mr(np(190)),mr(np(33)),nie,nen,nen1,ndf,numnp,numel,
     &              nummat)
        setvar = palloc(111,'TEMP1',0, 1)

c       Compute side lengths of mesh for tolerancing

        call pdxmsh(mr(np(190)),hr(np(43)))

c       Set up stress history addresses

        call sethis(mr(np(32)),mr(np(33)),mr(np(181)),nie,nen,nen1,
     &              numel,nummat,prt)

c       Initialize history database items: Include 2 x for DG interface

        if(nhmax.gt.0) then
          setvar = palloc( 50,'NH1  ', nhmax*2,  2)
          setvar = palloc( 51,'NH2  ', nhmax*2,  2)
        endif
        if(nh3max.gt.0) then
          setvar = palloc( 52,'NH3  ', nh3max*2, 2)
        endif

        hflgu  = .true.
        h3flgu = .true.
        ctan(1) = 1.d0
        ctan(2) = 0.d0
        ctan(3) = 0.d0

        nreg    = -1   ! All regions active

c       Set initial element sizes & call element module

        hsize (1) = 0.0d0
        hsize (2) = 0.0d0
        hksize(1) = 0.0d0
        hksize(2) = 0.0d0

c       call formfe(   u  ,   dr ,   dr ,   dr
        call formfe(np(40),np(26),np(26),np(26),
     &             .false.,.false.,.false.,.false.,14,1,numel,1)

        if(max(hsize(1),hsize(2)).gt.0.0d0) then
          write(iow,2028) hsize
        endif
        if(max(hksize(1),hksize(2)).gt.0.0d0) then
          write(iow,2029) hksize
        endif

      endif

c     Determine current profile

      if(tfl) then
        do j = 0,nneq-1
          mr(np(31)+j) = -abs(mr(np(31)+j+nneq))
        end do ! j

        mxpro = 0
        mxneq = 0
        do j = 1,4
          if(.not.tflp(j)) then
            if(prt.and.ior.lt.0) write(  *,2001) j
            if(prt)              write(iow,2001) j
            id31 = idpt(j)
            do i = 0,nneq-1
              mr(id31+i) = -abs(mr(np(31)+i+nneq))
            end do ! i
            call partpt(j,tflp(j),.false.)
            mxprop(j) = 0
            mxneqp(j) = 0

c           Set current profile

            call profil(mr(np(20+j)),mr(np(34)),mr(id31),
     &                  mr(np(33)),1,prt)
            call profil(mr(np(20+j)),mr(np(34)),mr(id31),
     &                  mr(np(33)),2,prt)
            nqp(j)    = neq
            nqr(j)    = neqr
            mxprop(j) = max(mxprop(j),(mr(np(20+j)+neq-1))*ipc)
            mxneqp(j) = max(mxneqp(j),neq*ipc)
            mxpro     = max(mxpro,mxprop(j))
            mxneq     = max(mxneq,mxneqp(j))
          endif
        end do ! j
        if(pfeap_blk) then  ! Form when all equations are included
          call pidreset(mr(id31))
        endif

        tfl = .false.

c       End of mesh message

        write(iow,2009)
        write(ilg,2010)

      endif

c     Command language module for establishing solution algorithm
c     N.B. Program starts in "partition 1"

      if(initf) then
        call partpt(1,tflp(1),.false.)
      endif
      call pmacr(initf)
      go to 1

c     Tie nodes within tolerance of one another
c     [tie] - merge parts   with common coordinates
c     [tie regi r1 r2]: tie regions  r1 to r2 if nodes have same coord
c     [tie mate m1 m2]: tie material m1 to m2 if nodes have same coord
c     [tie node n1 n2]: tie nodes    n1 to n2 if nodes have same coord
c     [tie line l1 m2]: tie material m1 to m2 if nodes have same coord
c     [tie coor c1 c2 c3]: tie x(c1,*) = c2 =/- c3 (tolerance)

500   call acheck(titl,yyy,15,90,90)
      titl( 1: 4) = yyy( 1: 4)
      titl(16:19) = yyy(16:19)
      setvar = vinput(yyy(31:90),60,td(1),4)

c     Retrieve current boundary connection status

      if(.not.tief) then
        setvar = palloc( 79,'IPOS ',numnp,  1)
        call pseqn(mr(np(79)),numnp)
c       Save untied data
        call ptieix(mr(np(31)+nneq),mr(np(33)),mr(np(190)),hr(np(27)),1)
        tief = .true.
      endif

c     Turn off all ties

      if(pcomp(titl(16:18),'off',3)) then
c       Retrieve untied data
        call ptieix(mr(np(31)+nneq),mr(np(33)),mr(np(190)),hr(np(27)),2)
        tiefl = .true.
        tief  = .false.

c     Tie line elements to regions

      elseif(pcomp(titl(16:19),'line',4)) then
        l2 = max(    1,min(nummat,int(td(1))))
        call ptiend(mr(np(32)),mr(np(33)),mr(np(78)),mr(np(79)),
     &              hr(np(43)),l2,nie,nen,nen1,ndm,numel)
      else

        if(pcomp(titl(16:19),'node',4)) then
          l1 = max(    1,int(td(1)))
          l2 = min(numnp,int(td(2)))
          l5    = 0
          td(2) = 0.0d0
          if(.not.gapfl) then
            gaptol = 1.d-3/sqrt(dble(max(1,numnp)))
          endif
          write(iow,2011) l1,l2
        elseif(pcomp(titl(16:19),'regi',4)) then
          l1 = 1
          l2 = numnp
          l3 = max(    0,int(td(1)))
          l4 = min(mxreg,int(td(2)))
          l5 =-1
          if(.not.gapfl) then
            gaptol = 1.d-3/sqrt(dble(max(1,numnp)))
          endif
          write(iow,2012) l3,l4
        elseif(pcomp(titl(16:19),'mate',4)) then
          l1 = 1
          l2 = numnp
          l3 = max(     1,int(td(1)))
          l4 = min(nummat,max(1,int(td(2))))
          l5 =-2
          if(.not.gapfl) then
            gaptol = 1.d-3/sqrt(dble(max(1,numnp)))
          endif
          write(iow,2013) l3,l4
        elseif(pcomp(titl(16:19),'coor',4)) then
          l1 = 1
          l2 = numnp
          l3 = 1
          l5 =-3
          do j = 4,1,-1
            td(j+1) = td(j)
          end do ! j
          if(.not.gapfl) then
            gaptol = 1.d-3/sqrt(dble(max(1,numnp)))
          endif
          write(iow,2021) gaptol,(td(j),j=2,ndm+1)
        elseif(pcomp(titl(16:19),'tol',3)  .or.
     &         pcomp(titl(16:19),'gap',3)) then
          if(td(1).eq.0.0d0) then
            gaptol = 1.d-3/sqrt(dble(max(1,numnp)))
            gapfl  = .false.
          else
            gaptol = td(1)
            gapfl  = .true.
          endif
          go to 1
        elseif(pcomp(titl(16:19),'face',4)) then
          l1 = nint(td(1))
          l2 = nint(td(2))
          l3 = max(0,min(1,nint(td(3))))
          write(iow,2022) wfac(l3+1),l1, wfac(l3+1),l2
          if(.not.gapfl) then
            gaptol = 1.d-3/sqrt(dble(max(1,numnp)))
          endif
          call ptiefac(hr(np(43)),mr(np(33)),mr(np(79)),
     &                 l1,l2,l3,gaptol)
          tfl = .true.
          go to 1
        else
          l1 = 1
          l2 = numnp
          l3 = 0
          l4 = 0
          l5 = nint(td(1))
          if(.not.gapfl) then
            gaptol = 1.d-3/sqrt(dble(max(1,numnp)))
          endif
          if(l5.gt.0) then
            write(iow,2014) l5,td(2)
          else
            write(iow,2015)
          endif
        endif
        setvar = palloc(111,'TEMP1',numnp, 1)
        setvar = palloc(112,'TEMP2',numnp, 1)

        call tienod(mr(np(33)),hr(np(43)),mr(np(79)),mr(np(111)),
     &              mr(np(112)),mr(np(78)),ndm,nen,nen1,
     &              numnp,numel,l1,l2,l3,l4,l5,gaptol,td(2))

        setvar = palloc(112,'TEMP2',0, 1)
        setvar = palloc(111,'TEMP1',0, 1)
      endif
      tfl = .true.
      go to 1

c     Error treatments

900   call  errclr ('PCONTR')
      call plstop()

910   if(ior.eq.icf) then
        call pincld('end')
        incf = .false.
        intx = evint
        cprt = evint
        go to 1
      endif

911   call  endclr ('PCONTR',titl)
      call plstop()

c     Input/output formats

1001  format(a)

1002  format(a4,11x,a)

1003  format(a4,11x,15i15)

1004  format(a4,2x,a12,i8,2l5)

1005  format(4f20.0)

1006  format(3i8,2l3)

2001  format(/5x,'P a r t i t i o n',i4)

2002  format(/'   N o d a l   P a r t i o n   D a t a'/
     &       3x,'Partition',10(i3,'-dof':)/(12x,10(i3,'-dof':)))

2003  format(/'   N o d a l   M a x i m u m   PDE   O r d e r'/
     &        '      ndf   Order'/ (i8,i12))

2004  format(/' *End of <FEAP> solution,  File: ',a/1x)

2005  format(/' *End of INCLUDE solution, File: ',a/1x)

2006  format(//'  **Variable Storage REAL**'///)

2007  format(//'  **Variable Storage COMPLEX**'///)

2008  format(/' ',a/)

2009  format(/'  ',28('-'),' END OF MESH INPUTS ',28('-')/)
2010  format(/'  ',33('-'),' END OF MESH INPUTS ',33('-')/)

2011  format(/5x,'Tie nodes from',i8,' to ',i8/1x)

2012  format(/5x,'Tie from region',i4,' to region',i4/1x)

2013  format(/5x,'Tie from material',i4,' to material',i4/1x)

2014  format(/5x,'Tie: direction =',i3,' X =',1p,1e12.5/1x)

2015  format(/5x,'Tie all nodes with common coordinates'/1x)

2016  format(/'   C o o r d i n a t e     L i n k s'/
     &        '      ndf        Link'/ (i8,i12))

2018  format(/'  ',33('-'),' START OF FEAP LOG ',34('-')/)

2021  format(/5x,'Tie all nodes with common coordinates'/
     &       10x,'Tol  = ',1p,1e12.4 /
     &       10x,'x_1  = ',1p,1e12.4:/
     &       10x,'x_2  = ',1p,1e12.4:/
     &       10x,'x_3  = ',1p,1e12.4:)

2022  format(/5x,'Tie 4-node elements with common faces'/
     &       10x,a,' 1 =',i4/ 10x,a,' 2 =',i4)

2023  format(//5x,'DOMAIN Information',//
     &        10x,'Number of partition nodes    : ',i7/
     &        10x,'Number of total problem nodes: ',i7,/
     &        10x,'Number of total problem eqns.: ',i7,/)

2024  format(/'   N o d a l   L a g r a n g e   M u l t p l i e r'/
     &        '      ndf  Multiplier'/'             Number'/(i8,i12))

2025  format(/'  ',28('-'),' END OF CONTACT INPUTS ',25('-')/)
2026  format(/'  ',33('-'),' END OF CONTACT INPUTS ',30('-')/)

2027  format(i10,10i7/10x,10i7)

2028  format(/5x,'E l e m e n t   S i z e   V a l u e s'/
     &       10x,'h-minimum =',1p,1e12.4/
     &       10x,'h-maximum =',1p,1e12.4)

2029  format(/5x,'E l e m e n t   K n o t   S i z e   V a l u e s'/
     &       10x,'hk-minimum =',1p,1e12.4/
     &       10x,'hk-maximum =',1p,1e12.4)

3000  format(/' *ERROR* PCONTR: Partition wrong: Npart =',i3)

3001  format(/' *ERROR* PCONTR: Attempt to solve problem before mesh'/
     &        '         input.  Check for error on FEAP or BINAry',
     &        ' record.'/1x)

3002  format(/' *ERROR* PCONTR: Can not do BATCH execution from this',
     &        ' mode.'/
     &        '         Do INTERACTIVE or put in INCLUDE file.'/1x)

3003  format(/' *ERROR* Include file: ',a,' does not exist')

      end
