c$Id:$
      subroutine contact (csw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add set of 'ifpck' to false (off)                03/05/2007
c       2. Add 'nummat' to 'rstprf' call                    21/07/2007
c       3. Only allocate 'IAD' once                         15/01/2008
c       4. Add 316 to set contact nodes for parallel        25/01/2013
c          Calls cidset to do surface node sets of np(111).
c       5. Add 317 to output contact data to files          04/02/2013
c       6. Add 318 to plot contact surface nodes            08/02/2013
c       7. Allocate np(224) under csw.eq.14                 09/02/2014
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0
c               Anna Haraldsson           January 1, 1998            1.1
c               Robert Taylor              March 15, 1998            1.2
c               Robert Taylor           February 25, 2001            1.3
c               Robert Taylor           December 31, 2002            1.4

c      Acronym: CONTact LIBrary

c      Purpose: Management of all contact activities

c      Inputs:
c         csw     - Contact switch

c      Outputs:
c                 - Perform requested activity
c-----[--.----+----.----+----.-----------------------------------------]
c     CONTACT SWITCH CONTROL TABLE

c          x   = defined in proper section
c          #   = defined in #      section
c          -   = use not allowed  -> return with no warning
c          .   = still undefined  -> return with no warning
c          I   = Internal use only


c       0 <= CSW <=  99         --> Call from FORMFE
c                                   Activity like for continuum element

c     100 <= CSW <= 199         --> Direct call for preliminary activity

c     200 <= CSW                --> Direct call for equivalent activity

c       0   -               200   x Show element info
c       1   x                       Input data
c       2   x                       Check data
c       3   2    103   x            Form stiffness / check geometry
c       4   -               204   2 Print stresses
c       5   -
c       6   -               206   x Form residual
c       7   -
c       8   .
c       9   .
c      10   x                       Perform augmentation
c      11   .
c      12   .
c      13   .
c      14   2                       Initialize history variables
c      15   .
c      16   .
c      17   .
c      18   .
c      19   .
c      20   .
c      21   .
c      22   .
c      23   .               223 --> Direct call for numerical tangent


c     300 <= CSW                --> Non standard calls

c     300   x                       Init contact flags for new problem
c     301   x                       Time step update
c     302   x                       Back-up to the beginning of the step
c     303   x                       Dump of contact arrays
c     304   x                       Set check of gap status
c     305   x                       Plot contact geometry
c     306   x                       Read data for restart
c     307   x                       Save data for restart
c     308   x                       Set profile & range to plot variable
c     309   x                       Switch contact on/off during solve.
c     310   x                       Reset penalty during solve.
c     311   x                       Compute average stresses & strains
c     312   x                       Adjust contact element on tied nodes
c     313   x                       Activate history variables
c     314   x                       Update element lagrange multipliers
c     315   x                       Output contact surface facet list
c     316   x                       Set contact element dofs active
c     317   x                       Output contact data sets
c     318   x                       Set contact nodes for plots


c     400 <= CSW                --> Special internal call

c     400   1                       Perform "once" activity for problem
c     403 103                       Reset profile for active contacts
c     408 308                       Plot contours of a contact variable

c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_contac.h'
      include  'c_dict.h'
      include  'c_ccp.h'
      include  'c_chp.h'
      include  'c_keyh.h'
      include  'c_tanfl.h'
      include  'c_values.h'
      include  'allotd.h'
      include  'cdata.h'
      include  'chdata.h'
      include  'compac.h'
      include  'compas.h'
      include  'complx.h'
      include  'idptr.h'
      include  'iodata.h'
      include  'iofile.h'
      include  'ndata.h'
      include  'part0.h'
      include  'print.h'
      include  'pview.h'
      include  'sdata.h'
      include  'setups.h'
      include  'ssolve.h'
      include  'umac1.h'

      include  'pointer.h'
      include  'comblk.h'

      logical   pcomp,setvar,palloc,errck,gettxtd,vinput,surfl
      character tx(2)*15,tname*5, ch4*4
      integer   csw, lcc0,kp,npair,kd,ncel,surtyp, typ
      real*8    td(15)

      save

c     Save debug information

      call cdebug0 ('CONTACT',csw)

c-----[--.----+----.----+----.-----------------------------------------]
c     Skip if contact is off

      if ((csw.ne.  1) .and. (csw.ne.300) .and.
     &    (csw.ne.309) .and. (csw.ne.400)) then
        if (.not.ifct) return
      endif

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from PCONTR for data input

      if (csw.eq.1) then

c       Set global flag (ON - OFF - DEBU)

        errck = gettxtd(tx,2,td,3,'skip')

c       Contact data exist but contact is set off

        if(pcomp(tx(2),'off',3)) then
          ifct = .false.

c       Wrong contact option

        elseif ((.not.pcomp(tx(2),'on',2))    .and.
     &          (.not.pcomp(tx(2),'    ',4))  .and.
     &          (.not.pcomp(tx(2),'debu',4))) then
          write (  *,3001) tx(2),'on','debug'
          write (iow,3001) tx(2),'on','debug'
          write (ilg,3001) tx(2),'on','debug'
          call plstop()

c       Contact data exist and contact active

        else
          ifct = .true.

c         Check debug mode

          if (pcomp(tx(2),'debu',4)) then
            ifdb = .true.
            indb = nint(td(1))
            call cdebug0 ('debf',nint(td(2)))
            call cdebug0 ('** PCONTR',0)
            call cdebug0 ('CONTACT',csw)
            call cdebug  (nint(td(3)),'debf',' ')
          endif

c         Initialize contact algorithm

          call continit ()

c         Rename pair element type

          uct = cis(typ(3,1))
          call celmt01 (ndm,ndf,hr,hr,0,npair,hr,hr,hr,hr,hr,mr,mr,
     &                  hr,hr,hr,hr,hr,ch4,ch4)
          cis(typ(3,1)) = uct

          uct = cis(typ(3,2))
          call celmt02 (ndm,ndf,hr,hr,0,npair,hr,hr,hr,hr,hr,mr,mr,
     &                  hr,hr,hr,hr,hr,ch4,ch4)
          cis(typ(3,2)) = uct

          uct = cis(typ(3,3))
          call celmt03 (ndm,ndf,hr,hr,0,npair,hr,hr,hr,hr,hr,mr,mr,
     &                  hr,hr,hr,hr,hr,ch4,ch4)
          cis(typ(3,3)) = uct

          uct = cis(typ(3,4))
          call celmt04 (ndm,ndf,hr,hr,0,npair,hr,hr,hr,hr,hr,mr,mr,
     &                  hr,hr,hr,hr,hr,ch4,ch4)
          cis(typ(3,4)) = uct

          uct = cis(typ(3,5))
          call celmt05 (ndm,ndf,hr,hr,0,npair,hr,hr,hr,hr,hr,mr,mr,
     &                  hr,hr,hr,hr,hr,ch4,ch4)
          cis(typ(3,5)) = uct

          uct = cis(typ(3,6))
          call celmt06 (ndm,ndf,hr,hr,0,npair,hr,hr,hr,hr,hr,mr,mr,
     &                  hr,hr,hr,hr,hr,ch4,ch4)
          cis(typ(3,6)) = uct

          uct = cis(typ(3,7))
          call celmt07 (ndm,ndf,hr,hr,0,npair,hr,hr,hr,hr,hr,mr,mr,
     &                  hr,hr,hr,hr,hr,ch4,ch4)
          cis(typ(3,7)) = uct

          uct = cis(typ(3,8))
          call celmt08 (ndm,ndf,hr,hr,0,npair,hr,hr,hr,hr,hr,mr,mr,
     &                  hr,hr,hr,hr,hr,ch4,ch4)
          cis(typ(3,8)) = uct

          uct = cis(typ(3,9))
          call celmt09 (ndm,ndf,hr,hr,0,npair,hr,hr,hr,hr,hr,mr,mr,
     &                  hr,hr,hr,hr,hr,ch4,ch4)
          cis(typ(3,9)) = uct

          uct = cis(typ(3,10))
          call celmt10 (ndm,ndf,hr,hr,0,npair,hr,hr,hr,hr,hr,mr,mr,
     &                  hr,hr,hr,hr,hr,ch4,ch4)
          cis(typ(3,10)) = uct

          uct = cis(typ(3,11))
          call celmt11 (ndm,ndf,hr,hr,0,npair,hr,hr,hr,hr,hr,mr,mr,
     &                  hr,hr,hr,hr,hr,ch4,ch4)
          cis(typ(3,11)) = uct

          uct = cis(typ(3,12))
          call celmt12 (ndm,ndf,hr,hr,0,npair,hr,hr,hr,hr,hr,mr,mr,
     &                  hr,hr,hr,hr,hr,ch4,ch4)
          cis(typ(3,12)) = uct

          uct = cis(typ(3,13))
          call celmt13 (ndm,ndf,hr,hr,0,npair,hr,hr,hr,hr,hr,mr,mr,
     &                  hr,hr,hr,hr,hr,ch4,ch4)
          cis(typ(3,13)) = uct

          uct = cis(typ(3,14))
          call celmt14 (ndm,ndf,hr,hr,0,npair,hr,hr,hr,hr,hr,mr,mr,
     &                  hr,hr,hr,hr,hr,ch4,ch4)
          cis(typ(3,14)) = uct

          uct = cis(typ(3,15))
          call celmt15 (ndm,ndf,hr,hr,0,npair,hr,hr,hr,hr,hr,mr,mr,
     &                  hr,hr,hr,hr,hr,ch4,ch4)
          cis(typ(3,15)) = uct

          uct = cis(typ(3,16))
          call celmt16 (ndm,ndf,hr,hr,0,npair,hr,hr,hr,hr,hr,mr,mr,
     &                  hr,hr,hr,hr,hr,ch4,ch4)
          cis(typ(3,16)) = uct

          uct = cis(typ(3,17))
          call celmt17 (ndm,ndf,hr,hr,0,npair,hr,hr,hr,hr,hr,mr,mr,
     &                  hr,hr,hr,hr,hr,ch4,ch4)
          cis(typ(3,17)) = uct

          uct = cis(typ(3,18))
          call celmt18 (ndm,ndf,hr,hr,0,npair,hr,hr,hr,hr,hr,mr,mr,
     &                  hr,hr,hr,hr,hr,ch4,ch4)
          cis(typ(3,18)) = uct

          uct = cis(typ(3,19))
          call celmt19 (ndm,ndf,hr,hr,0,npair,hr,hr,hr,hr,hr,mr,mr,
     &                  hr,hr,hr,hr,hr,ch4,ch4)
          cis(typ(3,19)) = uct

          uct = cis(typ(3,20))
          call celmt20 (ndm,ndf,hr,hr,0,npair,hr,hr,hr,hr,hr,mr,mr,
     &                  hr,hr,hr,hr,hr,ch4,ch4)
          cis(typ(3,20)) = uct

          call cnts2d  (ndm,ndf,hr,hr,0,npair,hr,hr,mr,mr,
     &                  hr,hr,hr,hr,ch4,ch4)

          call cnts3d  (ndm,ndf,hr,hr,0,npair,hr,mr,mr,
     &                  hr,hr,hr,hr,ch4,ch4)

          call cptpnd  (ndm,ndf,hr,hr,0,npair,hr,mr,mr,hr,hr,hr,hr,
     &                  ch4,ch4)

          call cntrnd  (ndm,ndf,hr,hr,0,npair,hr,hr,mr,mr,
     &                  hr,hr,hr,ch4,ch4)

          call ctied2d (ndm,ndf,hr,hr,0,npair,hr,hr,mr,mr,
     &                  hr,hr,hr,ch4,ch4)

c         Check amount of contact data

          call pnumc ()

c         Set up control table for commands

          call cccontab (lcc0)

c         Allocate memory for command control tables

          setvar = palloc( 131,'C0   ',lcc0,2)

c         Allocate memory of material data(max c_lmv=50 input constants)
c         N.B. Allocate at least one word to prevent screen warning from
c              PALLOC.

          setvar = palloc( 132,'CM   ',max(1,numcm*c_lmv),2)

c         Allocate memory for contact elements storage
c         REMARK - element storage extended after reading elements

          setvar = palloc( 150,'CTE15',numnp, 1)

c         Input contact data from file

          call pcont(hr(ccp(1)),hr(ccp(2)),hr(ccp(3)),hr(ccp(8)),
     &               hr(ccp(9)),hr(ccp(10)),hr(np(132)))

c         Extension of the element array

          setvar = palloc( 150,'CTE15',    0, 1)

c         Allocate history address correspondence vector for all pairs

          setvar = palloc( 134,'HIC  ',(c_lp1+c_lp3)*numcp,   1)

c         Pair loop

          do kp = 1,numcp
            npair = kp

c           Set default values for pair

            call defaultp (npair,hr(ccp(1)),hr(ccp(2)),hr(ccp(3)),
     &                     hr(np(132)))

c           Activate history variables for the specific type of contact
c           and determine history area for CH1, CH2, CH3

            call contlib (csw,npair,hr(ccp(1)),hr(ccp(2)),hr(ccp(3)),
     &                    mr(np(133)),hr(np(132)),hr(chp(1)),hr(chp(2)),
     &                    hr(chp(3)),mr(np(134)))
          end do ! kp

        endif

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from FORMFE for data check                      csw =   2
c     Called from FORMFE to compute stiffness and residual   csw =   3
c     Called from FORMFE to compute residual                 csw =   6
c     Called from FORMFE for augmentation                    csw =  10
c     Called from FORMFE for initialize history variables    csw =  14
c     Called from PMACR1 to print contact status             csw = 204
c     Called from PTIMPL to compute residual                 csw = 206
c     Called from PMACR1 to compute averages stresses        csw = 311
c     Called from UPDATE to update Lagrange multipliers      csw = 314

      elseif ((csw.eq.  2) .or.
     &        (csw.eq.  3) .or.
     &        (csw.eq.  6) .or.
     &        (csw.eq. 10) .or.
     &        (csw.eq. 14) .or.
     &        (csw.eq.204) .or.
     &        (csw.eq.206) .or.
     &        (csw.eq.311) .or.
     &        (csw.eq.314)) then

        do kp = 1, numcp
          npair = kp
          call contlib (csw,npair,hr(ccp(1)),hr(ccp(2)),hr(ccp(3)),
     &                  mr(np(133)),hr(np(132)),hr(chp(1)),hr(chp(2)),
     &                  hr(chp(3)),mr(np(134)))
        end do ! kp

c       Dump for debugging (only in init phase)

        if (csw.eq.14) then

c         Allocate array for Lagrange multipliers

          if(lagrm .and. np(224).eq.0) then
            setvar = palloc( 224,'IAD  ',3*numnp,1)
          endif

          call pnewpr()  ! Reset profile in case of changes

          call cdebug (0,'surf','after initialization')
          call cdebug (0,'mate','after initialization')
          call cdebug (0,'pair','after initialization')
          call cdebug (0,'ICS' ,'after initialization')
          call cdebug (0,'CM'  ,'after initialization')
          call cdebug (0,'CH'  ,'after initialization')
        endif

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from PMACR1 to set up geometrical variables     csw = 103
c     Called from OPTID  to set up geometrical variables     csw = 203

      elseif (csw.eq.103 .or. csw.eq.203) then

c       Set geometrical relationships

        do kp = 1,numcp
          npair = kp
          call contlib (103,npair,hr(ccp(1)),hr(ccp(2)),hr(ccp(3)),
     &                  mr(np(133)),hr(np(132)),hr(chp(1)),hr(chp(2)),
     &                  hr(chp(3)),mr(np(134)))
        end do ! kp

c       Count active contacts for compact storage and optimization

        iccom   = 1
        ncen    = 0
        numcels = 0
        do kp = 1,numcp
          npair = kp
          call contlib(403,npair,hr(ccp(1)),hr(ccp(2)),hr(ccp(3)),
     &                 mr(np(133)),hr(np(132)),hr(chp(1)),hr(chp(2)),
     &                 hr(chp(3)),mr(np(134)))
        end do ! kp

c       Assemble single array with active contact element connections

        if(numcels.gt.0) then

          if(lagrm) then
            ncen1 = ncen + 2
          else
            ncen1 = ncen
          endif

          setvar = palloc(168,'IXC  ', ncen1*numcels, 1)

          iccom   = 2
          numcels = 0
          do kp = 1,numcp
            npair = kp
            call contlib(403,npair,hr(ccp(1)),hr(ccp(2)),hr(ccp(3)),
     &                   mr(np(133)),hr(np(132)),hr(chp(1)),hr(chp(2)),
     &                   hr(chp(3)),mr(np(134)))
          end do ! kp
        endif

        if(csw.eq.203) return

c       Optimize profile and reset equation numbers

        if(optflg .or. lagrm) then

          if(optflg) call optid()

          call revopt(mr(np(89)),mr(np(89)+numnp),numnp)
          npair = ndf*numnp
          do kp = 0,npair-1
            mr(id31+kp) = mr(np(31)+kp+npair)
          end do ! kp

          call seteq(mr(id31),mr(np(99)),mr(np(182)),mr(np(100)),
     &               mr(np(101)),.false.)

c         Consider extra Lagrange multiplier equations

          if(np(224).eq.0) then
            setvar = palloc( 224,'IAD  ',3*numnp,1)
          endif

c         Zero IAD array

          do kp = 0,3*numnp-1
            mr(np(224)+kp) = 0
          end do ! kp

c         Set equations

          iccom   = 3
          do kp = 1,numcp
            npair = kp
            call contlib (403,npair,hr(ccp(1)),hr(ccp(2)),hr(ccp(3)),
     &                 mr(np(133)),hr(np(132)),hr(chp(1)),hr(chp(2)),
     &                 hr(chp(3)),mr(np(134)))
          end do ! kp

c         Adjust equation numbers when elements have no multiplers

          if(np(211).eq.0) then
            call newnum (mr(np(224)),mr(id31),mr(np(89)),
     &                   mr(np(89)+numnp),mr(np(99)),mr(np(101)))

c         Adjust equation numbers when elements have multiplers

          else
            call newnuml(mr(np(224)),mr(id31),mr(np(211)),mr(np(212)),
     &                   mr(np(168)),mr(np(89)+numnp),
     &                   mr(np(99)),mr(np(101)))
          endif

c         Reset size of solution array if necessary

          setvar = palloc(26,'DR   ',max(neq,numnp*max(ndf*ipc,ndm)),2)

          write(tname,'(a2,i1)') 'JP',npart
          setvar = palloc( 20+npart, tname, neq, 1)

        endif

c       Using standard Feap solvers

        if(solver) then

c         Form compact storage data

          if(compfl) then

c           Allocate storage arrays and retrieve existing ones

            setvar = palloc( 93,'OINC ', 2*neq+1, 1)
            setvar = palloc(111,'TEMP1', max(ndf*numnp,neq), 1)

c           Compute compressed profile array

            call pzeroi(mr(np(111)),max(ndf*numnp,neq))
            call elcnt(numel,nen,nen1,mr(id31),mr(np(33)),
     &                 mr(np(111)), 1)

c           Consider contact elements connected to each node

            if(numcels.gt.0) then
              call elcnt(numcels,ncen,ncen1,mr(id31),
     &                   mr(np(168)),mr(np(111)),-1)
            endif
            call sumcnt(mr(np(111)),neq,kp)

c           Set list of elements connected to nodes

            setvar = palloc(112,'TEMP2', kp, 1)
            call pelcon(numel,nen,nen1,mr(np(33)),mr(id31),
     &                  mr(np(111)),mr(np(112)),kp,1)
            if(numcels.gt.0) then
              call pelcon(numcels,ncen,ncen1,mr(np(168)),mr(id31),
     &                    mr(np(111)),mr(np(112)),kp, -1)
            endif

c           Allocate arrays for compressed profile solution

            if(np(94).ne.0) then
              setvar = palloc( 94,'OINO ',0, 1)
            endif
            setvar = palloc(113,'TEMP3',neq, 1)
            call comproa(numnp,nen,nen1,ndf,mr(np(33)),mr(id31),
     &                   mr(np(111)),mr(np(112)),mr(np(113)),kp,
     &                   kbycol,kdiag,kall)
            setvar = palloc(113,'TEMP3',  0, 1)
            setvar = palloc( 94,'OINO ', kp, 1)
            call comprob(numnp,nen,nen1,ndf,mr(np(33)),mr(id31),
     &                   mr(np(111)),mr(np(112)),mr(np(94)),mr(np(93)),
     &                   kbycol,kdiag,kall)
            kcmplx = kp

c           Delete unused arrays

            setvar = palloc(112,'TEMP2', 0, 1)
            setvar = palloc(111,'TEMP1', 0, 1)

c           Direct: Block profile

            write(tname,'(a4,i1)') 'TANG',npart
            if(ittyp.eq.-1) then
              ncel = kp

c           Direct: Block sparse (symmetric)

            elseif(ittyp.eq.-2) then

c             Destroy old tangent matrix

              if(np(npart).ne.0) then
                setvar = palloc(npart,tname, 0 , 2)
              endif
              if(np(67).ne.0) then
                setvar = palloc(67,'SPTAN', 0 , 1)
                nspo = 0
              endif

              call sortjc(mr(np(93)),mr(np(94)),neq)

c             Set 'p' and 'ip' to equation numbers

              setvar = palloc( 47, 'PNTER',neq,  1)
              setvar = palloc( 48, 'INVPT',neq,  1)
              do kd = 1,neq
                mr(np(47)+kd-1) = kd
                mr(np(48)+kd-1) = kd
              end do ! kd
              ncel = 0
              ddom = .true.

c           Direct: Profile (incore) & Other

            else
              ncel = 0
            endif

c           Allocate memory for added terms

            setvar = palloc(npart,tname, kp+neq+ncel, 2)

c           Reset pointers

            na     = np(npart)
            nnr    = kp  + neq
            nau    = na  + neq
            nal    = nau + ncel
          endif

c         Reset profile of continuum discretization

          call rstprf(mr(np(20+npart)),mr(np(34)),mr(id31),mr(np(33)),
     &                mr(np(32)),mr(np(99)),mr(np(100)),mr(np(101)),
     &                ndf,nen1,nen,neq,numnp,numel,nummat)

c         Reset profile for contact

          iccom   = 4
          do kp = 1,numcp
            npair = kp
            call contlib (403,npair,hr(ccp(1)),hr(ccp(2)),hr(ccp(3)),
     &                    mr(np(133)),hr(np(132)),hr(chp(1)),hr(chp(2)),
     &                    hr(chp(3)),mr(np(134)))
          end do ! kp

c         Set new pointer for profile of global stiffness

          call nwprof (mr(np(20+npart)),neq)

        endif ! solver = .true.

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from PMACR5 to show element information         csw = 200

      elseif (csw.eq.200) then
        write (iow,2000)
        if (ior.lt.0) write (*,2000)

c       Scan all available contact elements

        do kd = 1,c_ncel
          ncel = kd
          call contlib (csw,ncel,hr(ccp(1)),hr(ccp(2)),hr(ccp(3)),
     &                  mr(np(133)),hr(np(132)),hr(chp(1)),hr(chp(2)),
     &                  hr(chp(3)),mr(np(134)))
        end do ! kp

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from PCONTR to set flags for new problem        csw = 300

      elseif (csw.eq.300) then
        ifct    = .false.  ! Contact solution OFF
        ifdb    = .false.  ! Contact debuging OFF
        lagrm   = .false.  ! Lagrange multiplier solution OFF
        shakefl = .false.  ! Explicit rhake  flag
        rattlfl = .false.  ! Explicit rattle flag
        ifpck   = .false.  ! Initial penetration check OFF
        numcp   =  0       ! Number of contact pairs

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from PMACR2 to perform time step update         csw = 301

      elseif (csw.eq.301) then
        call creshis (hr(chp(1)),hr(chp(2)))

c       Give to the user the possibility to reset variables

        do kp = 1, numcp
          npair = kp
          call contlib (csw,npair,hr(ccp(1)),hr(ccp(2)),hr(ccp(3)),
     &                   mr(np(133)),hr(np(132)),hr(chp(1)),hr(chp(2)),
     &                   hr(chp(3)),mr(np(134)))
        end do ! kp

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from AUTBAC: Automatic back-up to step beginning csw = 302

      elseif (csw.eq.302) then
        call creshis (hr(chp(2)),hr(chp(1)))

c       Give user possibility to reset variables

        do kp = 1, numcp
          npair = kp
          call contlib (csw,npair,hr(ccp(1)),hr(ccp(2)),hr(ccp(3)),
     &                   mr(np(133)),hr(np(132)),hr(chp(1)),hr(chp(2)),
     &                   hr(chp(3)),mr(np(134)))
        end do ! kp

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from OUTARY for dump of contact arrays          csw = 303

      elseif (csw.eq.303) then
        kp = 0
        call cdebug (kp,yyy(16:19),'interactive')

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from PMACR3 to set flag for gap status check    csw = 304

      elseif (csw.eq.304) then
        if     (pcomp(yyy(16:19),'chec',4)) then
         ifistgn = .false.
        elseif (pcomp(yyy(16:19),'noch',4)) then
         ifistgn = .true.
        endif

        ifchist = .false.
        do kp = 1, numcp
          npair = kp
          call contlib (csw,npair,hr(ccp(1)),hr(ccp(2)),hr(ccp(3)),
     &                   mr(np(133)),hr(np(132)),hr(chp(1)),hr(chp(2)),
     &                   hr(chp(3)),mr(np(134)))
        end do ! kp

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from PPLOTF for plot of contact                 csw = 305

      elseif (csw.eq.305) then
        call pdefm(hr(np(43)),hr(np(40)),cs,ndm,ndf,numnp,hr(np(53)))
        call plopen
        do kp = 1, numcp
          npair = kp
          call contlib (csw,npair,hr(ccp(1)),hr(ccp(2)),hr(ccp(3)),
     &                   mr(np(133)),hr(np(132)),hr(chp(1)),hr(chp(2)),
     &                   hr(chp(3)),mr(np(134)))
        end do ! kp

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from RESTRT for read  data for restart          csw = 306
c     Called from RESTRT for write data for restart          csw = 307

      elseif ((csw.eq.306) .or. (csw.eq.307) .and. numcp.gt.0) then
        call crestrt (csw,hr(chp(1)),hr(chp(2)),hr(chp(3)))
        if(csw.eq.306) then
          write(iow,2001)
          if(ior.lt.0) write(*,2001)
        else
          write(iow,2002)
          if(ior.lt.0) write(*,2002)
        endif

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from PPLOTF for contour plot of variable        csw = 308

      elseif (csw.eq.308) then
        call pdefm(hr(np(43)),hr(np(40)),cs,ndm,ndf,numnp,hr(np(53)))
        call plopen

c       Set plot range

        do kp = 1, numcp
          npair = kp
          call contlib (csw,npair,hr(ccp(1)),hr(ccp(2)),hr(ccp(3)),
     &                   mr(np(133)),hr(np(132)),hr(chp(1)),hr(chp(2)),
     &                   hr(chp(3)),mr(np(134)))
        end do ! kp

c       Plot variable

        do kp = 1, numcp
          npair = kp
          call contlib (408,npair,hr(ccp(1)),hr(ccp(2)),hr(ccp(3)),
     &                   mr(np(133)),hr(np(132)),hr(chp(1)),hr(chp(2)),
     &                   hr(chp(3)),mr(np(134)))
        end do ! kp

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from PMACR3 to switch contact on/off            csw = 309
c     WARNING: Contact switched 'on' only if set at start of problem.

      elseif (csw.eq.309) then
        if    (pcomp(yyy(16:19),'on  ',4)) then
          ifct = .true.
          if(ior.lt.0) then
            write(*,2003)
          endif
        elseif(pcomp(yyy(16:19),'off ',4)) then
          ifct = .false.
          if(ior.lt.0) then
            write(*,2004)
          endif
        endif

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from PMACR3 to reset penalty                    csw = 310

      elseif (csw.eq.310) then
        if    (pcomp(yyy(16:19),'pena',4)) then
          errck = vinput(yyy(31:60),30,td,2)
          npair = nint(td(1))
          if(npair.gt.0 .and. npair.le.numcp) then
            if(ior.lt.0) write(*,2005) npair,td(2)
            cvalue(1) = td(2)
            cvalue(2) = dble(npair)
            call contlib (csw,npair,hr(ccp(1)),hr(ccp(2)),hr(ccp(3)),
     &                    mr(np(133)),hr(np(132)),hr(chp(1)),hr(chp(2)),
     &                    hr(chp(3)),mr(np(134)))
          endif
        elseif(pcomp(yyy(16:19),'fric',4)) then
          if(ior.lt.0) write(*,2006)
          iffron = .true.
          ifchist = .false.
          do kp = 1,numcp
            npair = kp
            call contlib (csw,npair,hr(ccp(1)),hr(ccp(2)),hr(ccp(3)),
     &                    mr(np(133)),hr(np(132)),hr(chp(1)),hr(chp(2)),
     &                    hr(chp(3)),mr(np(134)))
          end do ! kp
        elseif(pcomp(yyy(16:19),'nofr',4)) then
          if(ior.lt.0) write(*,2007)
          iffron = .false.
          ifchist = .false.
          do kp = 1,numcp
            npair = kp
            call contlib (csw,npair,hr(ccp(1)),hr(ccp(2)),hr(ccp(3)),
     &                    mr(np(133)),hr(np(132)),hr(chp(1)),hr(chp(2)),
     &                    hr(chp(3)),mr(np(134)))
          end do ! kp
        endif

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from PCONTR to check nodes eliminated by ties   csw = 312

      elseif (csw.eq.312) then

        call cksurf(hr(ccp(1)),mr(np(133)),mr(np(79)))

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from PCONTR to set storage for history variable csw = 313

      elseif (csw.eq.313) then

c         Pair loop

          surfl = .false.
          do kp = 1,numcp
            npair = kp
            call cptest(hr(ccp(3)),npair,surfl,surtyp)
          end do ! kp

c         Build list of facets connected ot each contact node

          if(surfl) then
            if(surtyp.eq.0) then   ! Used for current 3-d contact
              setvar = palloc( 192,'SURFP', numcs ,  1)
              call csurface0(hr(ccp(1)),mr(np(133)),mr(np(192)))
            elseif(surtyp.eq.1) then ! Form with split arrays
              setvar = palloc( 192,'SURFP', 3*numcs+3 ,  1)
              call csurface1(hr(ccp(1)),mr(np(133)),mr(np(192)))
            endif
          endif

c         Pair loop

          do kp = 1,numcp
            npair = kp

c           Activate history variables for the specific type of contact
c           and determine history area for CH1, CH2, CH3

            call contlib (csw,npair,hr(ccp(1)),hr(ccp(2)),hr(ccp(3)),
     &                    mr(np(133)),hr(np(132)),hr(chp(1)),hr(chp(2)),
     &                    hr(chp(3)),mr(np(134)))

c           Store history managment information

            call stohman (mr(np(134)),npair,hr(ccp(3)))
          end do ! kp

c         Scan all contact elements to reset flag for "only once action"

          do kd = 1,c_ncel
            ncel = kd
            call contlib (400,ncel,hr(ccp(1)),hr(ccp(2)),hr(ccp(3)),
     &                    mr(np(133)),hr(np(132)),hr(chp(1)),hr(chp(2)),
     &                    hr(chp(3)),mr(np(134)))
          end do ! kp

c         Allocation history array CH (ch1+ch2+ch3)

          setvar = palloc( 135,'CH   ',max(tlch1*2 + tlch3,1),   2)

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from POUTM to output contact lists                csw = 315

      elseif (csw.eq.315) then

        write(ios,2008)

c       Output surface data to saved mesh

        call coutm (hr(ccp(1)),hr(ccp(2)),hr(ccp(3)),mr(np(133)),
     &              mr(np(134)))

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from PCONTR to activate contact nodes             csw = 316

      elseif (csw.eq.316) then

        call cidset(hr(ccp(1)),mr(np(133)), mr(np(111)), ndm, ndf )

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from UOUTDOM to output parallel contact lists     csw = 317

      elseif (csw.eq.317) then

c       Output surface data to saved mesh

        call cpoutm (hr(ccp(1)),hr(ccp(2)),hr(ccp(3)),mr(np(133)),
     &               mr(np(134)))

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from PNUMNA to activate plot of contact nodes     csw = 318

      elseif (csw.eq.318) then

        call cidset(hr(ccp(1)),mr(np(133)), mr(np(111)), ndm, 0 )

c-----[--.----+----.----+----.-----------------------------------------]
c     Controlled not allowed or not necessary entries

      elseif ((csw.eq.0)  .or.
     &        (csw.eq.4)  .or.
     &        (csw.eq.5)  .or.
     &        (csw.eq.7)      ) then
        continue

c-----[--.----+----.----+----.-----------------------------------------]
c     For any other still not evaluated case let element decide

      else

        ifchist = .false.
        do kp = 1,numcp
          npair = kp
          call contlib (csw,npair,hr(ccp(1)),hr(ccp(2)),hr(ccp(3)),
     &                  mr(np(133)),hr(np(132)),hr(chp(1)),hr(chp(2)),
     &                  hr(chp(3)),mr(np(134)))
        end do ! kp

      endif

c     Formats

2000  format(/'   A v a i l a b l e    C o n t a c t   E l e m e n t s',
     &       /)

2001  format(10x,'Contact data input')
2002  format(10x,'Contact data output')

2003  format(/5x,'Contact is ON'/1x)
2004  format(/5x,'Contact is OFF'/1x)

2005  format(/5x,'Contact Penalty: Surface =',i5,
     &           ' Penalty =',1p,1e12.4/1x)

2006  format(/5x,'Contact friction is ON'/1x)
2007  format(/5x,'Contact friction is OFF'/1x)

2008  format(/'END'//'CONTACT')

3001  format(' *ERROR* CONTACT: Illegal Contact Option: ',a/
     &       '         OPTIONS: ',a,', ',a,' or blank.'/)

      end
