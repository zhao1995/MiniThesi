c$Id:$
      subroutine continit ()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Correct 'exte' label to BLENd not BLOCk          15/01/2007
c       2. Set number table rows from c_nr0 in 'c_0.h'      16/02/2007
c       3. Add initial penetration setting of coordinates   03/05/2007
c       4. Add axisymmetry option for node to segment       09/07/2010
c       5. Add 'nurb' and 'tspl'ine surface definitions     02/08/2011
c          Add 'side' to sub-command options.
c       6. Add 'norm'al, 'plus','minus' options to surface  02/10/2011
c       7. Add 'pert'urbed lagrangian solution option       18/10/2011
c       8. Add 'inte'rpolation type for nurbs, t-spline     06/20/2012
c       9. Add 'nfac' to prescribe NURB face number for 3-d 19/09/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0
c               R.L. Taylor              October 15, 2001            1.1

c      Acronym: CONTact INITialization

c      Purpose: Init variables and define command strings
c          Current Command List:

c      Inputs :

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_comnd.h'
      include  'c_contac.h'
      include  'c_dict.h'
      include  'c_dicti.h'
      include  'c_keyh.h'
      include  'c_tole.h'
      include  'iofile.h'
      include  'machnc.h'

      integer   typ,fep,opp,scp,sop,kr,kc,k,ncom,ncs

      save

      call cdebug0 ('  continit',0)

c-----[--.----+----.----+----.-----------------------------------------]
c     CONTACT COMMANDS

      cwd(1)  = 'surf'    ! SURFace  input data
      cwd(2)  = 'mate'    ! MATErial input data
      cwd(3)  = 'pair'    ! PAIR     input data
      cwd(4)  = 'auto'    ! AUTO
      cwd(5)  = 'read'    ! READ from file
      cwd(6)  = 'save'    ! SAVE data on file
      cwd(7)  = 'tran'    ! TRANsformation array
      cwd(8)  = 'help'    ! HELP
      cwd(9)  = 'unu1'    ! UNUsed 1
      cwd(10) = 'unu2'    ! UNUsed 2
      cwd(11) = 'unu3'    ! UNUsed 3
      cwd(12) = 'end'     ! END

c-----[--.----+----.----+----.-----------------------------------------]
c     CONTACT COMMANDS HELP FLAG

      ced(1)  = 3         ! SURF
      ced(2)  = 3         ! MATE
      ced(3)  = 3         ! PAIR
      ced(4)  = 0         ! AUTO
      ced(5)  = 0         ! SAVE
      ced(6)  = 0         ! HELP
      ced(7)  = 0         ! TRAN
      ced(8)  = 3         ! READ
      ced(9)  = 3         ! UNU1
      ced(10) = 3         ! UNU2
      ced(11) = 3         ! UNU3
      ced(12) = 3         ! END

c-----[--.----+----.----+----.-----------------------------------------]
c     CONTACT COMMANDS STRINGS QUANTIFICATION

c     Clean arrays

      do k = 1, c_ncc
        nty(k) = 0   ! # of defined TYPES
        nfe(k) = 0   ! # of defined FEATURES
        nop(k) = 0   ! max # of OPTIONS for each feature
        nsc(k) = 0   ! # of defined SUB-COMMANDS
        nso(k) = 0   ! max # of sub-command OPTIONS for each sub-command
      end do

c     For each command define the # of specified
c     TYPES
c       FEATURES,OPTIONS
c         SUB-COMMANDS,SUB-COMMAND OPTIONS

c     To add a definition increase the defined number, and then add
c     the string in the vector CIS

c-----[--.----+----.----+----.-----------------------------------------]
c     Command # 1 - SURF

      ncom = 1
      nty(ncom)  = 9   ! # of surface TYPES defined
      nfe(ncom)  = 0
      nop(ncom)  = 0
      nsc(ncom)  = 9   ! # of defined SUB-COMMANDS
      nso(ncom)  = 5   ! max # sub-command OPTIONS for each sub-command

c     Command # 2 - MATE

      ncom = 2
      nty(ncom)  = 3   ! # of material TYPES defined
      nfe(ncom)  = 0
      nop(ncom)  = 0
      nsc(ncom)  = 0
      nso(ncom)  = 0

c     Command # 3 - PAIR

      ncom = 3
      nty(ncom)  = 26   ! # of pair TYPES defined
      nfe(ncom)  = 11   ! # of defined FEATURES
      nop(ncom)  =  7   ! max # of OPTIONS for each feature
      nsc(ncom)  =  0
      nso(ncom)  =  0

c     Command # 4 - AUTO

      ncom = 4
      nty(ncom)  = 3   ! # of auto TYPES defined
      nfe(ncom)  = 0
      nop(ncom)  = 0
      nsc(ncom)  = 0
      nso(ncom)  = 0

c     Command # 5  - READ
c     Command # 6  - SAVE
c     Command # 7  - TRAN
c     Command # 8  - HELP
c     Command # 9  - UNU1
c     Command # 10 - UNU2
c     Command # 11 - UNU3
c     Command # 12 - END

c     Compute offsets to access CIS vector

      ofsfe = 0
      do k = 1, c_ncc
        ofsfe = ofsfe + nty(k)
      end do
      ofsop = ofsfe
      do k = 1, c_ncc
        ofsop = ofsop + nfe(k)
      end do
      ofssc = ofsop
      do k = 1, c_ncc
        ofssc = ofssc + nfe(k)*nop(k)
      end do
      ofsso = ofssc
      do k = 1, c_ncc
        ofsso = ofsso + nsc(k)
      end do

c     Check if string dictionary have enough space

      ncs = ofsso
      do k = 1, c_ncc
        ncs = ncs + nsc(k)*nso(k)
      end do

      if (ncs.gt.c_ncs) then
        write (  *,3001) ncs,c_ncs
        write (ilg,3001) ncs,c_ncs
        call plstop()
      endif

c-----[--.----+----.----+----.-----------------------------------------]
c     CONTACT COMMANDS STRINGS DEFINITION
c     How to define strings:
c     TYPE definitions

c     cis(typ(#command,#type)) = 'string'

c     FEATURE definitions

c     cis(fep(#command,#feature)) = 'string'

c     OPTIONS definitions

c     cis(opp(#command,#feature,#option)) = 'string'

c     SUB-COMMANDS definitions

c     cis(scp(#command,#sub-command)) = 'string'

c     SUB-COMMAND OPTIONS definitions

c     cis(sop(#command,,#sub-command,#sub-command option)) = 'string'

c-----[--.----+----.----+----.-----------------------------------------]

c     Initialize strings to blank

      do ncom = 1,4
        do k = 1,nty(ncom)
          cis(typ(ncom,k)) = '    '
        end do ! k
        do kc = 1,nfe(ncom)
          cis(fep(ncom,kc)) = '    '
          do k = 1,nop(ncom)
            cis(opp(ncom,kc,k)) = '    '
          end do ! k
        end do ! kc
        do kc = 1,nsc(ncom)
          cis(scp(ncom,kc)) = '    '
          do k = 1,nso(ncom)
            cis(sop(ncom,kc,k)) = '    '
          end do ! k
        end do ! kc
      end do ! ncom

c-----[--.----+----.----+----.-----------------------------------------]
c     Command # 1 - SURF
c     TYPES definition

      ncom = 1
      cis(typ(ncom,1))   = 'line'  ! LINE 2d elmt with 2 or more nodes
      cis(typ(ncom,2))   = 'tria'  ! TRIAngular 3d elmt: 3 or more nodes
      cis(typ(ncom,3))   = 'quad'  ! QUADrilateral 3d elmt: 8/more node
      cis(typ(ncom,4))   = 'beam'  ! BEAM 3d elmt with 2 or more nodes
      cis(typ(ncom,5))   = 'poin'  ! POINt element with 1 node.
      cis(typ(ncom,6))   = 'rigi'  ! RIGId element surface function
      cis(typ(ncom,7))   = 'nurb'  ! NURBs element surface function
      cis(typ(ncom,8))   = 'tspl'  ! TSPLine element surface function
      cis(typ(ncom,9))   = 'part'  ! PART name for elm surface fn

c     FEATURES definition
c     This definition used to allocate equivalent amount of columns
c     in command control tables

c     cis(fep(ncom,1))   = '    '

c     OPTIONS definitions

c     cis(opp(ncom,1,1))   = '    '

c     SUB-COMMANDS definitions
c     This definition used to switch data reading to corresponding
c     subroutine CREL01 ..... CREL10
c     10 subroutines are actually available

      cis(scp(ncom,1))   = 'face'  !  FACEt, manual input
      cis(scp(ncom,2))   = 'bloc'  !  BLOCk elem automatic generation
      cis(scp(ncom,3))   = 'blen'  !  BLENd elem automatic generation
      cis(scp(ncom,4))   = 'regi'  !  REGIon elem automatic generation
      cis(scp(ncom,5))   = 'func'  !  FUNCtion for rigid surfaces
      cis(scp(ncom,6))   = 'side'  !  SIDE number for nurbs patches
      cis(scp(ncom,7))   = 'norm'  !  NORMal direction over-ride
      cis(scp(ncom,8))   = 'nfac'  !  FACE number for nurbs 3D Block
      cis(scp(ncom,9))   = 'patc'  !  PATCH area for nurbs blocks

c     SUB-COMMANDS OPTIONS definitions

c     cis(sop(ncom,1,1))   = '    '

      cis(sop(ncom,2,1))   = 'gap'   ! BLOCk    GAP set    option
      cis(sop(ncom,2,2))   = 'segm'  ! BLOCk    SEGMent    option
      cis(sop(ncom,2,3))   = 'pola'  ! BLOCk    POLAr      option
      cis(sop(ncom,2,4))   = 'cart'  ! BLOCk    CARTesian  option
      cis(sop(ncom,2,5))   = 'regi'  ! BLOCk    REGIon     option

      cis(sop(ncom,3,1))   = 'gap'   ! BLENd    GAP set    option
      cis(sop(ncom,3,2))   = 'segm'  ! BLENd    SEGMent    option
      cis(sop(ncom,3,3))   = 'exte'  ! BLEND    EXTErnal   option

      cis(sop(ncom,5,1))   = 'cyli'  ! FUNCtion CYLInder   option
      cis(sop(ncom,5,2))   = 'sphe'  ! FUNCtion SPHEre     option
      cis(sop(ncom,5,3))   = 'cart'  ! FUNCtion CARTesian  option
      cis(sop(ncom,5,4))   = 'plan'  ! FUNCtion PLANe      option
      cis(sop(ncom,5,5))   = 'poly'  ! FUNCtion POLYnomial option

      cis(sop(ncom,7,1))   = 'plus'  ! PART     NORMal     PLUS
      cis(sop(ncom,7,2))   = 'minu'  ! PART     NORMal     MINUS
c-----[--.----+----.----+----.-----------------------------------------]
c     Command # 2 - MATE
c     TYPES definition
c     This definition used to switch data reading to corresponding
c     subroutine CRMAT01 ..... CRMAT10
c     10 subroutines are actually available

      ncom = 2
      cis(typ(ncom,1))   = 'stan'  ! STANdard "rigid" matl: Coulomb fr
      cis(typ(ncom,2))   = 'nlfr'  ! "rigid" material: NonLinear FRic
      cis(typ(ncom,3))   = 'user'  ! User defined material: data in cm0

c     FEATURES definition
c     This definition used to allocate an equivalent amount of columns
c     in the command control tables

c     cis(fep(ncom,1))   = '    '

c     OPTIONS definitions

c     cis(opp(ncom,1,1))   = '    '

c     SUB-COMMANDS definitions

c     cis(scp(ncom,1))   = '    '

c     SUB-COMMANDS OPTIONS definitions

c     cis(sop(ncom,1,1))   = '    '

c-----[--.----+----.----+----.-----------------------------------------]
c     Command # 3 - PAIR
c     TYPES definition
c     This definition used to switch the solution to the corresponding
c     contact element CELMT01 ..... CELMT20
c     10 contact elements are actually available

      ncom = 3
      cis(typ(ncom, 1))   = 'cel1'  ! User contact CELMT01
      cis(typ(ncom, 2))   = 'cel2'  ! User contact CELMT02
      cis(typ(ncom, 3))   = 'cel3'  ! User contact CELMT03
      cis(typ(ncom, 4))   = 'cel4'  ! User contact CELMT04
      cis(typ(ncom, 5))   = 'cel5'  ! User contact CELMT05
      cis(typ(ncom, 6))   = 'cel6'  ! User contact CELMT06
      cis(typ(ncom, 7))   = 'cel7'  ! User contact CELMT07
      cis(typ(ncom, 8))   = 'cel8'  ! User contact CELMT08
      cis(typ(ncom, 9))   = 'cel9'  ! User contact CELMT09
      cis(typ(ncom,10))   = 'ce10'  ! User contact CELMT10

      cis(typ(ncom,11))   = 'ce11'  ! User contact CELMT11
      cis(typ(ncom,12))   = 'ce12'  ! User contact CELMT12
      cis(typ(ncom,13))   = 'ce13'  ! User contact CELMT13
      cis(typ(ncom,14))   = 'ce14'  ! User contact CELMT14
      cis(typ(ncom,15))   = 'ce15'  ! User contact CELMT15
      cis(typ(ncom,16))   = 'ce16'  ! User contact CELMT16
      cis(typ(ncom,17))   = 'ce17'  ! User contact CELMT17
      cis(typ(ncom,18))   = 'ce18'  ! User contact CELMT18
      cis(typ(ncom,19))   = 'ce19'  ! User contact CELMT19
      cis(typ(ncom,20))   = 'ce20'  ! User contact CELMT20

      cis(typ(ncom,21))   = 'ntos'  ! Node  To Segment
      cis(typ(ncom,22))   = 'ptop'  ! Point To Point
      cis(typ(ncom,23))   = 'nton'  ! Node  To Node
      cis(typ(ncom,24))   = 'ptor'  ! Point To Rigid
      cis(typ(ncom,25))   = 'ntor'  ! Node  To Rigid
      cis(typ(ncom,26))   = 'tied'  ! Tied surface to surface

c     FEATURES definition
c     This definition used to allocate an equivalent amount of columns
c     in the command control tables

      cis(fep(ncom,1))   = 'swit'  ! Switch on/off for features
      cis(fep(ncom,2))   = 'solm'  ! SOLution Method
      cis(fep(ncom,3))   = 'deta'  ! DETection Algorithm
      cis(fep(ncom,4))   = 'mate'  ! MATErials for the pair
      cis(fep(ncom,5))   = 'augm'  ! AUGMentation
      cis(fep(ncom,6))   = 'tole'  ! TOLerances of the Pair
      cis(fep(ncom,7))   = 'adhe'  ! ADHEsion limit (0 = infinite)
      cis(fep(ncom,8))   = 'pene'  ! PENEtration checking (default off)
      cis(fep(ncom,9))   = 'axis'  ! AXISym. anal. on/off (default off)
      cis(fep(ncom,10))  = 'inte'  ! INTErpolation type (default lagr)
      cis(fep(ncom,11))  = 'quad'  ! QUADrature for surface n1 x n2

c     OPTIONS definitions

      cis(opp(ncom,1,1))   = 'off'   ! SWIT - Switch off part/all pair
      cis(opp(ncom,1,2))   = 'on'    ! SWIT - Switch on  part/all pair
      cis(opp(ncom,1,3))   = 'timf'  ! SWIT - Switch on/off with time fn

      cis(opp(ncom,2,1))   = 'pena'  ! SOLM - PENAlty solution method
      cis(opp(ncom,2,2))   = 'lagm'  ! SOLM - LAGrange Multiplier method
      cis(opp(ncom,2,3))   = 'croc'  ! SOLM - CROss-Constraints method
      cis(opp(ncom,2,4))   = 'cons'  ! SOLM - CONStraint elimination
      cis(opp(ncom,2,5))   = 'shak'  ! SOLM - SHAKe  explicit algorithm
      cis(opp(ncom,2,6))   = 'ratt'  ! SOLM - RATTle explicit algorithm
      cis(opp(ncom,2,7))   = 'pert'  ! SOLM - PERTurb lagrangian method

      cis(opp(ncom,3,1))   = 'basi'  ! DETA - BASIc mode (default)
      cis(opp(ncom,3,2))   = 'semi'  ! DETA - SEMI-fixed mode
      cis(opp(ncom,3,3))   = 'rigi'  ! DETA - RIGId contact surface

      cis(opp(ncom,4,1))   = '    '  ! MATE -

      cis(opp(ncom,5,1))   = 'off'   ! AUGM - Set augmentation off
      cis(opp(ncom,5,2))   = 'basi'  ! AUGM - BASIc mode (default)
      cis(opp(ncom,5,3))   = 'hset'  ! AUGM - History SET augmentation
      cis(opp(ncom,5,4))   = 'lise'  ! AUGM - LImit Stiff. Enhance mode
      cis(opp(ncom,5,5))   = 'smau'  ! AUGM - SMArt mode U??

      cis(opp(ncom,6,1))   = 'none'  ! TOLE - No argument
      cis(opp(ncom,6,2))   = 'pene'  ! TOLE - PENEtration
      cis(opp(ncom,6,3))   = 'open'  ! TOLE - OPENing of gap
      cis(opp(ncom,6,4))   = 'outs'  ! TOLE - OUT of Segment

      cis(opp(ncom,7,1))   = 'infi'  ! ADHE - infinite
      cis(opp(ncom,7,2))   = 'stre'  ! ADHE - stress value given

      cis(opp(ncom,8,1))   = 'on  '  ! PENE - checking on
      cis(opp(ncom,8,2))   = 'off '  ! PENE - checking off

      cis(opp(ncom,9,1))   = 'on  '  ! PENE - checking on
      cis(opp(ncom,9,2))   = 'off '  ! PENE - checking off

      cis(opp(ncom,10,1))  = 'lagr'  ! INTE - type (default)
      cis(opp(ncom,10,2))  = 'nurb'  ! INTE - NURBS
      cis(opp(ncom,10,3))  = 'tspl'  ! INTE - Tspline

c     SUB-COMMANDS definitions

c     cis(scp(ncom,1))   = '    '

c     SUB-COMMANDS OPTIONS definitions

c     cis(sop(ncom,1,1))   = '    '

c-----[--.----+----.----+----.-----------------------------------------]
c     Command # 4 - AUTO
c     TYPES definition
c     This definition is used to perform auto generation of surface data

      ncom = 4
      cis(typ(ncom,1))   = 'surf'  ! Generate surfaces automatically
      cis(typ(ncom,2))   = 'pair'  ! Generate pair table automatically
      cis(typ(ncom,3))   = '    '  ! Equivalent to 'surf'

c     FEATURES definition
c     This definition used to allocate an equivalent amount of columns
c     in the command control tables

c     cis(fep(ncom,1))   = '    '

c     OPTIONS definitions

c     cis(opp(ncom,1,1))   = '    '

c     SUB-COMMANDS definitions

c     cis(scp(ncom,1))   = '    '

c     SUB-COMMANDS OPTIONS definitions

c     cis(sop(ncom,1,1))   = '    '

c-----[--.----+----.----+----.-----------------------------------------]
c     CONTACT OFFSETS for surfaces and materials data

      ofssurf = 1
      ofsmate = 1
      isgp1   = 1
      isgp3   = 1

c-----[--.----+----.----+----.-----------------------------------------]
c     CONTACT COMMAND CONTROL TABLES definitions
c     To store new data increase existing values or define new ones

c     Command table structure

c     1 column  for each feature listed
c     N columns for the user (specify here)
c     M columns for the system (specify here)
c       (the first one for the type declaration data)
c     Clean arrays

      do k = 1, c_ncc
        nc0(k) = 0           ! # of positive columns
        of0(k) = 0           ! offset for first system column
        nuc(k) = 0           ! # of user columns
        n0c(k) = 0           ! # of system columns
      end do

c     # of rows for all the control tables

      nr0 = c_nr0

c     # of columns storage for the features defined above by vector nfe

c     # of columns storage for the user

      nuc(1) = 1             ! surf
      nuc(2) = 1             ! mate
      nuc(3) = 1             ! pair

c     # of columns storage for the system

      n0c(1) = 1             ! surf
      n0c(2) = 1             ! mate
      n0c(3) = 1             ! pair

c-----[--.----+----.----+----.-----------------------------------------]
c     HISTORY VARIABLES
c     Set default string for history variables in the dictionary
c     for CH1 & CH2, and CH3

      do kr = 1, c_lp1
        do kc = 1, c_ncel
          w1(kr,kc) = '--------'
        end do
      end do
      do kr = 1, c_lp3
        do kc = 1, c_ncel
          w3(kr,kc) = '--------'
        end do
      end do

c-----[--.----+----.----+----.-----------------------------------------]
c     SET YOUR ESTIMATE FOR ZERO HERE
c     REMARK - if you want to avoid trouble
c              chose it bigger than the machine zero

      zero = 100.d0*epmac

c-----[--.----+----.----+----.-----------------------------------------]
c     Init length of history variables vector

      tlch1 = 0
      tlch3 = 0

c-----[--.----+----.----+----.-----------------------------------------]
c     Set geom ck flag true (check open/close at each iteration

      ifistgn = .true.
      ifchist = .true.
c-----[--.----+----.----+----.-----------------------------------------]

3001  format ('*ERROR* CONTINIT:'/
     &        '        Requested size of string dictionary',i5/
     &        '        Available                          ',i5/
     &        '        Change parameter c_ncs in "c_0.h" include file')

      end
