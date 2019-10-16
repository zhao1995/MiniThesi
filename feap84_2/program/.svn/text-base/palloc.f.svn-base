c$Id:$
      logical function palloc(num,name,length,precis)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add Parallel arrays                              02/11/2006
c       2. Move 'NFORC' from position 244 to 256            02/11/2006
c       3. Initialize ilist to 600, set nlist > list        19/01/2007
c       4. Remove extra blank character after XADJN         17/04/2007
c       5. Add 257, 'ELINK' alloction                       25/12/2007
c       6. Add 258, 'GUVAL' alloction                       27/03/2008
c       7. Add 259, 'RVELM' alloction                       27/03/2008
c       8. Add 260,261, 'FRVEL', 'SRVEL' allocation         27/07/2008
c       9. Add 262, 'IDSAV' allocation                      27/07/2008
c      10. Add 263, 'NURBS' allocation                      11/11/2008
c      11. Add 264, 'LNURB' allocation                      25/11/2008
c      12. Add 265 to 267 'LDTAB/LDNOD/LDVAL' allocation    03/01/2009
c      13. Add 268 'SPINS' for spin data                    09/03/2009
c      14. Add 269 'RVEMA' for rve material numbers         12/04/2009
c          Add 270 'RVESD' for rve send order
c      15. Add 271-275 for KNOTS,NSIDE,NTRAN,NBSID,LKSID    15/04/2009
c      16. Add 276-280 for I_LIN,X_lin,U_LIN,N_LIN,T_LIN    07/07/2009
c      17. Add 281-282 for S_LIN and P_LIN stress plots     27/07/2009
c      18. Add 283     for E_LIN of eigenvectors            29/07/2009
c          Add 284     for ID_LN of equations.
c      19. Add 285     for LKNOTs definition (T-spline)     17/08/2009
c      20. Add 286     for NELEments for F_bar submesh use  21/08/2009
c      21. Add 287-297 for T-Spline and Plot arrays         11/11/2010
c      22. Move 271-175 to 298-302; add ESPIN & ESPTR       03/01/2011
c      23. Add 303-305 for History plot arrays              04/01/2012
c      24. Add 206     for H_LIN nurbs history plots        12/01/2012
c      25. Add 307 'IPTAN' for sparse integer array         01/05/2012
c      26. Add 308-310 currently unused parameters          20/07/2012
c          Add 311-320 for temporary plot variables
c      27. Add 308 'KNOTL', 309 'SIDEL' for NURBS meshes    12/11/2012
c          310 'BLOKL', 311 KTNUM and 312 NBLOK integer data
c          MOVE 'PTEMn' up 10 places to 321-330.
c      28. Add 274,275 for'TRIAD','LTRIA' for 3-d bc rots   11/02/2013
c      29. Add 331 to 333 for Hill-Mandel arrays            13/04/2013
c      30. Add 319 'KNOTE' for knot element data            18/11/2013
c      31. Add 334 & 335 'CPTMX' & 'C_EMX' for mixed model  21/11/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Define, delete, or resize a dictionary entry.
c               Pointer defined for integer (single) and real
c               (double precision arrays).

c      Inputs:
c         num        - Entry number for array (see below)
c         name       - Name of array          (see below)
c         length     - Length of array defined: =0 for delete
c         precis     - Precision of array: 1 = integers; 2 = reals

c      Outputs:
c         np(num)    - Pointer to first word of array in memory.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    list
      parameter (list = 335)

      include   'allotd.h'
      include   'allotn.h'
      include   'pointer.h'

      logical    ualloc
      character  name*(*)
      integer    i, num,length,precis

      save

c     Set active arrays for FEAP into list

      data   (nlist(i),i=1,list)/

     &         'TANG1', 'TANG2', 'TANG3', 'TANG4', 'UTAN1', 'UTAN2',
     &         'UTAN3', 'UTAN4', 'CMAS1', 'CMAS2', 'CMAS3', 'CMAS4',

     &         'LMAS1', 'LMAS2', 'LMAS3', 'LMAS4', 'DAMP1', 'DAMP2',
     &         'DAMP3', 'DAMP4', 'JP1  ', 'JP2  ', 'JP3  ', 'JP4  ',

c     Solution arrays

c              'TANG1',     !     1: Symm. tangent, partition 1
c              'TANG2',     !     2: Symm. tangent, partition 2
c              'TANG3',     !     3: Symm. tangent, partition 3
c              'TANG4',     !     4: Symm. tangent, partition 4

c              'UTAN1',     !     5: Unsym tangent, partition 1
c              'UTAN2',     !     6: Unsym tangent, partition 2
c              'UTAN3',     !     7: Unsym tangent, partition 3
c              'UTAN4',     !     8: Unsym tangent, partition 4

c              'CMAS1',     !     9: Consist. mass, partition 1
c              'CMAS2',     !    10: Consist. mass, partition 2
c              'CMAS3',     !    11: Consist. mass, partition 3
c              'CMAS4',     !    12: Consist. mass, partition 4

c              'LMAS1',     !    13: Diagonal mass, partition 1
c              'LMAS2',     !    14: Diagonal mass, partition 2
c              'LMAS3',     !    15: Diagonal mass, partition 3
c              'LMAS4',     !    16: Diagonal mass, partition 4

c              'DAMP1',     !    17: Symm. damping, partition 1
c              'DAMP2',     !    18: Symm. damping, partition 2
c              'DAMP3',     !    19: Symm. damping, partition 3
c              'DAMP4',     !    20: Symm. damping, partition 4

c              'JP1  ',     !    21: Profile pointer, partition 1
c              'JP2  ',     !    22: Profile pointer, partition 2
c              'JP3  ',     !    23: Profile pointer, partition 3
c              'JP4  ',     !    24: Profile pointer, partition 4

c     Mesh arrays

     &         'D    ', 'DR   ', 'F    ', 'F0   ', 'FPRO ', 'FTN  ',
     &         'ID   ', 'IE   ', 'IX   ', 'LD   ', 'P    ', 'S    ',
     &         'SLODI', 'T    ', 'TL   ', 'U    ', 'UL   ', 'VEL  ',
     &         'X    ', 'XL   ', 'ANG  ', 'ANGL ',


c              'D    ',     !    25: Material parameters
c              'DR   ',     !    26: Residual/reactions

c              'F    ',     !    27: Nodal load/displacement, current
c              'F0   ',     !    28: Nodal load/displace pattern, base
c              'FPRO ',     !    29: DOF proportional load numbers
c              'FTN  ',     !    30: Nodal load/displacement, current

c              'ID   ',     !    31: Equation numbers/boundary conds
c              'IE   ',     !    32: Element assembly information
c              'IX   ',     !    33: Element connection data

c              'LD   ',     !    34: Element local/global eq numbers

c              'P    ',     !    35: Element vector

c              'S    ',     !    36: Element array
c              'SLODI',     !    37: Surface load data

c              'T    ',     !    38: Nodal temperatures
c              'TL   ',     !    39: Nodal temperaturese, element

c              'U    ',     !    40: Nodal solutions/increments
c              'UL   ',     !    41: Nodal solutions/increments,element

c              'VEL'  ,     !    42: Nodal transient solution values

c              'X    ',     !    43: Nodal coordinates
c              'XL   ',     !    44: Nodal coordinates, element

c              'ANG  ',     !    45: Nodal angles
c              'ANGL ',     !    46: Nodal angles, element

c     Sparse solution data

     &         'PNTER', 'INVPT',

c              'PNTER',     !    47: Pointer array
c              'INVPT',     !    48: Inverse pointer array

c     History data

     &         'H    ', 'NH1  ', 'NH2  ', 'NH3  ',

c              'H    ',     !    49: Element history parameters
c              'NH1  ',     !    50: Element history data at t_n
c              'NH2  ',     !    51: Element history data at t_n+1
c              'NH3  ',     !    52: Element history data, time ind

c     Plot data

     &         'CT   ', 'FCIX ', 'FCZM ', 'FIDP ', 'NDER ', 'NDNP ',
     &         'NPAX ', 'NDNS ', 'OUTL ', 'SYMM ', 'TECON', 'TEFAC',
     &         'TENRM', 'VISN ',


c              'CT   ',     !    53: Plot deformed coordinate storage.
c              'FCIX ',     !    54: Face nodal connection list.
c              'FCZM ',     !    55: Face z-sort coordinates.
c              'FIDP ',     !    56: Faces for each node.
c              'NDER ',     !    57: Error estimator projections.
c              'NDNP ',     !    58: Stress projection array.
c              'NPAX ',     !    59: Principal axis projections.
c              'NDNS ',     !    60: Local stress projection array.
c              'OUTL ',     !    61: Outline construction.
c              'SYMM ',     !    62: Symmetric reflections table.
c              'TECON',     !    63: 3D edge definitions.
c              'TEFAC',     !    64: 3D edge factors.
c              'TENRM',     !    65: 3D edge normals
c              'VISN ',     !    66: Visible face list

c     Other data

     &         'SPTAN', 'AUR  ', 'BFGD ', 'BFGO ', 'BFGS ', 'BFGT ',
     &         'BFGV ', 'BFGW ', 'EIGE ', 'EVAL ', 'EVEC ', 'EXTND',
     &         'IPOS ', 'JPR  ', 'MO   ', 'MR   ', 'MT   ', 'MU1  ',
     &         'MU2  ', 'NDAM ', 'NMAS ', 'NSTI ', 'NREN ', 'OINMC',
     &         'OINMO', 'OINB ', 'OINC ', 'OINO ',


c              'SPTAN',     !    67: Sparse Tangent Array
c              'AUR  ',     !    68: Preconditioner Array

c              'BFGD ',     !    69: BFGS working vector
c              'BFGO ',     !    70: BFGS working vector
c              'BFGS ',     !    71: BFGS working vector
c              'BFGT ',     !    72: BFGS working vector
c              'BFGV ',     !    73: BFGS vectors, U
c              'BFGW ',     !    74: BFGS vectors, W

c              'EIGE ',     !    75: Element eigenpairs
c              'EVAL ',     !    76: Subspace eigenvalues
c              'EVEC ',     !    77: Subspace eigenvectors
c              'EXTND',     !    78: External nodes

c              'IPOS ',     !    79: Tie node list

c              'JPR  ',     !    80: Preconditioner column pointers

c              'MO   ',     !    81: Rotational DOF option
c              'MR   ',     !    82: Rotational DOF lambda
c              'MT   ',     !    83: Rotational DOF thickness
c              'MU1  ',     !    84:                             (mu1)
c              'MU2  ',     !    85:                             (mu2)

c              'NDAM ',     !    86: Nodal damping
c              'NMAS ',     !    87: Nodal mass
c              'NSTI ',     !    88: Nodal stiffness
c              'NREN ',     !    89: Nodal renumbering list

c              'OINMC',     !    90: Consistent mass equation pointers
c              'OINMO',     !    91: Consistent mass entries/equation

c              'OINB ',     !    92: Block information for out-of-core
c              'OINC ',     !    93: Sparse equation pointers
c              'OINO ',     !    94: Sparse entries/equation

c     Rigid Body data

     &         'RCG  ', 'REQRB', 'REVO ', 'RINER', 'RIRB ', 'RIXT ',
     &         'RJNT ', 'RJNX ', 'RJTU ', 'RLAMB', 'RLIST', 'RLOAD',
     &         'RMASS', 'RUROT', 'REXMS', 'REXIN',


c              'RCG  ',     !    95: Rigid body: center of mass
c              'REQRB',     !    96: Rigid body: equation numbers
c              'REVO ',     !    97: Rigid body: revolutes
c              'RINER',     !    98: Rigid body: inertia dyadic
c              'RIRB ',     !    99: Rigid body: rigid body nodes
c              'RIXT ',     !   100: Rigid body: rigid node indicators
c              'RJNT ',     !   101: Rigid body: joint data
c              'RJNX ',     !   102: Rigid body: joint coordinates
c              'RJTU ',     !   103: Rigid body: joint solution values
c              'RLAMB',     !   104: Rigid body: rotation lambda values
c              'RLIST',     !   105: Rigid body:                (nrlist)
c              'RLOAD',     !   106: Rigid body: loads
c              'RMASS',     !   107: Rigid body: mass
c              'RUROT',     !   108: Rigid body: rotational solutions

c              'REXMS',     !   109: Rigid body: Explicit mass coupling
c              'REXIN',     !   110: Rigid body: Explicit 6x6 inertias

c     Temporay arrays

     &         'TEMP1', 'TEMP2', 'TEMP3', 'TEMP4', 'TEMP5', 'TEMP6',
     &         'TEMP7', 'TEMP8', 'TEMP9', 'TEMP0',

c              'TEMP1',     !   111:  Temporary array
c              'TEMP2',     !   112:  Temporary array
c              'TEMP3',     !   113:  Temporary array
c              'TEMP4',     !   114:  Temporary array
c              'TEMP5',     !   115:  Temporary array
c              'TEMP6',     !   116:  Temporary array
c              'TEMP7',     !   117:  Temporary array
c              'TEMP8',     !   118:  Temporary array
c              'TEMP9',     !   119:  Temporary array
c              'TEMP0',     !   120:  Temporary array

c     Proportional loading table arrays (Type 2 prop. loads)

     &         'PROP0', 'PROP1', 'PROP2',

c              'PROP0',     !   121:  Prop. load offset table
c              'PROP1',     !   122:  Prop. load values table
c              'PROP2',     !   123:  Prop. load temporary use

c     Multiple support base excitations

     &         'PROBS', 'NUBAS', 'MASBS', 'PHIBS',

c              'PROBS',     !   124:  Base proportional factors
c              'NUBAS',     !   125:  Base pattern indicators
c              'MASBS',     !   126:  Static Mass projections
c              'PHIBS',     !   127:  Static modes

c     Follower nodal loads

     &         'APLOT', 'FOLLI', 'FOLLR',

c              'APLOT',     !   128:  Tag active plot elements
c              'FOLLI ',    !   129:  Follower force nodes
c              'FOLLR ',    !   130:  Follower force values

c     Contact arrays

     &         'C0   ', 'CM   ', 'ICS  ', 'HIC  ', 'CH   ',


c              'C0   ',     !   131:  Command control table      (ncc0)
c              'CM   ',     !   132:  Material table              (ncm)
c              'ICS  ',     !   133:  List of nodes for geometry (nics)
c              'HIC  ',     !   134:  History correspond vector  (nhic)
c              'CH   ',     !   135:  History values      (ch1,ch2,ch3)

c     Contact temporary arrays

     &         'CTEM1', 'CTEM2', 'CTEM3', 'CTEM4', 'CTEM5','CTEM6',
     &         'CTEM7', 'CTEM8', 'CTEM9', 'CTE10', 'CTE11','CTE12',
     &         'CTE13', 'CTE14', 'CTE15',

c              'CTEM1',     !   136:  Contact temporary
c              'CTEM2',     !   137:  Contact temporary
c              'CTEM3',     !   138:  Contact temporary
c              'CTEM4',     !   139:  Contact temporary
c              'CTEM5',     !   140:  Contact temporary
c              'CTEM6',     !   141:  Contact temporary
c              'CTEM7',     !   142:  Contact temporary
c              'CTEM8',     !   143:  Contact temporary
c              'CTEM9',     !   144:  Contact temporary
c              'CTE10',     !   145:  Contact temporary
c              'CTE11',     !   146:  Contact temporary
c              'CTE12',     !   147:  Contact temporary
c              'CTE13',     !   148:  Contact temporary
c              'CTE14',     !   149:  Contact temporary
c              'CTE15',     !   150:  Contact temporary

c     User temporary arrays

     &         'USER1', 'USER2', 'USER3', 'USER4', 'USER5', 'USER6',
     &         'USER7', 'USER8', 'USER9', 'USER0',

c              'USER1',     !   151:  Temporary array
c              'USER2',     !   152:  Temporary array
c              'USER3',     !   153:  Temporary array
c              'USER4',     !   154:  Temporary array
c              'USER5',     !   155:  Temporary array
c              'USER6',     !   156:  Temporary array
c              'USER7',     !   157:  Temporary array
c              'USER8',     !   158:  Temporary array
c              'USER9',     !   159:  Temporary array
c              'USER0',     !   160:  Temporary array

c     Blending arrays

     &         'BNODE', 'BSIDE', 'BTRAN', 'BLEND', 'BFACE', 'BNILR',

c              'BNODE',     !   161:  Super nodes for blending functions
c              'BSIDE',     !   162:  Sides for blending functions
c              'BTRAN',     !   163:  Transformations for blends
c              'BLEND',     !   164:  Blending function storage
c              'BFACE',     !   165:  Blending function storage
c              'BNILR',     !   166:  Blending layer storage

c     Rigid array

     &         'RLINK',

c              'RLINK',     !   167:  Rigid body link definitions.

c     Contact element connection array (total active)

     &         'IXC  ',

c              'IXC  ',     !   168:  Contact connection array

c     Control arrays for modal analyses

     &         'MCTRL', 'CCTRL', 'KCTRL',

c              'MCTRL',     !   169:  Mass      control array
c              'CCTRL',     !   170:  Damping   control array
c              'KCTRL',     !   171:  Stiffness control array

     &         'SVDA ', 'SVDV ', 'SVDW ', 'SVDX ',

c              'SVDA ',     !   172:  Singular valued decomp: A
c              'SVDV ',     !   173:  Singular valued decomp: V
c              'SVDW ',     !   174:  Singular valued decomp: W
c              'SVDX ',     !   175:  Singular valued decomp: X

c     Modal/Rigid node number array

     &         'IMODF', 'AFD  ', 'AFL  ', 'AFU  ', 'BFORC',

c              'IMODF',     |   176:  Modal equation numbers
c              'AFD  ',     |   177:  Modal stiffness diagonal
c              'AFL  ',     |   178:  Modal lower stiffness
c              'AFU  ',     |   179:  Modal upper stiffness
c              'BFORC',     |   180:  Modal force vector

c     Additional rigid body and modal arrays

     &         'RBEN ','RBOU ','UMOD ',

c              'RBEN',      |   181:  Rigid body number of element
c              'RBOU',      |   182:  Rigid body boundary restraints
c              'UMOD',      |   183:  Modal displacement value

c     Modal data

     &         'CTEMP', 'KTEMP', 'MTEMP', 'FSMOD', 'YYMOD', 'WBASE',

c              'CTEMP',     !   184: Modal damping parameters
c              'KTEMP',     !   185: Modal stiffness parameters
c              'MTEMP',     !   186: Modal mass parameters
c              'FSMOD',     !   187: Modal forces                (mfs)
c              'YYMOD',     !   188: Modal solution parameters   (my)
c              'WBASE',     !   189: Modal base solution parameters

c     Node type data

     &         'NDTYP',

c              'NDTYP'      !   190: Node type tags

c     Contact variables for surface descriptors

     &         'KNOTN', 'SURFP', 'INSEG', 'CNSEG', 'PNSEG', 'XISEG',

c              'KNOTN',     !   191: Node - surface listing
c              'SURFP',     !   192: Surface points
c              'INSEG',     !   193: In segment indicator
c              'CNSEG',     !   194: Contact flag
c              'PNSEG',     !   195: Point   flag
c              'XISEG',     !   196: Surface coordinates

c     Stress projection arrays

     &         'NS1  ', 'NS2  ',

c              'NS1  ',     !   197:                             (ns1)
c              'NS2  ',     !   198:                             (ns2)

c     Beam surface plot arrays

     &         'MXBEA', 'XBEAM', 'SBEAM', 'WBEAM',

c              'MXBEA',     !   199: Mesh for surface mesh of beams
c              'XBEAM',     !   200: Coordinates for surface mesh
c              'SBEAM',     !   201: Stress for surface mesh
c              'WBEAM',     !   202: Weights for surface mesh

c     Consistent damping arrays

     &         'OINDC', 'OINDO',

c              'OINDC',     !   203: Consistent damping eq pointers
c              'OINDO',     !   204: Consistent damping entry/eqn

c     Slave boundary array

     &         'NSLAV',

c              'NSLAV',     !   205: Slave node numbers

c     Normal vector

     &         'NORMV','NSCR ','VTILD','DELTX',

c              'NORMV',     !   206: Normal vector (shell use)
c              'NSCR ',     !   207: TRI2D size projection data.
c              'VTILD',     !   208: Broyden vector 1
c              'DELTX',     !   209: Broyden vector 2

c     Interface storage and Lagrange multiplier

     &         'MATCH','LAGRE','LAGRN','ULAGR','HINTE','HINT1',
     &         'HINT2','HINT3',

c              'MATCH',     !   210: Lagrange multiplier solutions
c              'LAGRE',     !   211: Lagrange multiplier equations
c              'LAGRN',     !   212: Lagrange multiplier equations
c              'ULAGR',     !   213: Lagrange multiplier solutions
c              'HINTE',     !   214: Interface history variables
c              'HINT1',     !   215: Interface history variables
c              'HINT2',     !   216: Interface history variables
c              'HINT3',     !   217: Interface history variables

c     Zienkiewicz-Zhu Projector arrays

     &         'ZZIB ','ZZIP ',

c              'ZZIB ',     !   218: Zienkiewicz-Zhu boundary nodes
c              'ZZIP ',     !   219: Zienkiewicz-Zhu active nodes

c     Auto contact pointers

     &         'ACON2','ASLD2','ACIQ2','ACON3',

c              'ACON2',     !   220: Autocon array: length = numnp
c              'ASLD2',     !   221: Autocon slideline: lg = 2*numnp+8
c              'ACIQ2',     !   222: Autocon IP   :     lg = ip(numnp)
c              'ACON3',     !   223: Autocon array :    lg = ip(numnp)*2

c     Contact lagrange multipliers

     &         'IAD  ',

c              'IAD  ',     !   224: Contact lagrange multiplier array

c     User solver pointers

     &         'USOL1','USOL2','USOL3','USOL4','USOL5','USOL6',
     &         'USOL7','USOL8','USOL9','USOL0',

c              'USOL1',     !   225: User solver array
c              'USOL2',     !   226: User solver array
c              'USOL3',     !   227: User solver array
c              'USOL4',     !   228: User solver array
c              'USOL5',     !   229: User solver array
c              'USOL6',     !   230: User solver array
c              'USOL7',     !   231: User solver array
c              'USOL8',     !   232: User solver array
c              'USOL9',     !   233: User solver array
c              'USOL0',     !   234: User solver array

c     Diagonal scaling array

     &         'DSCA1','DSCA2','DSCA3','DSCA4',

c              'DSCA1',     !   235: Reciprocal to diagonal sqare root
c              'DSCA2',     !   236: Reciprocal to diagonal sqare root
c              'DSCA3',     !   237: Reciprocal to diagonal sqare root
c              'DSCA4',     !   238: Reciprocal to diagonal sqare root

c     Interface type array

     &         'ITYPE','IEDOF',

c              'ITYPE',     !   239: Interface types
c              'IEDOF',     !   240: Interface types

c     Surface load real data

     &         'SLODR','EULER','LEULR',

c              'SLODR',     !   241: Surface load real data
c              'EULER',     !   242: Euler angles
c              'LEULR',     !   243: Euler angles, element

c     Parallel solver arrays

     &         'GN   ','EQ   ','DNNZ ','ONNZ ','GETP ',
     &         'GETV ','SENP ','SENV ',

c              'GN   ',     !   244: Local to Global node mapping
c              'EQ   ',     !   245: Local node to Global eq mapping
c              'DNNZ ',     !   246: Number non-zero diag-block entries
c              'ONNZ ',     !   247: No. non-zero off-diag-block entries
c              'GETP ',     !   248: Get data partition pointers
c              'GETV ',     !   249: Get data node values
c              'SENP ',     !   250: Send data partition pointers
c              'SENV ',     !   251: Send data node values

c     Mesh partioner (METIS/PARMETIS) arrays

     &         'XADJN','NODG ','NODPT','VXDST',

c              'XADJN',     !   252: Pointer for nodal adjacency list
c              'NODG ',     !   253: Nodal adjacency list
c              'NODPT',     !   254: Nodal partition numbers
c              'VXDST',     !   255: Distribution array

c     Nodal follower forces

     &         'NFORC','ELINK','GUVAL',

c              'NFORC',     !   256: Nodal follower forces
c              'ELINK',     !   257: Edge link direction list
c              'GUVAL',     !   258: Global equation values

c     Representative volume element for multi-scale analysis

     &         'RVELM','FRVEL','SRVEL',

c              'RVELM'      !   259: Representative volume elements
c              'FRVEL'      !   260: Deformation gradient for RVE
c              'SRVEL'      !   261: Stress and tangent from RVE

     &         'IDSAV',

c              'IDSAV'      !   262: Eq numbers without multipliers

c     NURB coordinate storage

     &         'NURBS','LNURB',

c              'NURBS'      !   263: Nurb weights
c              'LNURB'      !   264: Nurb local weights

c     Load table storage

     &         'LDTAB','LDNOD','LDVAL','SPINS',

c              'LDTAB'      !   265: Load table pointers
c              'LDNOD'      !   266: Load node list
c              'LDVAL'      !   267: Load node values
c              'SPINS'      !   268: Spin node values

c     Representative Volume Element storage

     &         'RVEMA','RVESD' ,

c              'RVEMA'      !   269: List of material numbers for RVE
c              'RVESD'      !   270: List of material numbers for RVE

c     E-spin data

     &         'ESPIN','ESPTR','CEPTR','TRIAD','LTRIA',

c              'ESPIN'      !   271: List of nodes to spin
c              'ESPTR'      !   272: Pointer to start of spin node lists
c              'CEPTR'      !   273: Pointer array ofr extraction matrix

c     3-d Boundary rotation triad

c              'TRIAD'      !   274: 3-d boundary triads global
c              'LTRIA'      !   275: 3-d boundary triads element

c     Plot arrays for NURBS/T-spline projections

     &         'I_LIN','X_LIN','U_LIN','N_LIN','T_LIN','S_LIN',
     &         'P_LIN','E_LIN','ID_LN','LKNOT','NELEM',

c              'I_LIN',     !   276: NURB 3-d graphics elements
c              'X_LIN',     !   277: NURB 3-d graphics coordinatates
c              'U_LIN',     !   278: NURB 3-d graphics displacements
c              'N_LIN',     !   279: NURB 3-d graphics renumber lists
c              'T_LIN',     !   280: NURB 3-d graphics node type
c              'S_LIN',     !   281: NURB 3-d graphics stresses
c              'P_LIN',     !   282: NURB 3-d graphics stresses
c              'E_LIN',     !   283: NURB 3-d graphics eigenvectors
c              'ID_LN',     !   284: NURB 3-d graphics id-array
c              'LKNOT',     !   285: NURB local knot definitions
c              'NELEM',     !   286: NURB submesh nodes

c     T-Spline storage for extraction operators and Bezier extraction

     &         'KNOTP','KNOTV','C_E  ','RC_E ','IC_E ','P_BEZ',
     &         'X_BEZ','W_BEZ','IXBEZ',

c              'KNOTP'      !   287: NURB knot pointer
c              'KNOTV'      !   288: NURB knot values
c              'C_E  '      !   289: NURB Extraction operator
c              'RC_E '      !   290: NURB Extraction operator
c              'IC_E '      !   291: NURB Extraction pointer
c              'P_BEZ'      !   292: NURB Bezier orders
c              'X_BEZ'      !   293: NURB Bezier nodes
c              'W_BEZ'      !   294: NURB Bezier nodes
c              'IXBEZ'      !   295: NURB Bezier elements

c     PLot arrays for NURBS/T-spline transient

     &         'V_LIN','EC_E ' ,

c              'V_LIN',     !   296: NURB 3-d graphics displacements
c              'EC_E ',     !   297: NURB 3-d graphics displacements

c     Set memory for nurb nodes, sides, and blends

     &         'KNOTS','NSIDE','KTDIV','NBSID','LKSID',

c              'KNOTS'      !   298: List of knot vectors for NURBS
c              'NSIDE'      !   299: List of side vectors for NURBS
c              'KTDIV'      !   300: Knot Div for NURB mesh builds
c              'NBSID'      !   301: NURB Block side numbers
c              'LKSID'      !   302: NURB Block length and knot numbers

c     Plot arrays for NURBS/T-spline history variables

     &         'HPLTB','HSELM','HDNP ','H_LIN',

c              'HPLTB'      !   303: History plot table
c              'HSELM'      !   304: Element array for history plots
c              'HDNP '      !   305: Global  array for history plots
c              'H_LIN'      !   306: Global  array for history plots

c     Sparse solver integer array

     &         'IPTAN',

c              'IPTAN',     !   307: Sparse integer tangent array

c     NURBS integer array size arrays

     &         'KNOTL', 'SIDEL', 'BLOKL', 'KTNUM', 'NBLOK',
     &         'NESID',

c              'KNOTL',     !   308:  KNOT (lknot)
c              'SIDEL',     !   309:  SIDE (lside,kside)
c              'BLOKL',     !   310:  BLOCK (nblkdm,nurmat,nuregn->nblk)
c              'KTNUM',     !   311:  KNOT (knotnum,knotlen)
c              'NBLOK',     !   312:  NURB (nblksd)
c              'NESID',     !   313:  NURB (nepatch)

c     NURBS Unused

     &         'UNUR1', 'UNUR2', 'UNUR3', 'UNUR4', 'UNUR5',

c              'UNUR1',     !   314:  NURB Unused
c              'UNUR2',     !   315:  NURB Unused
c              'UNUR3',     !   316:  NURB Unused
c              'UNUR4',     !   317:  NURB Unused
c              'UNUR5',     !   318:  NURB Unused

c     Diagonal Stiffness for nodal stiffness/damper/mass

     &         'KNOTE', 'DTANG',

c              'KNOTE',     !   319:  NURB element size
c              'DTANG',     !   320:  Diagonal for nodal LHS

c     Plot temporary arrays

     &         'PTEM1', 'PTEM2', 'PTEM3', 'PTEM4', 'PTEM5',
     &         'PTEM6', 'PTEM7', 'PTEM8', 'PTEM9', 'PTEM0',

c              'PTEM1',     !   321: Plot temporary variable
c              'PTEM2',     !   322: Plot temporary variable
c              'PTEM3',     !   323: Plot temporary variable
c              'PTEM4',     !   324: Plot temporary variable
c              'PTEM5',     !   325: Plot temporary variable
c              'PTEM6',     !   326: Plot temporary variable
c              'PTEM7',     !   327: Plot temporary variable
c              'PTEM8',     !   328: Plot temporary variable
c              'PTEM9',     !   329: Plot temporary variable
c              'PTEM0',     !   330: Plot temporary variable

c     Hill-Mandel arrays

     &         'HILLI', 'HILLG', 'HILLX',

c              'HILLI',     !   331: Hill-Mandel IXL array
c              'HILLG',     !   332: Hill-Mandel   G array
c              'HILLX',     !   333: Hill-Mandel  XS array

c     Mixed model extraction arrays

     &         'CPTMX', 'C_EMX'/

c              'CPTMX',     !   334: Mixed model extraction pointer
c              'C_EMX',     !   335: Mixed model extraction operator

      if(num.eq.0) then

c       Zero pointer arrays

        do i = 1,num_nps
          np(i) = 0
        end do ! i

        do i = 1,num_ups
          up(i) = 0
        end do ! i

        do i = list+1,600
          nlist(i) = '     '
        end do ! i

        do i = 1,600
          ilist(1,i) = 0
          ilist(2,i) = 0
        end do ! i
        llist  =  num_nps

      endif

c     Check user allocations then do allocation operation

      palloc = ualloc(num-llist,name,length,precis)

      end
