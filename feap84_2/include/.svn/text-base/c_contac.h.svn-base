
                             ! MAIN CONTACT COMMON
c     Modfifications:                                   Date: dd/mm/yyyy
c     1.  Add parameter 'ifpck' for penetration check         03/05/2008

      integer  cck,numcs,numcm,numcp,ofssurf,ofsmate,isgp1,isgp3,indb
      logical  ifct,ifdb,ifistgn,ifchist,iffron,lagrm,ifpck
      common /c_contac/
     &  cck(c_ncc),          ! Contact Commands Kounter (# read)
     &  numcs,               ! # Contact Surfaces  (= cck(1))
     &  numcm,               ! # Contact Materials (= cck(2))
     &  numcp,               ! # Contact Pairs     (= cck(3))
     &  ofssurf,             ! counting offset for surface data
     &  ofsmate,             ! counting offset for material data
     &  isgp1,               ! CH1 & CH2 counting offset
     &  isgp3,               ! CH3       counting offset
     &  indb,                ! Debug level
     &  ifct,                ! ConTact Flag, if true contact is active
                             ! set by CONT,ON - CONT,OFF
     &  ifdb,                ! # DeBug Flag, if true debug is active
                             ! set by CONT,DEBU, output units 98, 99
     &  ifistgn,             ! Geometrical check flag
     &  ifchist,             ! CHange of geometry status flag
     &  iffron,              ! FRiction ON flag
     &  lagrm,               ! LAGRange Multiplier flag
     &  ifpck                ! Penetration ChecK flag
