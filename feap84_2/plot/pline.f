c$Id:$
      subroutine pline(x,ie,ix,id,ip,numnp,numel,ndm,
     &                 nen1,nen,nie,ct,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Pass ct(1) to xpline                             01/05/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Mesh plot and outline driver

c      Inputs:
c         x(ndm,*)  - Nodal coordinates
c         ie(nie,*) - Assembly data for material sets
c         ix(nen1,*)- Element nodal connection list
c         id(*)     - Number of elements attached to nodes
c         ip(8,*)   - Sorted element order for each quadrant
c         numnp     - Number of nodes in mesh
c         numel     - Number of elements in mesh
c         ndm       - Dimesion of x array
c         nen1      - Dimesion of ix array
c         nie       - Dimesion of ie array
c         ct(*)     - Color changing
c         isw       - Flag, Plot mesh if true, otherwise outline

c      Outputs:
c         none      - Plot output to screen/file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pointer.h'
      include  'comblk.h'

      logical   isw, setvar,palloc
      integer   numnp,numel,ndm,nen1,nen,nie
      integer   ie(*),ix(*),id(*),ip(*)
      real*8    x(*),ct(*)

      save

c     Determine necessary storage for mesh lines and allocate storage

      call xcompp(ie,ix,id,nie,nen1,nen,numnp,numel)

      setvar = palloc(112,'TEMP2',id(numnp+1),1)

c     Draw mesh

      call xpline(x,ie,ix,id,mr(np(112)),ip,numnp,numel,ndm,
     &            nen1,nen,nie,ct(1),isw)

c     Delete storage for mesh lines

      setvar = palloc(112,'TEMP2',0,1)

      end
