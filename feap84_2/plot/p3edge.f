c$Id:$
      subroutine p3edge(iface,x,ndm,nface,ifc1)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute data to identify edges for 3-d objects

c      Inputs:
c         iface(ifc1,*) - Face node numbers, matl sets, regions, etc.
c         x(ndm,*)   - Nodal coordinates for mesh
c         ndm        - Dimension of x array
c         nface      - Number of faces

c      Outputs:
c         none       - Through pointers to arrays
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'comblk.h'
      include  'pointer.h'

      logical   setvar,palloc
      integer   ndm,nface,ifc1
      integer   iface(ifc1,nface)
      real*8    x(ndm,numnp)

      save

c     Allocate memory for computing faces and normals

      setvar = palloc( 64,'TEFAC',numnp+1,  1)
      setvar = palloc( 65,'TENRM',max(nface,numnp)*3,2)

      call p3elfa(x,mr(np(64)),hr(np(65)),iface,nface,ndm,ifc1,numnp)

c     Allocate memory for computing elements surrounding a face

      setvar = palloc( 63,'TECON',mr(np(64)+numnp), 1)

      call p3econ(iface,mr(np(64)),mr(np(63)),nface,ifc1,numnp)

      end
