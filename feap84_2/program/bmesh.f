c$Id:$
      subroutine bmesh(ie,d,id,x,ix,f,t,an,
     &                 ndd,nie,ndm,ndf,nen,nen1,numnp,numel,nummat)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:        Read mesh data from binary mode data file.

c      Inputs:
c         ndd        - Number of entries in d array/material set
c         nie        - Dimension of ie array
c         ndm        - Dimension of x array
c         ndf        - Dimension of id and f arrays
c         nen        - Number of nodes/element
c         nen1       - Dimension of ix array
c         numnp      - Number of nodes
c         numel      - Number of elements
c         nummat     - Number of material sets

c      Outputs:
c         ie(nie,*)  - Assembly information for material sets
c         d(*)       - Material set parameters
c         id(ndf,*)  - Nodal equation numbers/b.c. codes
c         x(ndm,*)   - Nodal coordinates
c         ix(nen1,*) - Element nodal connection data
c         f(ndf,*,2) - Nodal loads/displacements
c         t(*)       - Nodal temperatures
c         an(*)      - Nodal angles for sloping b.c.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iodata.h'
      include  'mdata.h'
      include  'pdata6.h'

      integer   ndd,nie,ndm,ndf,nen,nen1,numnp,numel,nummat, i,nn
      integer   ie(nie,nummat),id(ndf,numnp),ix(nen1,numel)
      real*8    d(ndd,nummat),x(ndm,numnp),f(ndf,numnp,2),t(numnp)
      real*8    an(numnp)

      save

c     Read binary mode data from file

      read(ios) ((x(i,nn),i=1,ndm),nn=1,numnp)
      read(ios) ((ix(i,nn),i=1,nen),ix(nen1,nn),nn=1,numel)
      read(ios) ((id(i,nn),i=1,ndf),nn=1,numnp)
      read(ios) (( f(i,nn,1),i=1,ndf),nn=1,numnp)
      read(ios) (( f(i,nn,2),i=1,ndf),nn=1,numnp)
      read(ios) (  t(nn)  ,nn=1,numnp)
      read(ios) ( an(nn)  ,nn=1,numnp)
      read(ios) ((ie(i,nn),i=1,nie),nn=1,nummat)
      read(ios) (( d(i,nn),i=1,ndd),nn=1,nummat)
      read(ios) ia,inord,ipord

c     Reset B.C. codes

      do nn = 1,numnp
        do i = 1,ndf
          if(id(i,nn).gt.0) then
            id(i,nn) = 0
          else
            id(i,nn) = 1
          endif
        end do ! i
      end do ! nn

      end
