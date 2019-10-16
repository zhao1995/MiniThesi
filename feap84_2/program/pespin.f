c$Id:$
      subroutine pespin(x,isp,ndtyp,ndm,numnp,nspin,prt,prth)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    03/01/2011
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set edge spin conditions.

c      Inputs:
c         x(ndm,*)   - Nodal coordinates of mesh
c         ndtyp(*)   - Node type (>=0 exist; <0 tied)
c         ndm        - Spatial dimension of mesh
c         numnp      - Number of nodes in mesh
c         prt        - Output generated results if true
c         prth       - Output title/header data if true

c      Outputs:
c         isp(*)     - List of spin nodes for edge 'i'
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'
      include  'pointer.h'
      include  'comblk.h'

      logical   prt,prth
      integer   ndm,nsp,numnp,nspin, i,j,n, ndtyp(*), isp(numnp)
      real*8    pdiff, dx,x0,xx, gap, x(ndm,numnp),xc(2)

      save

c     Read input of boundary edge for restraints

      gap   = 1.0d-3/sqrt(dble(max(1,numnp)))
      nsp   = 0
      do j = 1,ndm
        xc(j) = 0.0d0
      end do ! j
      i  = 3
      x0 = 5.0d0
      dx = pdiff(x,i,ndm,numnp)*gap
      do n = 1,numnp
        xx = abs(x(i,n)-x0)
        if(ndtyp(n).ge.0 .and. xx.le.dx) then
          nsp      = nsp + 1
          isp(nsp) = n
        endif
      end do ! n
c     go to 100

      call prtitl(prth)
      if(prt) then
        write(iow,2000)
        if(ior.lt.0) then
          write(*,2000)
        endif
      endif

      do n = 1,nsp
        if(prt .and. ndtyp(n).ge.0) then
          write(iow,2001) n,isp(n)
          if(ior.lt.0) then
            write(*,2001) n,isp(n)
          endif
        endif
      end do ! n

      nspin = nsp

      call pldtabl(mr(np(265)),1,1,1,4)
      call pldtabl(mr(np(265)),1,1,2,4)
      call pldtabl(mr(np(265)),1,nspin,2,2)

c     Formats

2000  format('  E d g e    S p i n    N o d e s'/4x,'Spin',4x,'Node')
2001  format(2i8)

      end
