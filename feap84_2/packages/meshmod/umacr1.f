c$Id:$
      subroutine umacr1(lct,ctl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c     Author: R.L. Taylor                                      8/14/2006
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Convert each 8-node quadrilateral into a 10-node quad.

c      Inputs:
c         lct       - Command character parameters
c         ctl(3)    - Command numerical parameters

c      Outputs:
c         File:     - 'Coord_10' : file containing nodal coordinates
c                     for added node 10 values.
c                   - 'Elemt_10' : file containing nodal connections
c                     for added node 10 of elements that originally had
c                     8 nodes.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'cdat1.h'
      include  'iofile.h'
      include  'umac1.h'
      include  'sdata.h'

      include  'pointer.h'
      include  'comblk.h'

      logical   pcomp
      integer   nn
      character lct*15
      real*8    ctl(3)

      save

c     Set command word

      if(pcomp(uct,'mac1',4)) then      ! Usual    form
        uct = 'ten8'                    ! Specify 'ten8' convert
      elseif(urest.eq.1) then           ! Read  restart data

      elseif(urest.eq.2) then           ! Write restart data

c     Set coordinates and element array for node 11

      else

c       Check for available space in ix array

        if(nen.lt.8) then
          write(iow,2001) ' *ERROR*  Need at least 8-nodes for',
     &                    ' maximum number of nodes on each element.'
          write(ilg,2001) ' *ERROR*  Need at least 8-nodes for',
     &                    ' maximum number of nodes on each element.'
          if(ior.lt.0) then
            write(*,2001) ' *ERROR*  Need at least 8-nodes for',
     &                    ' maximum number of nodes on each element.'
          endif
        else
          write(iow,2002)
          if(ior.lt.0) then
            write(*,2002)
          endif
          nn = numnp
          call ckixsz8(mr(np(32)),mr(np(33)),nie,nen,nen1,numel, nn)
          nn = numnp
          call ckset8(mr(np(32)),mr(np(33)),nie,nen,nen1,numel,ndm,
     &                hr(np(43)), nn)
        endif

      end if

c     I/O formats

2001  format(a)

2002  format(' ---> Convert 8-node quads to 10-node ones'/
     &       '      Nodal   coordinates in file = Coor_10'/
     &       '      Element connections in file = Elmt_10'/)

      end

      subroutine ckixsz8(ie,ix,nie,nen,nen1,numel, nn)

c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Add node 9,10 to 8-node quadrilaterals

c      Inputs:
c         ie(nie,*)   - Element property array
c         ix(nen1,*)  - Element connection array
c         nie         - Dimension for ie array
c         nen         - Maximum number of nodes/element
c         nen1        - Dimension for ix array
c         numel       - Number of elements

c      Outputs:
c         nn          - Number of node: counter for adding nodes 9 & 10
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'iodata.h'

      character  eformat*27, size*2
      integer    nie,nen,nen1,numel, nn, n,i,j
      integer    ie(nie,*),ix(nen1,*)

c     Open file for saving element connections

      open(unit=ios,file='Elmt_10',status='unknown',form='formatted')

c     Build output format for elements to fit data

c     Format:    1...5...10...15...20...25...30
      eformat = '(i10,i6,i10,13i10:/(16i10))'

      n = nint(log10(dble(nn+numel))) + 2
      write(size,'(i2)') n
      eformat( 3: 4) = size
      eformat(10:11) = size
      eformat(16:17) = size
      eformat(24:25) = size

c     Compute and output element list

      write(ios,'(a)') 'ELEMents'

c     Compute element list

      do n = 1,numel

c       Check for a surface element

        if(ie(1,ix(nen1,n)).eq.2) then
          do i = nen,1,-1
            if(ix(i,n).gt.0) then
              exit
            endif
          end do ! i
          if(i.eq.8) then
            nn = nn + 1
            write(ios,eformat) n,0,ix(nen1,n),(ix(j,n),j=1,8),nn,nn+1
            nn = nn + 1
          else
            write(ios,eformat) n,0,ix(nen1,n),(ix(j,n),j=1,i)
          endif
        endif
      end do ! n

c     Close file

      write(ios,'(a)') ' '
      close(ios,status='keep')

      end

      subroutine ckset8(ie,ix,nie,nen,nen1,numel, ndm,x, nn)

c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Set coordinate values for nodes 9 & 10 on quadrilateral
c                elements

c      Inputs:
c         ie(nie,*)   - Element property array
c         ix(nen1,*)  - Element connection array
c         nie         - Dimension for ie array
c         nen         - Maximum number of nodes/element
c         nen1        - Dimension for ix array
c         numel       - Number of elements
c         ndm         - Spatial dimension of nodal coordinate array
c         x(ndm,*)    - Nodal coordinate array existing

c      Outputs:
c         nn          - Number of node: counter for adding nodes 9 & 10
c         x(ndm,*)    - Nodal coordinate array augmented by internal
c                       value at nodes 9 & 10 of a quadrilateral
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'iodata.h'

      character  nformat*19, size*2
      integer    nie,nen,nen1,numel,ndm, nn, n,i
      integer    ie(nie,*),ix(nen1,*)
      real*8     x(ndm,*), xx(ndm)

c     Open file for saving coordinate values

      open(unit=ios,file='Coor_10',status='unknown',form='formatted')

c     Build output format for coordinates to fit data

c     Format:    1...5...10...15...20
      nformat = '(i10,i6,1p,3e16.8)'

      n = nint(log10(dble(nn+numel))) + 2
      write(size,'(i2)') n
      nformat( 3: 4) = size

      write(ios,'(a)') 'COORdinates all'
      do n = 1,nn
        write(ios,nformat) n,0,(x(i,n),i=1,ndm)
      end do ! n

      do n = 1,numel

c       Check for a quadrilateral element

        if(ie(1,ix(nen1,n)).eq.2) then
          do i = nen,1,-1
            if(ix(i,n).gt.0) then
              exit
            endif
          end do ! i
          if(i.eq.8) then
            nn = nn + 1
            do i = 1,ndm
              xx(i) = (4.0d0*(x(i,ix( 4,n))
     &                      + x(i,ix( 5,n))
     &                      + x(i,ix( 6,n)))
     &                      - x(i,ix( 1,n))
     &                      - x(i,ix( 2,n))
     &                      - x(i,ix( 3,n)))/9.0d0
            end do ! i
            write(ios,nformat) nn,0,(xx(i),i=1,ndm)
            nn = nn + 1
            write(ios,nformat) nn,0,(xx(i),i=1,ndm)
          endif
        endif
      end do ! n

c     Close file

      write(ios,'(a)') ' '
      close(ios,status='keep')

      end
