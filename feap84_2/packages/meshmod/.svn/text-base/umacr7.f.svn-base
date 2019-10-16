c$Id:$
      subroutine umacr7(lct,ctl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c     Author: R.L. Taylor                                      8/14/2006
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Convert each 10-node tetrahedron to an 11-node
c                tetrahedron.

c      Inputs:
c         lct       - Command character parameters
c         ctl(3)    - Command numerical parameters

c      Outputs:
c         File:     - 'coord_11'   : file containing nodal coordinates
c                     for added node 11 values.
c                   - 'element_11' : file containinc nodal connections
c                     for added node 11 of elements that originally had
c                     10 nodes.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'cdat1.h'
      include  'iofile.h'
      include  'umac1.h'
      include  'sdata.h'

      include  'pointer.h'
      include  'comblk.h'

      logical   pcomp, setval,palloc
      integer   nn
      character lct*15
      real*8    ctl(3)

      save

c     Set command word

      if(pcomp(uct,'mac7',4)) then      ! Usual    form
        uct = 'ten1'                    ! Specify 'ten1'1 convert
      elseif(urest.eq.1) then           ! Read  restart data

      elseif(urest.eq.2) then           ! Write restart data

c     Set coordinates and element array for node 11

      else

c       Check for available space in ix array

        if(nen.lt.11) then
          write(iow,2001) ' *ERROR*  Need at least 11-nodes for',
     &                    ' maximum number of nodes on each element.'
          write(ilg,2001) ' *ERROR*  Need at least 11-nodes for',
     &                    ' maximum number of nodes on each element.'
        else
          write(iow,2002)
          if(ior.lt.0) then
            write(*,2002)
          endif
          nn = numnp
          call ckixsz(mr(np(32)),mr(np(33)),nie,nen,nen1,numel, nn)
          setval = palloc(43,'X    ',nn*ndm,2)
          nn = numnp
          call cksetx(mr(np(32)),mr(np(33)),nie,nen,nen1,numel,ndm,
     &                hr(np(43)), nn)
        endif

      end if

c     I/O formats

2001  format(a)

2002  format(' ---> Convert 10-node tets to 11-node tets'/
     &       '      Nodal   coordinates in file = coord_11'/
     &       '      Element connections in file = element_11'/)

      end

      subroutine ckixsz(ie,ix,nie,nen,nen1,numel, nn)

c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Add node 11 to 10-node tets

c      Inputs:
c         ie(nie,*)   - Element property array
c         ix(nen1,*)  - Element connection array
c         nie         - Dimension for ie array
c         nen         - Maximum number of nodes/element
c         nen1        - Dimension for ix array
c         numel       - Number of elements

c      Outputs:
c         nn          - Number of node: counter for adding node 11
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'iodata.h'

      character  eformat*27, size*2
      integer    nie,nen,nen1,numel, nn, n,i
      integer    ie(nie,*),ix(nen1,*)

c     Compute element list

      do n = 1,numel

c       Check for a solid element

        if(ie(1,ix(nen1,n)).eq.3) then
          do i = nen,9,-1
            if(ix(i,n).gt.0) then
              exit
            endif
          end do ! i
          if(i.eq.10) then
            nn       = nn + 1
            ix(11,n) = nn
          endif
        endif
      end do ! n

c     Open file for saving element connections

      open(unit=ios,file='element_11',status='unknown',form='formatted')

c     Build output format for elements to fit data

c     Format:    1...5...10...15...20...25...30
      eformat = '(i10,i6,i10,13i10:/(16i10))'

      n = nint(log10(dble(max(nn,numel)))) + 2
      write(size,'(i2)') n
      eformat( 3: 4) = size
      eformat(10:11) = size
      eformat(16:17) = size
      eformat(24:25) = size

c     Compute and output element list

      write(ios,'(a)') 'ELEMents all'
      do n = 1,numel
        write(ios,eformat) n,0,ix(nen1,n),(ix(i,n),i=1,nen)
      end do ! n
      write(ios,'(a)') ' '

c     Close file

      close(ios,status='keep')

      end

      subroutine cksetx(ie,ix,nie,nen,nen1,numel, ndm,x, nn)

c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Set coordinate values for node 11 on tetrahedral
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
c         nn          - Number of node: counter for adding node 11
c         x(ndm,*)    - Nodal coordinate array augmented by internal
c                       value at node 11 of tetrahedron
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'iodata.h'

      character  nformat*19, size*2
      integer    nie,nen,nen1,numel,ndm, nn, n,i
      integer    ie(nie,*),ix(nen1,*)
      real*8     x(ndm,*)

      do n = 1,numel

c       Check for a solid element

        if(ie(1,ix(nen1,n)).eq.3) then
          do i = nen,10,-1
            if(ix(i,n).gt.0) then
              exit
            endif
          end do ! i
          if(i.eq.11) then
            nn       = nn + 1
            do i = 1,ndm
              x(i,nn) = 0.250d0*(x(i,ix( 5,n))
     &                         + x(i,ix( 6,n))
     &                         + x(i,ix( 7,n))
     &                         + x(i,ix( 8,n))
     &                         + x(i,ix( 9,n))
     &                         + x(i,ix(10,n)))
     &                - 0.125d0*(x(i,ix( 1,n))
     &                         + x(i,ix( 2,n))
     &                         + x(i,ix( 3,n))
     &                         + x(i,ix( 4,n)))
            end do ! i
          endif
        endif
      end do ! n

c     Open file for saving coordinate values

      open(unit=ios,file='coord_11',status='unknown',form='formatted')

c     Build output format for coordinates to fit data

c     Format:    1...5...10...15...20
      nformat = '(i10,i6,1p,3e16.8)'

      n = nint(log10(dble(nn))) + 2
      write(size,'(i2)') n
      nformat( 3: 4) = size

      write(ios,'(a)') 'COORdinates all'
      do n = 1,nn
        write(ios,nformat) n,0,(x(i,n),i=1,ndm)
      end do ! n
      write(ios,'(a)') ' '

c     Close file

      close(ios,status='keep')

      end
