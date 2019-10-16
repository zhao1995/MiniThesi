c$Id:$
      subroutine uplot9(ctl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c      1. Change call palloc to setvar = palloc             14/08/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Plot interface to PARAVIEW

c      Inputs:
c         ctl(3)    - Command numerical parameters

c      Outputs:
c         N.B.  Users are responsible for command actions.  See
c               programmers manual for example.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'evdata.h'
      include  'iofile.h'
      include  'pointer.h'
      include  'sdata.h'
      include  'strnum.h'
      include  'umac1.h'

      include  'comblk.h'

      logical   pcomp, setvar, palloc
      character lct*15
      integer   i,ii,jj, node, plu,ix(32)
      real*8    ctl(3)

      save

      data      plu / 99 /

c     Set command word

      if(pcomp(uct,'plt9',4)) then      ! Usual    form
        uct = 'pvie'                    ! Specify 'name' ''

      else                              ! Perform user operation

        lct = 'PARAFEAP_00.vtu'
        node = nint(ctl(1))
        if(node.ge.1 .and. node.lt.10) then
          write(lct(11:11),'(i1)') node
        elseif(node.ge.10 .and. node.lt.100) then
          write(lct(10:11),'(i2)') node
        endif
        write(*,*) 'Saving PARAVIEW data to ',lct
        open(unit=plu,file=lct,access='sequential')

c       Write out top header
        write(plu,1000)

        write(plu,1020) numnp,numel               ! Start Mesh/Piece Section
        write(plu,1010) '<Points>'                ! Start Point/Node data

        write(plu,1030) 'Float64','nodes',3       ! ndm
        do i = 1,numnp
          write(plu,2000) (hr(npxx+(i-1)*ndm+(ii-1)),ii = 1,ndm)
     &                   ,(0.0d0,ii = ndm+1,3)
        end do ! i
        write(plu,1010) '</DataArray>'            ! Close Node data

        write(plu,1010) '</Points>'               ! Close Points section

        write(plu,1010) '<Cells>'                 ! Start Cell Section
        write(plu,1030) 'Int32','connectivity',1  ! Start Elements

c       Offsets memory allocation

        setvar = palloc(111,'TEMP1',numel+1,1)

        mr(np(111)) = 0;
        do i = 1,numel
          jj            = 0
          do ii = 1,nen
            node = mr(np(33) + ii-1 + nen1*(i-1))
            if (node .ne. 0) then
              jj     = jj + 1
              ix(jj) = node - 1
            endif
          end do ! ii
          mr(np(111)+i) = mr(np(111)+i-1) + jj
          write(plu,2010) (ix(ii),ii=1,jj)
        end do ! i

        write(plu,1010) '</DataArray>'           ! Close Elements

        write(plu,1030) 'Int32','offsets',1      ! Start Offsets
        write(plu,2010) (mr(np(111)+i), i = 1,numel)

        write(plu,1010) '</DataArray>'            ! Close Offsets

        write(plu,1030) 'UInt8','types',1                  ! Start Element types
        do i = 1,numel,10
          do ii = 1,10
            if (mr(np(111)+i)-mr(np(111)+i-1) .eq. 2) then     ! 2 node line
              ix(ii) = 3
            elseif (mr(np(111)+i)-mr(np(111)+i-1) .eq. 3) then ! 3 node triangle
              ix(ii) = 5
            elseif (mr(np(111)+i)-mr(np(111)+i-1) .eq. 4) then ! 4 node quad
              ix(ii) = 9
            elseif (mr(np(111)+i)-mr(np(111)+i-1) .eq. 8) then ! 8 node brick
              ix(ii) = 12
            endif
          end do ! ii
          write(plu,2010) (ix(jj),jj=1,min(10,numel-i+1))
        end do ! i
        setvar = palloc(111,'TEMP1',0,1)

        write(plu,1010) '</DataArray>'              ! Close Element types
        write(plu,1010) '</Cells>'                  ! Close Cell Section

        write(plu,1010) '<PointData>'               ! Start Point Data

        write(plu,1030) 'Float64','Displ',ndf        ! Start Displacements
        do i = 1,numnp
          write(plu,2000) (hr(npuu+(i-1)*ndf+ii),ii = 0,ndf-1)
        end do ! i
        write(plu,1010) '</DataArray>'               ! Close Displacements

        if(np(42).ne.0) then
          write(plu,1030) 'Float64','Velocity',ndf   ! Start Velocity
          do i = 1,numnp
            write(plu,2000) (hr(npud+(i-1)*ndf+ii),ii = 0,ndf-1)
          end do ! i
          write(plu,1010) '</DataArray>'             ! Close Velocity

          write(plu,1030) 'Float64','Acceleration',ndf ! Start Acceleration
          do i = 1,numnp
            write(plu,2000) (hr(npud+nneq+(i-1)*ndf+ii),ii = 0,ndf-1)
          end do ! i
          write(plu,1010) '</DataArray>'              ! Close Acceleration
        endif

        if(abs(istv).gt.0) then
          write(plu,1030) 'Float64','Stress',abs(istv)   ! Start Stresses
          do i = 1,numnp
            write(plu,2000) (hr(npnp+(i-1) + ii*numnp),ii=1,abs(istv))
          end do
          write(plu,1010) '</DataArray>'              ! Close Stresses

          write(plu,1030) 'Float64','PStress',7       ! Start Principal Stresses
          do i = 1,numnp
            write(plu,2000) (hr(nper+(i-1) + ii*numnp),ii=1,7)
          end do ! i
          write(plu,1010) '</DataArray>'              ! Close Stresses

        else
          write(*,*) ' No stresses output to Paraview file'
        endif

        write(plu,1010) '</PointData>'              ! Close Point Data Section

        write(plu,1010) '</Piece>'                  ! Close Mesh/Piece

c       Close the XML file

        write(plu,1010) '</UnstructuredGrid> </VTKFile>'
        close(plu, status = 'keep')

      endif

c     I/O Formats

1000  format('<?xml version="1.0"?>',/
     &       '<VTKFile type="UnstructuredGrid" version="0.1">',/
     &       '<UnstructuredGrid>')

1010  format(a)

1020  format('<Piece NumberOfPoints="',i10,
     &       '" NumberOfCells="',i10,'">')

1030  format('<DataArray type="',a,'" Name="',a,
     &       '" NumberOfComponents="',i2,'" format="ascii">')

2000  format(1p,6e13.5)
2010  format(10i8)

      end
