c$Id:$
      subroutine umacr9(lct,ctl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c      1. Change call palloc to setvar = palloc             14/08/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Interface to PARAVIEW

c      Inputs:
c         lct       - Command character parameters
c         ctl(3)    - Command numerical parameters

c      Outputs:
c         N.B.  Users are responsible for command actions.  See
c               programmers manual for example.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'iofile.h'
      include  'pointer.h'
      include  'sdata.h'
      include  'strnum.h'
      include  'umac1.h'

      include  'comblk.h'

      logical   pcomp, setvar, palloc
      character lct*15, parafile*24,parext*4
      integer   i,ii,node, plu
      real*8    ctl(3)

      save

      data      plu / 99 /

c     Set command word

      if(pcomp(uct,'mac9',4)) then      ! Usual    form
        uct = 'pexp'                    ! Specify 'name'
      elseif(urest.eq.1) then           ! Read  restart data

      elseif(urest.eq.2) then           ! Write restart data

      else                              ! Perform user operation

        parext         = 'vtu '
        parafile       = '    '
        parafile(1:15) = lct(1:15)
        call addext(parafile,parext,24,4)

        write(*,*) 'Saving PARAVIEW data to ',parafile
        open(unit=plu,file=parafile,access='sequential')

c       Write out top header
        write(plu,1000)

        write(plu,1020) numnp,numel               ! Start Mesh/Piece Section
        write(plu,1010) '<Points>'                ! Start Point/Node data

        write(plu,1030) 'Float64','nodes',ndm
        do i = 1,numnp
          do ii = 1,ndm
            write(plu,2000) hr(np(43)+(i-1)*ndm+(ii-1))
          end do
        end do
        write(plu,1010) '</DataArray>'            ! Close Node data

        write(plu,1010) '</Points>'               ! Close Points section

        write(plu,1010) '<Cells>'                 ! Start Cell Section
        write(plu,1030) 'Int32','connectivity',1  ! Start Elements

c       Offsets memory allocation
        setvar = palloc(111,'TEMP1',numel+1,1)
        mr(np(111)) = 0;

        do i = 1,numel
          mr(np(111)+i) = mr(np(111)+i-1)
          do ii = 1,nen
            node = mr(np(33) + ii-1 + nen1*(i-1))
            if (node .ne. 0) then
              write(plu,2010) node-1
              mr(np(111)+i) = mr(np(111)+i) + 1
            endif
          end do
        end do

        write(plu,1010) '</DataArray>'           ! Close Elements

c       Output offsets for elements

        write(plu,1030) 'Int32','offsets',1      ! Start Offsets
        do i = 1,numel
          write(plu,2010) mr(np(111)+i)
        end do
        write(plu,1010) '</DataArray>'            ! Close Offsets

c       Output element connectivity type

        write(plu,1030) 'UInt8','types',1                  ! Start Element types
        do i = 1,numel
         if (mr(np(111)+i)-mr(np(111)+i-1) .eq. 2) then       ! 2 node line
           write(plu,2010) 3
         elseif (mr(np(111)+i)-mr(np(111)+i-1) .eq. 3) then   ! 3 node triangle
           write(plu,2010) 5
         elseif (mr(np(111)+i)-mr(np(111)+i-1) .eq. 4) then   ! 4 node quad
           write(plu,2010) 9
         elseif (mr(np(111)+i)-mr(np(111)+i-1) .eq. 8) then   ! 8 node brick
           write(plu,2010) 12
         endif
        end do

        write(plu,1010) '</DataArray>'              ! Close Element types
        write(plu,1010) '</Cells>'                  ! Close Cell Section

c       Delete memory used for offsets

        setvar = palloc(111,'TEMP1',0,1)

c       Output displacements

        write(plu,1010) '<PointData>'               ! Start Point Data

        write(plu,1030) 'Float64','Displ',ndf        ! Start Displacements
        do i = 1,numnp
          do ii =0,ndf-1
            write(plu,2000) hr(np(40)+(i-1)*ndf+ii)
          end do ! ii
        end do ! i
        write(plu,1010) '</DataArray>'              ! Close Displacements

c       Output stresses

        if(abs(istv).gt.0) then

          write(plu,1030) 'Float64','Stress',abs(istv)   ! Start Stresses
          do i = 1,numnp
            do ii =1,abs(istv)
              write(plu,2000) hr(np(58)+(i-1) + ii*numnp)
            end do
          end do
          write(plu,1010) '</DataArray>'              ! Close Stresses

          write(plu,1030) 'Float64','PStress',7       ! Start Principal Stresses
          do i = 1,numnp
            do ii =1,7
              write(plu,2000) hr(np(57)+(i-1) + ii*numnp)
            end do
          end do
          write(plu,1010) '</DataArray>'              ! Close Stresses
        else
          write(*,*) ' No stresses output to PARAVIEW file'
        endif

        write(plu,1010) '</PointData>'              ! Close Point Data Section

        write(plu,1010) '</Piece>'                  ! Close Mesh/Piece

c       Close the XML file

        write(plu,1010) '</UnstructuredGrid> </VTKFile>'
        close(plu, status = 'keep')

      endif

c     Formats

1000  format('<?xml version="1.0"?>',/
     &       '<VTKFile type="UnstructuredGrid" version="0.1">',/
     &       '<UnstructuredGrid>')

1010  format(a)

1020  format('<Piece NumberOfPoints="',i10,
     &       '" NumberOfCells="',i10,'">')

1030  format('<DataArray type="',a,'" Name="',a,
     &       '" NumberOfComponents="',i2,'" format="ascii">')

2000  format(1p,1e15.5,' ',$)
2010  format(i6,' ',$)

      end
