c$Id:$
      subroutine ptieix(id,ix,ntyp,f, isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Initilize 'ixfl' to blank                        02/03/2009
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Save element data in file 'IXxxxxx'

c      Inputs:
c         id(ndf,*)  - Boundary codes before ties
c         ix(nen1,*) - Element list to output (isw = 1)
c         ntyp(*)    - Node types
c         f(ndf,*,2) - Nodal force/displacements
c         nen1       - First  dimension
c         numel      - Number of elements
c         ndf        - Dof's/node
c         numnp      - Number of nodes

c      Outputs:
c         id(ndf,*)  - Boundary codes before ties
c         ix(nen1,*) - Element list to  input (isw = 2)
c         ntyp(*)    - Node types
c         f(ndf,*,2) - Nodal force/displacements
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'comfil.h'
      include   'cdata.h'
      include   'iodata.h'
      include   'sdata.h'
      include   'pointer.h'
      include   'comblk.h'

      logical    exst
      logical    setvar, palloc
      character  ixfil*129,fmac*6
      integer    isw, n,j
      integer    id(ndf,numnp),ix(nen1,numel),ntyp(numnp)
      real*8     f(ndf,numnp*2)

      fmac       = 'IXfile'
      ixfil      = ' '              ! Clear any erroneous characters
      ixfil(1:2) = 'IX'
      j          = index(finp,' ')
      ixfil(3:j) = finp(2:j-1)

      call opnfil('IXfile',ixfil,-1,ios,exst)

      if(isw.eq.1) then
        write(ios,'(10i8)') ((id(j,n),j=1,ndf),n=1,numnp)
        write(ios,'(10i8)') ((ix(j,n),j=1,nen1),n=1,numel)
        write(ios,'(10i8)') (ntyp(n),n=1,numnp)
        write(ios,'(1p,4e24.15)') ((f(j,n),j=1,ndf),n=1,numnp*2)
      else

        setvar = palloc(111,'TEMP1',numel*nen,1)
        call pcopyix(ix,mr(np(111)),nen1,nen,numel)
        read (ios,'(10i8)') ((id(j,n),j=1,ndf),n=1,numnp)
        read (ios,'(10i8)') ((ix(j,n),j=1,nen1),n=1,numel)
        read (ios,'(10i8)') (ntyp(n),n=1,numnp)
        read (ios,'(1p,4e24.15)') ((f(j,n),j=1,ndf),n=1,numnp*2)
        call psetu(ix,mr(np(111)),nen1,nen,numel,hr(np(40)),ndf)
        setvar = palloc(111,'TEMP1',0,1)
      endif

      close(ios,status = 'keep')

      end

      subroutine pcopyix(ix,ixc,nen1,nen,numel)
      implicit   none
      integer    nen1,nen,numel, n,j
      integer    ix(nen1,numel), ixc(nen,numel)

      do n = 1,numel
        do j = 1,nen
          ixc(j,n) = ix(j,n)
        end do ! j
      end do ! n

      end

      subroutine psetu(ix,ixc,nen1,nen,numel,u,ndf)

      implicit   none

      integer    nen1,nen,numel, ndf, n,i,j
      integer    ix(nen1,numel), ixc(nen,numel)
      real*8     u(ndf,*)
      do n = 1,numel
        do j = 1,nen
          if(ixc(j,n).ne.ix(j,n) .and.
     &       ixc(j,n).ne.0       .and.
     &       ix (j,n).ne.0       ) then  ! Node numbering is changed
            do i = 1,ndf
              u(i,ix(j,n)) = u(i,ixc(j,n))
            end do ! i
          endif
        end do ! j
      end do ! n

      end
