c$Id:$
      subroutine pcrdrd(x,ndtyp,prt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Input all coordinates

c      Inputs:
c         prt      - Print results if true

c      Outputs:
c         x(ndm,*) - Nodal coordinates
c         ndtyp(*) - Active node identifier
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'chdata.h'
      include   'comfil.h'
      include   'sdata.h'
      include   'iofile.h'
      include   'ioincl.h'
      include   'trdata.h'

      logical    prt, pcomp, incfl
      character  inam*4,fnam*15,gnam*4
      integer    i,j,n,noge, ibc, ndtyp(numnp)
      real*8     xx(3),x(ndm,numnp)

c     Check for an include or nogeneration record
c     [include filename <noge> - Input from filename
c     [nogeneration            - Omit generation field

      read(ior,1000) record
      call pstrip(xxx,record,1)
      call acheck(xxx,yyy,15,80,80)
      read(yyy,1002)  inam,fnam,gnam
      if(pcomp(inam,'incl',4)) then
        call pincld(fnam)
        noge  = 1
        incfl = .true.
        if(pcomp(gnam,'noge',4)) then
          noge = 0
        endif
      elseif(pcomp(inam,'noge',4)) then
        noge  = 0
        incfl = .false.
      else
        backspace(ior)
        noge  = 1
        incfl = .false.
      endif

c     List directed input of all the data

      do n = 1,numnp
        if(noge.eq.0) then
          read(ior,*) j,    (xx(i),i=1,ndm)
        else
          read(ior,*) j,ibc,(xx(i),i=1,ndm)
        endif
        do i = 1,ndm
          x(i,j) = xx(i)
        end do ! i
        ndtyp(j) = 0
      end do ! n

c     Transform if necessary

      if((tr(1,1)+tr(2,2)+tr(3,3) .ne. 3.d0) .or.
     &   (max(abs(xr(1)),abs(xr(2)),abs(xr(3))).ne.0.0d0) ) then
        do n = 1,numnp
          do i = 1,ndm
            xx(i) = x(i,n)
          end do ! i
          do i = 1,ndm
            x(i,n) = xr(i)
            do j = 1,ndm
              x(i,n) = x(i,n) + tr(i,j)*xx(j)
            end do ! j
          end do ! i
        end do ! n
      endif

c     Output computed data

      if(prt) then
        call prtitl(prt)
        write(iow,2000)  (i,i=1,ndm)
        do n = 1,numnp
          write(iow,2001) n,(x(i,n),i=1,ndm)
        end do ! n
      endif

      irecrd(isf) = irecrd(isf) + numnp

c     Clear the include file if necessary

      if(incfl) then
        call pincld('end')
      endif

c     Formats

1000  format(a)
1002  format(a4,11x,a15,a4)

2000  format(5x,'Nodal Coordinates'//6x,'Node',3(i6,' Coord':))
2001  format(i10,1p,3e12.4)

      end
