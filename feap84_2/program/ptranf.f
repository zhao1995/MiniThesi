c$Id:$
      subroutine ptranf(tx,prt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add check on zero input of transformation        21/08/2007
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set coordinate transformation array and offsets
c               x_new = tr * x_old + xr

c      Inputs:
c         tx*(*)   - Character identification for transforming
c         prt      - Logical print control flag

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'iofile.h'
      include   'trdata.h'

      logical    prt, pcomp, errck, pinput
      character  tx*(*)
      integer    j
      real*8     td(3),tro(3,3),xro(3)

      save

c     [tran],<inc> - Specify coordinate transformation array

c     Incremental rotation update: tr_new = tinc * tr_old
c                                  xr_new = tinc * xr_old + xr_inc

      if(pcomp(tx,'inc',3)) then
        do j = 1,3
          tro(j,1) = tr(j,1)
          tro(j,2) = tr(j,2)
          tro(j,3) = tr(j,3)
          xro(j)   = xr(j)
        end do ! j
  1     errck = pinput(td,3)
        if(errck) go to 1
        if(abs(td(1))+abs(td(2))+abs(td(3)).eq.0.0d0) then
          write(ilg,3000)
          write(iow,3000)
          call plstop()
        endif
        do j = 1,3
          tr(1,j) = td(1)*tro(1,j) + td(2)*tro(2,j) + td(3)*tro(3,j)
        end do ! j
        xr(1) = td(1)*xro(1) +  td(2)*xro(2) +  td(3)*xro(3)
        errck = pinput(td,3)
        if(abs(td(1))+abs(td(2))+abs(td(3)).eq.0.0d0) then
          write(ilg,3000)
          write(iow,3000)
          call plstop()
        endif
        do j = 1,3
          tr(2,j) = td(1)*tro(1,j) + td(2)*tro(2,j) + td(3)*tro(3,j)
        end do ! j
        xr(2) = td(1)*xro(1) +  td(2)*xro(2) +  td(3)*xro(3)
        errck = pinput(td,3)
        if(abs(td(1))+abs(td(2))+abs(td(3)).eq.0.0d0) then
          write(ilg,3000)
          write(iow,3000)
          call plstop()
        endif
        do j = 1,3
          tr(3,j) = td(1)*tro(1,j) + td(2)*tro(2,j) + td(3)*tro(3,j)
        end do ! j
        xr(3) = td(1)*xro(1) +  td(2)*xro(2) +  td(3)*xro(3)
        errck = pinput(td,3)
        do j = 1,3
          xr(j) = xr(j) + td(j)
        end do ! j

c     Total rotation set: tr_new & xr_new set from inputs

      else
  2     errck = pinput(xr,3)
        if(errck) go to 2
        if(abs(xr(1))+abs(xr(2))+abs(xr(3)).eq.0.0d0) then
          write(ilg,3000)
          write(iow,3000)
          call plstop()
        endif
        do j = 1,3
          tr(1,j) = xr(j)
        end do ! j
        errck = pinput(xr,3)
        if(abs(xr(1))+abs(xr(2))+abs(xr(3)).eq.0.0d0) then
          write(ilg,3000)
          write(iow,3000)
          call plstop()
        endif
        do j = 1,3
          tr(2,j) = xr(j)
        end do ! j
        errck = pinput(xr,3)
        if(abs(xr(1))+abs(xr(2))+abs(xr(3)).eq.0.0d0) then
          write(ilg,3000)
          write(iow,3000)
          call plstop()
        endif
        do j = 1,3
          tr(3,j) = xr(j)
        end do ! j
        errck = pinput(xr,3)
      endif
      trdet = tr(1,1)*(tr(2,2)*tr(3,3) - tr(2,3)*tr(3,2))
     &      + tr(1,2)*(tr(2,3)*tr(3,1) - tr(2,1)*tr(3,3))
     &      + tr(1,3)*(tr(2,1)*tr(3,2) - tr(2,2)*tr(3,1))
      if(prt) then
        call mprint(tr,3,3,3,'Coord. T_matrix')
        call mprint(xr,1,3,1,'Coord. X_vector')
      endif

c     Format

3000  format(/' *ERROR in TRANsform* Entries are all zero')

      end
