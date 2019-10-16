c$Id:$
      subroutine plinlod(ntyp,vv,x,f,prt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add 'ntyp' test to ensure node is active         21/04/2007
c       2. Increase 'ndm + ndf' to 31                       05/07/2007
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Compute values of nodal loads for line load

c     Inputs:
c        ntyp(*)  - Type of loading
c        vv(3)    - Orientation vector for line loading
c        x(ndm,*) - Nodal coordinate array
c        prt      - Print flag: Output values if 'true'

c     Outputs:
c        f(ndf,*) - Nodal loading for a line load
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'iofile.h'
      include   'sdata.h'

      logical    prt, errck, tinput, flag
      character  num*15
      integer    i,j,n
      integer    ntyp(*)
      real*8     x(ndm,*), f(ndf,*)
      real*8     len,lenx,lenn,gap, xi
      real*8     td(31),xx(3,2),forc(12,2)
      real*8     dx(12),ds(3),n1(3),n2(3),vv(3), g(3)

      save

      data      gap  / 1.0d-6 /

c     Line loading

      len   = 0.0d0
      do i = 1,ndm
        len = len + vv(i)*vv(i)
      end do ! i
      if(len.gt.0.0d0) then
        len = 1.d0/sqrt(len)
        do i = 1,ndm
          vv(i) = vv(i)*len
        end do ! i
      else
        write(iow,*) ' ERROR: Reference vector zero'
        call plstop()
      endif

c     Input end of line segment and load values
      do j= 1,2
        errck = tinput(num,1,td,min(ndm+ndf,15))
        if(ndm+ndf.gt.15) then
          errck = tinput(num,1,td(16),ndm+ndf-15)
        endif
        do i = 1,ndm
          xx(i,j) = td(i)
        end do ! i
        do i = 1,ndf
          forc(i,j) = td(i+ndm)
        end do ! i
      end do ! j

      do i = 1,ndf
        forc(i,2) = forc(i,2) - forc(i,1)
      end do ! i

c     Compute unit vector along input segment

      lenx = 0.0d0
      do i = 1,ndm
        ds(i) = xx(i,2) - xx(i,1)
        lenx   = lenx + ds(i)*ds(i)
      end do
      if(lenx.gt.0.0d0) then
        lenx = sqrt(lenx)
        len  = 1.d0/lenx
        do i = 1,ndm
          n1(i) = ds(i)*len
        end do ! i
      else
        write(iow,*) ' ERROR: Line has zero length'
        call plstop()
      endif

c     Compute cross product

      if(ndm.eq.3) then
        g(1) = n1(2)*vv(3) - n1(3)*vv(2)
        g(2) = n1(3)*vv(1) - n1(1)*vv(3)
        g(3) = n1(1)*vv(2) - n1(2)*vv(1)
        len  = g(1)*g(1) + g(2)*g(2) + g(3)*g(3)
        if(len.gt.0.0d0) then
          len = 1.d0/sqrt(len)
          do i = 1,3
            g(i) = g(i)*len
          end do ! i
        else
          write(iow,*) ' ERROR: Normal parallel to line'
          call plstop()
        endif
      elseif(ndm.eq.2) then
        g(1) =  n1(2)
        g(2) = -n1(1)
        g(3) =  0.0d0
      endif

c     Loop over nodes to find ones on segment

      if(prt) then
        call prtitl(.true.)
        write(iow,2000) (i,i=1,ndf)
        if(ior.lt.0) then
          write(*,2000) (i,i=1,ndf)
        endif
      endif
      do n = 1,numnp

        if(ntyp(n).ge.0) then
c         Compute unit vector from line segment to node

          flag = .false.
          lenn = 0.0d0
          do i = 1,ndm
            dx(i) = x(i,n) - xx(i,1)
            lenn  = lenn + dx(i)*dx(i)
          end do ! i
          if(lenn.gt.0.0d0) then
            lenn = sqrt(lenn)
            len  = 1.d0/lenn
            do i = 1,ndm
              n2(i) = dx(i)*len
            end do ! i

c           Compute distance from node to line segment

            len = 0.0d0
            do i = 1,ndm
              len = len + g(i)*n2(i)
            end do ! i
            if(abs(len).le.gap) then
              flag = .true.
            endif

c         Node at end point

          else
            flag = .true.
          endif

c         Interpolate force vector

          if(flag) then
            xi = lenn/lenx
            do i = 1,ndf
              dx(i)  = forc(i,1) + xi*forc(i,2)
              f(i,n) = f(i,n) + dx(i)
            end do ! i
            if(prt) then
              write(iow,2001) n,(dx(i),i=1,ndf)
              if(ior.lt.0) then
                write(*,2001) n,(dx(i),i=1,ndf)
              endif
            endif
          endif ! flag
        endif ! ntyp
      end do ! n

c     Formats

2000  format(5x,'Nodal Forces From Line Inputs'//
     &       4x,'Node',6(i6,'-Force'):/(8x,6(i6,'-Force')))
2001  format(i8,1p,6e12.4:/(8x,1p,6e12.4))

      end
