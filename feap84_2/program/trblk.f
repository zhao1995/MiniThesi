c$Id:$
      subroutine trblk(nr,xl,ixl,x,ix,ndm,nod1,nuel1,nen1,ma,ntyp,
     &                 ctype,prt,prth)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Generate triangular block of 3-node triangular elements

c      Inputs:
c         nr        - Number elements in 1-local coordinate dir.
c         xl(ndm,*) - Block nodal coordinate array
c         ixl(*)    - Block nodal connection list
c         ndm       - Spatial dimension of mesh
c         nod1      - Initial node number for block
c         nuel1     - Initial element number for block
c         nen1      - Dimension of ix array
c         ma        - Material set number for block
c         ntyp      - Element type: 1 = 3-node; 2 = 6-node
c         ctype     - Input coordinate types
c         prt       - Output generated data if true
c         prth      - Output title/header data if true

c      Outputs:
c         x(ndm,*)  - Nodal coordinates for block
c         ix(*)     - Element nodal connection list for block
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'cdat2.h'
      include  'iofile.h'
      include  'trdata.h'

      logical   prt,prth, pcomp
      character xh*6, ctype*15
      integer   ni,nn,mct,k,i,j,i1,j1,nei,n1,n2,n3
      integer   nr,ndm,nod1,nuel1,nen1,ma,ntyp, ixl(1),ix(nen1,1)
      real*8    dl, rr,sn,cn, xl(3,1),x(ndm,1),xx(3),tshp(6),el(3)

      save

      data      xh/' coord'/

      ni  = nr
      mct = 0

      do i=1,3
        j = mod(i,3) + 1
        if((ixl(i+3)).ne.0) then
          do i1=1,3
            xl(i1,i+3) = xl(i1,i+3) - 0.5d0*(xl(i1,i) + xl(i1,j))
          end do ! i1
        endif
        xx(i) = 0.0d0
      end do ! i

c     Generate nodes

      nn = nod1
      dl  = 1.0d0/ni
      do i=0,ni
        el(3) = i*dl

        do j=0,ni-i
          el(2) = j*dl
          el(1) = 1.0d0 - el(3) - el(2)

c         Form shape functions

          do i1=1,3
            j1 = mod(i1,3) + 1
            tshp(i1)   = el(i1)
            tshp(i1+3) = 4.0d0*el(i1)*el(j1)
          end do ! i1

          if(nn.gt.numnp) then
            write(*,*) ' Trying to generate node',nn
            call plstop()
          endif
          do i1 = 1,ndm
            xx(i1) = 0.0d0
            do k=1,6
              xx(i1) = xx(i1) + tshp(k)*xl(i1,k)
            end do ! k
          end do ! i1
          if(pcomp(ctype,'pola',4) .or. pcomp(ctype,'cyli',4)) then
            call pdegree(xx(2), sn,cn)
            rr    = xx(1)
            xx(1) = x0(1) + rr*cn
            xx(2) = x0(2) + rr*sn
          endif
          do k = 1,ndm
            x(k,nn) = xr(k)+tr(k,1)*xx(1)+tr(k,2)*xx(2)+tr(k,3)*xx(3)
          end do ! k
          if(prt) then
            mct = mct + 1
            if(mod(mct,50).eq.1) then
              call prtitl(prth)
              write(iow,2003) (k,xh,k=1,ndm)
              if(ior.lt.0) write(*,2003) (k,xh,k=1,ndm)
            endif
            write(iow,2004) nn,(x(k,nn),k=1,ndm)
            if(ior.lt.0) write(*,2004) nn,(x(k,nn),k=1,ndm)
          endif
          nn = nn + 1
        end do ! j
      end do ! i

c     Generate elements

      nei = nuel1
      n1 = nod1
      do i=1,ni,ntyp
        do j=1,ni-i+1,ntyp
          if(ntyp.eq.1) then
            n2           = n1 + ni - i + 2
            ix(1,nei)    = n1 + j - 1
            ix(2,nei)    = n1 + j
            ix(3,nei)    = n2 + j - 1
            ix(nen1,nei) = ma
            if(j.lt.(ni-i+1)) then
              nei = nei + 1
              ix(1,nei)    = n2 + j - 1
              ix(2,nei)    = n1 + j
              ix(3,nei)    = n2 + j
              ix(nen1,nei) = ma
            endif
          elseif(ntyp.eq.2) then
            n3           = n1 + ni - i + 2
            n2           = n3 + ni - i + 1
            ix(1,nei)    = n1 + j - 1
            ix(4,nei)    = n1 + j
            ix(2,nei)    = n1 + j + 1
            ix(6,nei)    = n3 + j - 1
            ix(5,nei)    = n3 + j
            ix(3,nei)    = n2 + j - 1
            ix(nen1,nei) = ma
            if(j.lt.(ni-i)) then
              nei          = nei + 1
              ix(1,nei)    = n2 + j - 1
              ix(4,nei)    = n3 + j
              ix(2,nei)    = n1 + j + 1
              ix(5,nei)    = n3 + j + 1
              ix(6,nei)    = n2 + j
              ix(3,nei)    = n2 + j + 1
              ix(nen1,nei) = ma
            endif
          endif
          nei = nei + 1
        end do ! j
        n1 = n2
      end do ! i

c     Formats

2003  format(/'  N o d a l   C o o r d i n a t e s'//6x,'Node',5(i7,a6))

2004  format(i10,5f13.4)

      end
