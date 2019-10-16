c$Id:$
      subroutine pcurve(x,id,xs,is,numsd,isd,ndm,ndf,numnp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Move tr,xr into and out of trb for call          11/06/2011
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Input boundary conditions, loads and displacements
c               along curve specified by a 'side'

c      Inputs:
c         x(ndm,*)   - Nodal coordinates
c         is(isd,*)  - Side descriptors
c         xs(3,*)    - Super node coordinates
c         numsd      - Number of sides to search
c         isd        - Dimension of is array.
c         ndm        - Number space dimensions
c         ndf        - Number dof/node
c         numnp      - Number nodes in mesh

c      Scratch:

c      Outputs:
c         id(ndf,*)  - Boundary condition restraint conditsions
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'bdata.h'
      include  'dstars.h'
      include  'iofile.h'
      include  'trdata.h'
      include  'p_point.h'
      include  'pointer.h'
      include  'comblk.h'

      logical   errck,palloc, pcomp,tinput,sidfl,flag
      character text*15
      integer   numsd,isd,ndm,ndf,numnp
      integer   n, i, ii,jj,is(isd,*), id(ndf,*)
      integer   s1,s2, ni,ns, node
      real*8    dd,ddmin
      real*8    x(ndm,*), xs(3,*)
      real*8    td(16), trb(3,4)

      save

c     Loop through records

      write(iow,2000) head
      text = 'start'
      do while (.not.pcomp(text,'    ',4))
        errck = tinput(text,1,td,15)
        if(.not.pcomp(text,'    ',4)) then
          s1 = nint(td(1))
          s2 = nint(td(2))
          ni = nint(td(3))
          write(iow,2001) text(1:10),s1,s2,ni
          sidfl = .false.
          do i = 1,numsd
            if    (is(2,i).eq.s1 .and. is(3,i).eq.s2) then
              ii    = i
              sidfl = .true.
              exit
            elseif(is(3,i).eq.s1 .and. is(2,i).eq.s2) then
              ii    = -i
              sidfl = .true.
              exit
            endif
          end do ! i

          if(sidfl) then
            ns = abs(ii)
            do jj = 13,4,-1
              if(is(jj,ns).ne.0) then
                exit
              endif
            end do ! jj
            jj = jj - 1
            errck = palloc ( 111, 'TEMP1', (ni+1)     , 2)
            errck = palloc ( 112, 'TEMP2', jj*ndm     , 2)
            errck = palloc ( 113, 'TEMP3', (ni+1)*ndm , 2)
            do i = 1,3
              do ii = 1,3
                trb(i,ii) = tr(i,ii)
              end do ! ii
              trb(i,4) = xr(i)
            end do ! i
            call pside1(ni,xs,trb,ii,is(2,ns),jj,ndm,hr(np(111)),
     &                  hr(np(112)), hr(np(113)), is(1,ns))
            do i = 1,3
              do ii = 1,3
                tr(i,ii) = trb(i,ii)
              end do ! ii
              xr(i) = trb(i,4)
            end do ! i

            do i = 1,ni+1
              point = np(113) + ndm*(i - 1)  - 1
              flag = .true.
              do n = 1,numnp
                dd = 0.0d0
                do ii = 1,ndm
                  dd = dd + (x(ii,n) - hr(point+ii))**2
                end do ! ii
                if(flag) then
                  ddmin = dd
                  node  = n
                  flag  = .false.
                elseif(dd.lt.ddmin) then
                  ddmin = dd
                  node  = n
                endif
              end do ! n
              if(pcomp(text,'boun',4)) then
                do ii = 1,ndf
                  id(ii,node) = abs(id(ii,node)) + nint(td(ii+3))
                end do ! ii
              endif
            end do ! i

c           Destroy temporary arrays

            errck = palloc ( 113, 'TEMP3',     0      , 2)
            errck = palloc ( 112, 'TEMP2',     0      , 2)
            errck = palloc ( 111, 'TEMP1',     0      , 2)
          endif
        endif
      end do ! while

c     Formats

2000  format(1x,19a4,a3//'   C U R V E    C o n d i t i o n s'//
     &       12x,'Type',4x,'S_node 1',4x,'S_node 2',2x,'Increments')
2001  format(5x,a10,3i12)

      end
