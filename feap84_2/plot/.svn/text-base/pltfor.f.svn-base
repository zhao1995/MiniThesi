c$Id:$
      subroutine pltfor(x,f,id,ip,ndm,ndf,numnp,n1,ct, isgn)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c      1. Change fms to fms(1)                              15/08/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Draw vectors for forces on mesh

c      Inputs:
c         x(ndm,*)  - Nodal coordinates for mesh
c         f(ndf,*)  - Nodal forces
c         id(ndf,*) - Boundary condition indicator array
c         ip(*)     - Active node indicator
c         ndm       - Dimension of x array
c         ndf       - Dimension of f and id arrays
c         numnp     - Number of nodes in mesh
c         n1        - Flag, place tip at node if > 0
c         isgn      - Flag, plot displacements if > 1

c      Outputs:
c         none      - Plot outputs to screen/file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'
      include  'pdata1.h'
      include  'pdata4.h'
      include  'pdatay.h'
      include  'pdatxt.h'
      include  'pointer.h'
      include  'rigid1.h'
      include  'comblk.h'

      logical   vfl,zoom, fdis, frig
      integer   ndm,ndf,numnp,n1,isgn,nsy, i,j,n
      real*8    fm,dx1,dx2,dx3,d, ct

      integer   id(ndf,*),ip(*)
      real*8    dd(3),xx(3,4),x(ndm,*),f(ndf,*), tbuf(2), fms(1)

      save

c     Compute longest vector

      fdis   = isgn.gt.1
      if(ct.le.0.0d0) then
        frig   = fdis
        fm   = 0.d0
        do n = 1,numnp
          if(ip(n).gt.0 .and. mr(npty+n-1).ge.0
     &                  .and. zoom(x(1,n),ndm)) then
            if(rbody) then
              frig = fdis .or. mr(np(100)+n-1).ne.0
            endif
            d = 0.d0
            do i = 1,min(ndm,ndf)
              j = pdf(i)
              if(j.gt.0 .and. j.le.ndf) then
                if(id(j,n).gt.0 .or. frig .or. isgn.lt.0) then
                  d = d + f(j,n)**2
                end if
              end if
            end do ! i
            fm = max(fm,d)
          endif
        end do ! n

c       Parallel exchange to set common scale

        fms(1) = fm
        call pfeapsr(fms,tbuf, 1, .false.)
        fm = max(fm,fms(1))


c     Set to specified length

      else
        fm = ct*ct
      endif

c     Zero length vectors

      if(fm.le.0.0d0) then
        if(iow.lt.0) write(*,2000)

c     Compute vector at each node

      else

        do i = 1,3
          xx(i,1) = 0.0d0
          xx(i,2) = 0.0d0
          xx(i,3) = 0.0d0
          xx(i,4) = 0.0d0
        end do ! i

        frig = fdis
        fm   = isgn*sqrt(fm)*scale*40.d0
        do nsy = 1,nsym

          do n = 1,numnp
            if(ip(n).gt.0 .and. zoom(x(1,n),ndm)) then
              if(rbody) then
                frig = mr(np(100)+n-1).ne.0
              endif
              vfl = .false.
              do i = 1,3
                dd(i) = 0.0d0
              end do ! i
              do i = 1,min(ndm,ndf)
                j = pdf(i)
                if(j.gt.0 .and. j.le.ndf) then
                  if( ((id(j,n).gt.0) .or. frig .or. isgn.lt.0)
     &                .and. (f(j,n).ne.0.0d0) .or. fdis ) then
                    dd(i) = f(j,n)
                    vfl = .true.
                  endif
                endif
              end do ! i
              if(vfl) then
                dd(1) = dd(1)/fm
                dd(2) = dd(2)/fm
                if(ndm.ge.3) then
                  dd(3)   = dd(3)/fm
                  xx(3,1) = x(3,n)
                  xx(3,2) = xx(3,1) + dd(3)
                  xx(3,3) = xx(3,2) -.6d0*dd(3) + .2d0*(dd(1)+dd(2))
                  xx(3,4) = xx(3,2) -.6d0*dd(3) - .2d0*(dd(1)+dd(2))
                endif
                xx(1,1) = x(1,n)
                xx(2,1) = x(2,n)
                xx(1,2) = xx(1,1) + dd(1)
                xx(2,2) = xx(2,1) + dd(2)
                xx(1,3) = xx(1,2) -.6d0*dd(1) - .2d0*(dd(2)+dd(3))
                xx(2,3) = xx(2,2) -.6d0*dd(2) + .2d0*(dd(1)+dd(3))
                xx(1,4) = xx(1,2) -.6d0*dd(1) + .2d0*(dd(2)+dd(3))
                xx(2,4) = xx(2,2) -.6d0*dd(2) - .2d0*(dd(1)+dd(3))

c               Plot vector (noting correct symmetry)

                do i = 1,4
                  do j = 1,3
                    xx(j,i) = (xx(j,i) - xsyc(j))*xsym(j,nsy) + xsyc(j)
                  end do ! j
                end do ! i
                if(n1.eq.0) then
                  dx1 = 0.d0
                  dx2 = 0.d0
                  dx3 = 0.d0
                else
                  dx1 = xx(1,1) - xx(1,2)
                  dx2 = xx(2,1) - xx(2,2)
                  dx3 = xx(3,1) - xx(3,2)
                endif
                call plotl(xx(1,1)+dx1,xx(2,1)+dx2,xx(3,1)+dx3,3)
                call plotl(xx(1,2)+dx1,xx(2,2)+dx2,xx(3,2)+dx3,2)
                call plotl(xx(1,3)+dx1,xx(2,3)+dx2,xx(3,3)+dx3,2)
                call plotl(xx(1,4)+dx1,xx(2,4)+dx2,xx(3,4)+dx3,2)
                call plotl(xx(1,2)+dx1,xx(2,2)+dx2,xx(3,2)+dx3,2)
              endif
            endif
          end do ! n
        end do ! nsy
      endif

2000  format('  Zero values acting on mesh ')

      end
