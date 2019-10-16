c$Id:$
      subroutine perspz(x,ix, ip, nen1,nen,ndm,numnp,numel)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Introduce 'ic' to do sort                        09/01/2012
c          Remove ze and zn from argument list
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Sort perspective projection of coordinates by z-values

c      Inputs:
c         x (ndm,*) - Nodal coordinates
c         ix(nen1,*)- Element connection list
c         ndm       - 1st dimension of x
c         nen1      - 1st dimension of ix
c         nen       - Number nodes connected to elements
c         numel     - Total number of elements
c         numnp     - Total number of nodal points

c      Outputs:
c         ze(numel,1) - Perspective z-coordinates
c         zn(numnp,1) - Perspective z-coordinates(destroyed)
c         ip(8,*)   - Element order for z-coordinates
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'
      include  'pbody.h'
      include  'ppers.h'
      include  'pdatay.h'
      include  'pointer.h'
      include  'comblk.h'

      include  'p_int.h'

      integer   i,ii, n,nen1,nen,ndm,numnp,numel,nsy,nne
      integer   ix(nen1,numel), ip(8,numel), ic(numel)
      real*8    t(3),x(ndm,numnp)
      real*4    zn(numnp),ze(numel)
      real*4    zmax,zpmin(8),zpmax(8)

      save

c     Loop over data points and find projection

      do nsy = 1,nsym

c       Reset block sorting array - isym

        isym(nsy) = nsy
        lsym = isym(nsy)

c       Reflect coordinates and compute depth coordinate - zn

        call pltsym(x,ndm,numnp,lsym)
        zpmax(lsym) = 0.0e0
        zpmin(lsym) = 0.0e0
        fp(1)       = npty - 1
        do n=1,numnp
          if(mr(fp(1)+n).ge.0) then
            do i=1,3
              t(i) = x(i,n) - e(i)
            end do ! i
            zn(n) = -real(xlbda(3,1)*t(1)
     &                  + xlbda(3,2)*t(2)
     &                  + xlbda(3,3)*t(3))
            zpmin(lsym) = max(zpmin(lsym),zn(n))
          else
            zn(n) = 0.0e0
          endif
        end do ! n

c       Search visible faces for depth sort

        do n = 1,nfac(lsym)
          nne = ip(lsym,n)
          if(nne.gt.0 .and. nne.le.numel) then
            if(ix(nen1-1,n).ge.nreg1 .and. ix(nen1-1,n).le.nreg2 .and.
     &                                     ix(nen1,nne).gt.0) then
              zmax  = 0.0e0
              do i = 1,nen
                ii = abs(ix(i,nne))
                if( ii.gt.0 .and. ii.le.numnp) then
                  zmax = max(zmax,zn(ii))
                endif
              end do ! i
              ze(nne)     = zmax
              zpmax(lsym) = max(zpmax(lsym),ze(nne))
              zpmin(lsym) = min(zpmin(lsym),ze(nne))
            endif
          endif
        end do ! n

c       Sort element plot order array 'ip' to produce hidden surface

        call merges ( -1, ze, 8, ip(lsym,1), nfac(lsym), ic )

c       Return coordinates to original values

        call pltsym(x,ndm,numnp,lsym)
      end do ! nsy

c     Sort blocks so panels on reflected mesh are hidden

      call merges ( -1, zpmin, 1, isym, nsym, ic )

      end
