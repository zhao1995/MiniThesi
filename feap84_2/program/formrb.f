c$Id:$
      subroutine formrb(ul,xl,tl,p,s,ie,d,x,ix,rben,rlam,eqrb,
     &                  nrbody,neqrb)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute rigid body properties in initial configuration

c      Inputs:
c         ul(ndf,*)  - Element solution array
c         xl(ndm,*)  - Element nodal coordinate array
c         tl(*)      - Element nodal temperature array
c         p(*)       - Element lumped mass array
c         s(nst,nst) - Element consistent mass array
c         ie(nie,*)  - Assembly information for material sets
c         d(*)       - Material set parameters
c         x(ndm,*)   - Nodal coordinates
c         ix(nen1,*) - Element nodal connection array
c         eqrb(*)    - Type of dynamic update for each rigid body
c         nrbody     - Number of rigid bodies
c         neqrb      - Type of update for current rigid body

c      Outputs:
c         rlam(*)    - Initial rotational quantities (Lambda)
c                      Also returns mass and inertia properties through
c                      pointers
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'cdat1.h'
      include  'eldata.h'
      include  'evdata.h'
      include  'iofile.h'
      include  'pointer.h'
      include  'sdata.h'
      include  'comblk.h'

      include  'p_int.h'

      integer   i, j, ii, irb, nrbody, neqrb
      integer   ie(nie,*),ix(nen1,numel),rben(*),eqrb(nrbody)
      real*8    ul(ndf,nen),xl(ndm,nen),tl(nen),p(nst),s(nst,nst)
      real*8    d(ndd,*),x( ndm,numnp), xmascc(3), jmascc(3,3), xcg(3)
      real*8    rlam(9,6,nrbody), tmasc

      save

c     Zero rigid body arrays

      call pzero(hr(np(107)),nrbody)
      call pzero(hr(np( 98)),nrbody*9)
      call pzero(hr(np( 95)),nrbody*33)

      call pzero(xcg , 3 )
      call pzero(rlam,54*nrbody)

c     Loop over rigid bodies

      do irb = 1, nrbody

c       Set lambda to unit matrix

        if(neqrb.eq.-6 .or. neqrb.eq.-7) then

          do i = 1,3
            rlam(1,i,irb) = 1.d0
            rlam(5,i,irb) = 1.d0
            rlam(9,i,irb) = 1.d0
          end do ! i

c       Set lambda and relative rotation to unit quaternion

        else

          do i = 1,3
            rlam(4,i,irb) = 1.0d0
            rlam(8,i,irb) = 1.0d0
          end do ! i

        endif

c       Set integration type: quaternion update

        eqrb(irb) = neqrb

        do n = 1,numel

c         Check if element is a part of rigid body,irb

          if(rben(n).eq.irb) then

            call pzero(s,nst*nst)
            call pzero(p,nst)
            call pzero(xl,ndm*nen)
            call pzero(ul,ndf*nen)
            call pzero(tl,nen)

            do i = 1,nen
              ii = ix(i,n)
              if(ii.gt.0) then
                nel = i
                do j = 1,ndm
                  xl(j,i) = x(j,ii)
                end do ! j
              end if
            end do ! i

c           Determine material number for element
            ma = ix(nen1,n)

c           Determine material type for element
            iel = ie(nie-1,ma)

c           Compute element mass matrix

            imtyp = 1
            iel   = ie(nie-1,ma)
            call elmlib(d(1,ma),ul,xl,ix(1,n),tl,s,p,
     &                  ndf,ndm,nst,iel,5)

c           Compute element mass, x-bar*mass, inertia tensor

            call rmas3d(xl,s,xmascc,jmascc,xcg,tmasc,ndf,ndm,nel,nst)

c           Assemble terms

            call rassbl(hr(np(107)),hr(np(95)),hr(np(98)),
     &                  tmasc,xmascc,jmascc,irb,ndm)

          end if
        end do ! n

      end do ! irb

c     Compute center of mass and inertia tensor at center of mass

      do i = 1,nrbody
        call rmascg(hr(np(107)),hr(np(95)),hr(np(98)),xcg,i,ndm)
      end do ! i

c     Output rigid body mass and center of mass

      fp(1) = np(95) - 1
      write(iow,2000) (i,i=1,ndm)
      do i = 0,nrbody-1
        write(iow,2001) i+1,hr(np(107)+i),(hr(fp(1)+j),j=1,ndm)
        fp(1) = fp(1) + 33
      end do ! i

c     Output rigid body inertia tensor

      write(iow,2002)
      do i = 0,nrbody-1
        fp(1) = np(98) + 9*i - 1
        if(ndm.eq.2) then
          write(iow,2001) i+1,hr(fp(1)+9)
        elseif(ndm.eq.3) then
          fp(1) = np(98) + 9*i - 1
          write(iow,2003) i+1,(hr(fp(1)+j),j=1,9)
        endif
      end do ! i

c     Formats

2000  format(/'     Body Total Mass',3('  x-',i1,' center':))

2001  format(/i8,1p,4e12.4)

2002  format(/'     Body       Inertia Tensor')

2003  format(/i8,1p,3e12.4/(8x,1p,3e12.4))

      end
