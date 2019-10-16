c$Id:$
      subroutine pnforc(id,fpro,cforc,jd,a,al,dr,x,u,
     &                  alfl,aufl,bfl,dfl,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Assemble nodal radial follower forces

c      Inputs:
c         id(*)     - Equation numbers of active dof
c         cforc(*)  - Force values
c         fpro(*)   - DOF proportional load values
c         jd(*)     - Column pointers
c         x(*)      - Nodal coordinate values
c         u(*)      - Solution values
c         alfl      - Assemble lower part if true
c         aufl      - Assemble upper part if true
c         bfl       - Assemble residual if true
c         dfl       - Flag, assemble uncompressed residual if true
c         isw       - Switch on action to perform

c      Outputs:
c         a(*)      - Diagonal/upper part of tangent array
c         al(*)     - Lower part of tangent array
c         dr(*)     - Residual array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'eldata.h'
      include  'sdata.h'
      include  'prld1.h'
      include  'prlod.h'


      logical   alfl,aufl,bfl,dfl,gfl,rel
      integer   nov,isw
      integer   id(ndf,numnp),fpro(ndf,numnp),jd(*)
      real*8    cforc(ndf,numnp)
      real*8    a(*),al(*),dr(*),x(ndm,*),u(ndf,*)

      integer   i,nn, nss
      integer   ld(3), ix(1)
      real*8    rad,cs,sn
      real*8    s(3,3),p(3), xx(3),ff(3)

      save

c     Zero unused part of s-array

      do i = 1,3
        s(i,3) = 0.0d0
        s(3,i) = 0.0d0
      end do ! i

c     Modify tangent: diagonal only

      if((isw.eq.3) .or. (isw.eq.6)) then

        nss = min(ndm,ndf,3)   ! Allow for only 3-force components

        do nn = 1,numnp
          nov = 0                ! Counts elements when isw.eq.8
          rel   = .false.
          gfl   = .true.
          nel   = 1
          ix(1) = nn
          do i = 1,nss
            ld(i) = id(i,nn)          ! Set assembly locations
            xx(i) = x(i,nn) + u(i,nn) ! Set deformed coordinates
            ff(i) = cforc(i,nn)       ! Force comp (need to multiply)
          end do ! i
          do i = nss+1,3
            ld(i) = 0                 ! Mask assembly for 2-d problems
            xx(i) = 0.0d0
            ff(i) = 0.0d0
          end do ! i

c         Multiply by correct proportional load

          do i = 1,nss
            if(fpro(i,nn).gt.0) then
              ff(i) = ff(i)*prldv(i)
            else
              ff(i) = ff(i)*prop
            endif
          end do ! i

          if(max(abs(ff(1)),abs(ff(2)),abs(ff(3))).gt.0.0d0) then

c           Geometric values

            rad = sqrt(xx(1)**2 + xx(2)**2)
            if(rad.gt.0.0d0) then
              cs  = xx(1)/rad
              sn  = xx(2)/rad
              rad = 1.d0/rad
            else
              cs  = 1.0d0
              sn  = 0.0d0
              rad = 1.0d0
            endif

c           Form residual

            p(1)   = cs*ff(1) - sn*ff(2)
            p(2)   = sn*ff(1) + cs*ff(2)
            p(3)   = ff(3)

c           Form tangent

            s(1,1) = -sn*p(2)*rad
            s(1,2) =  cs*p(2)*rad
            s(2,1) =  sn*p(1)*rad
            s(2,2) = -cs*p(1)*rad

c           Assemble

            call passble(s,p,ld,ix,jd,a,al,dr,
     &                   alfl,aufl,bfl,dfl,gfl,rel, 3,nov, isw)
          endif
        end do ! nn

c     Compute momentum and energy

      elseif(isw.eq.13) then


      endif

      end
