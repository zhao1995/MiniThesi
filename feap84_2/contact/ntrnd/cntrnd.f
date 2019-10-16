c$Id:$
      subroutine cntrnd(ndm,ndf,x,u,csw,npair,cs02,cp0,
     &                  ix1,ix2,ch1,ch2,ch3,ww1,ww3)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Robert L. Taylor        November 19, 2001            1.0

c      Acronym: Contact DRIVER for 2D slave node to rigid surface

c      Purpose: Management of specific contact formulation

c      Inputs :
c         ndm     - Space dimension of mesh
c         ndf     - Number dof/node
c         x(*)    - Nodal coordinates
c         u(*)    - Current nodal solution vectors
c         csw     - Contact switch
c         npair   - # of current pair
c         cs02(*) - Contact surface 2 control data
c         cp0(*)  - Contactpair control data
c         ix1(*)  - Element nodal connection list for surface 1
c         ix2(*)  - Element nodal connection list for surface 2

c      Outputs:
c         ch1(*)  - Contact history variables (old)
c         ch2(*)  - Contact history variables (current)
c         ch3(*)  - Contact history variables (static)
c         ww1(*)  - Dictionary of variables for CH1 & CH2
c         ww3(*)  - Dictionary of variables for CH3
c                 - Data exchange with main program subroutine calls
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_comnd.h'
      include  'c_contac.h'
      include  'c_geom.h'
      include  'c_keyh.h'
      include  'c_pair.h'
      include  'augdat.h'
      include  'iofile.h'
      include  'ptdat4.h'
      include  'pointer.h'
      include  'comblk.h'

      logical   ifprt,once, compg
      character ww1(*)*(*),ww3(*)*(*)
      integer   ndm,ndf, csw,npair
      integer   ix1(dnope1,*),ix2(dnope2,*), ixl(1),ida(3),ilm(3)
      integer   ke,kn,nel1,nod1,ns,kset,istgn,fel,lel,n1,i,j
      integer   dopt, lcnt,liter
      real*8    cs02(nr0,n0c1:*),cp0(nr0,n0c3:*), cn(4)
      real*8    ch1(lh1,*),ch2(lh1,*),ch3(lh3,*), x(ndm,*),u(ndf,*)
      real*8    cxs(3),tanm(4,4),resv(4)

      save

c     Set active dof and dof order (idl length = ndf)

      data      ida     /1,2,3/
      data      ilm     /0,0,0/

      call cdebug0 ('    cntrnd',csw)

c-----[--.----+----.----+----.-----------------------------------------]
      if (csw.eq.0) then
        once = .true.

c     Called from PCONTR for activation of requested history variables

      elseif (csw.eq.1) then

c       Once actions

        if (once) then
          once = .false.

c         Load dictionary of history variables

          call defhvptr (ww1,ww3)

c         Print a warning if unsymmetric solver is needed

          if (iffric.eq.1) then
            write (*,3000)
          endif
        endif

c-----[--.----+----.----+----.-----------------------------------------]
c       N.B. Initial check automatically carried out with csw = 14

      elseif (csw.eq.2) then
        continue

c-----[--.----+----.----+----.-----------------------------------------]
c       Called from FORMFE to compute stiffness and residual

      elseif (csw.eq.3) then

c       Loop over all surface 1 elements

        if(max(cp0(4,0),abs(cp0(5,0))).eq.0.0d0) then
          cn(1) = cp0(3,2)
          cn(2) = 0.0d0
        else
          cn(1) = cp0(4,0)
          cn(2) = cp0(5,0)
        endif
        do ke = 1,neps1
          nel1 = ke
          do kn = 1,nope1
            nod1 = kn

c           Skip if node just checked

            if (dnope1.eq.1 .or. nope1.eq.4 .or.
     &          (nod1.eq.1) .or. (ix1(dnope1,nel1).eq.0)) then

              kset = nel1
              if (nope1.eq.4) then
                kset = 4*(nel1-1)+nod1
              elseif (nod1.eq.2) then
                kset = neps1 + 1
              endif

c             Form stiffness and residual

              istgn = nint(ch2(p1(1),kset))

              if (istgn .gt. 0) then

c               Call material law

                ns     = ix1(nod1,nel1)
                cxs(1) = x(1,ns) + u(1,ns)
                cxs(2) = x(2,ns) + u(2,ns)
                if(ndm.eq.3) then
                  cxs(3) = x(3,ns) + u(3,ns)
                else
                  cxs(3) = 0.0d0
                endif

                call criggap(cs02,ch2(1,kset),cn,cxs,ndm,resv,tanm, 2)

c               Get involved dof

                ixl(1) = ix1(nod1,nel1)

c               Assemble  stiffness and residual

                if(ifsolm.eq.1) then
                  ilm(1) = 0
                  n1     = 0
                elseif(ifsolm.eq.2) then
                  ilm(1) = ixl(1)
                  n1     = 1
                endif
                call constass (ixl,ida,1,ndm,ilm,n1,n1,4,tanm,resv)

c               Dump on file of stiffness matrix

                if (ifdb) then
                  write (99,*) 'ixl',ixl
                  write (99,4000) 'RESV', (resv(i),i=1,4)
                  write (99,4000) 'TANM',((tanm(j,i),i=1,4),j=1,4)
                endif
              endif

            endif
          end do ! kn
        end do ! ki

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from FORMFE to compute residual                    :   6

      elseif (csw.eq.6) then

c       Contact geometry
c       Loop over all surface 1 elements

        if(max(cp0(4,0),abs(cp0(5,0))).eq.0.0d0) then
          cn(1) = cp0(3,2)
          cn(2) = 0.0d0
        else
          cn(1) = cp0(4,0)
          cn(2) = cp0(5,0)
        endif
        do ke = 1,neps1
          nel1 = ke
          do kn = 1,nope1
            nod1 = kn

c           Skip if node just checked

            if (dnope1.eq.1 .or. nope1.eq.4 .or.
     &          (nod1.eq.1) .or. (ix1(dnope1,nel1).eq.0)) then

c             Search for closest master node (global search

              kset = nel1
              if (nope1.eq.4) then
                kset = 4*(nel1-1)+nod1
              elseif (nod1.eq.2) then
                kset = neps1 + 1
              endif

c             Form residual

              istgn = nint(ch2(p1(1),kset))
              if (istgn .gt. 0) then

                ns = ix1(nod1,nel1)
                cxs(1) = x(1,ns) + u(1,ns)
                cxs(2) = x(2,ns) + u(2,ns)
                if(ndm.eq.3) then
                  cxs(3) = x(3,ns) + u(3,ns)
                else
                  cxs(3) = 0.0d0
                endif

                call criggap(cs02,ch2(1,kset),cn,cxs,ndm,resv,tanm, 2)

c               Get involved dof

                ixl(1) = ix1(nod1,nel1)

c               Clean for security stiffness

                do j = 1,4
                  do i = 1,4
                    tanm(i,j) = 0.0d0
                  end do ! i
                end do ! j

c               Assemble  residual

                if(ifsolm.eq.1) then
                  ilm(1) = 0
                  n1     = 0
                elseif(ifsolm.eq.2) then
                  ilm(1) = ixl(1)
                  n1     = 1
                endif
                call constass (ixl,ida,1,ndm,ilm,n1,n1,4,tanm,resv)

              endif
            endif
          end do ! kn
        end do ! ke

c-----[--.----+----.----+----.-----------------------------------------]
c     Augmentation

      elseif (csw.eq.10) then

c       Check if augmentation is active

        if (ifaugm.ne.1) then

c         Perform update of variables after the last iteration
c         Loop to update geometrical variables and contact forces

          if(max(cp0(4,0),abs(cp0(5,0))).eq.0.0d0) then
            cn(1) = cp0(3,2)
            cn(2) = 0.0d0
          else
            cn(1) = cp0(4,0)
            cn(2) = cp0(5,0)
          endif
          do ke = 1,neps1
            nel1 = ke
            do kn = 1,nope1
              nod1 = kn

c             Skip if node just checked

              if (dnope1.eq.1 .or. nope1.eq.4 .or.
     &            (nod1.eq.1) .or. (ix1(dnope1,nel1).eq.0)) then

c               Set number

                kset = nel1
                if (nope1.eq.4) then
                  kset = 4*(nel1-1)+nod1
                elseif (nod1.eq.2) then
                  kset = neps1 + 1
                endif

c               Update augmented force

                istgn = nint(ch2(p1(1),kset))
                if (istgn .gt. 0) then
                  ch2(p1(151),kset) = ch2(p1(151),kset)
     &                              + cn(1)*ch2(p1(2),kset)
                  augg              = max(abs(ch2(p1(2),kset)),augg)
                endif

              endif

            end do ! kn
          end do ! ke

        endif

c-----[--.----+----.----+----.-----------------------------------------]
c     History variables initialization

      elseif (csw.eq.14) then

c       Loop over all surface 1 elements

        if(max(cp0(4,0),abs(cp0(5,0))).eq.0.0d0) then
          cn(1) = cp0(3,2)
          cn(2) = 0.0d0
        else
          cn(1) = cp0(4,0)
          cn(2) = cp0(5,0)
        endif
        do ke = 1,neps1
          nel1 = ke
          do kn = 1,nope1
            nod1 = kn

c           Skip if node just checked

            if (dnope1.eq.1 .or. nope1.eq.4 .or.
     &          (nod1.eq.1) .or. (ix1(dnope1,nel1).eq.0)) then

c             Coordinates of master node

c             ns = ix1(nod1,nel1)
c             cxs(1) = x(1,ns) + u(1,ns)
c             cxs(2) = x(2,ns) + u(2,ns)
c             if(ndm.eq.3) then
c               cxs(3) = x(3,ns) + u(3,ns)
c             else
c               cxs(3) = 0.0d0
c             endif

c             Search for closest master node (global search)

c             kset = nel1
c             if (nope1.eq.4) then
c               kset = 4*(nel1-1)+nod1
c             elseif (nod1.eq.2) then
c               kset = neps1 + 1
c             endif

c             Compute # of active elements
c             (for Lagrangian Multipliers, special augmentation, ...)

c             istgn = ch2(p1(1),kset)
c             if (istgn.gt.0) then
                nacte = nacte + 1
c             endif
            endif
          end do ! kn
        end do ! ke

c       Save total active elements

        cp0(11,-1) = nacte

c-----[--.----+----.----+----.-----------------------------------------]
c     CSW = 103: Called from PMACR1 to check geometry of contacts
c     CSW = 304: Called from PMACR3 to check geometry of contacts

      elseif (csw.eq.103 .or. csw.eq.304) then

c       Zero counter for active elements

        nacte = 0

c       Set flag for geometry check

        dopt  = nint(cp0(2,3)) ! Detection method
        lcnt  = lcnt + 1
        liter = nint(max(2.d0,cp0(3,3)))

        if(dopt.le.1) then
          compg = (csw.eq.103 .and.      ifistgn) .or.
     &            (csw.eq.304 .and. .not.ifistgn)
          lcnt  = 0
        elseif(dopt.gt.1 .and. lcnt.le.liter) then
          compg = .true.
        else
          compg = .false.
        endif
        if(max(cp0(4,0),abs(cp0(5,0))).eq.0.0d0) then
          cn(1) = cp0(3,2)
          cn(2) = 0.0d0
        else
          cn(1) = cp0(4,0)
          cn(2) = cp0(5,0)
        endif

c       Loop over all surface 1 elements

        if(compg) then
          do ke = 1,neps1
            nel1 = ke
            do kn = 1,nope1
              nod1 = kn

c             Skip if node just checked

              if (dnope1.eq.1 .or. nope1.eq.4 .or.
     &            (nod1.eq.1) .or. (ix1(dnope1,nel1).eq.0)) then

c               Search for closest master node (global search)

                kset = nel1
                if (nope1.eq.4) then
                  kset = 4*(nel1-1)+nod1
                elseif (nod1.eq.2) then
                  kset = neps1 + 1
                endif

                ns     = ix1(nod1,nel1)
                cxs(1) = x(1,ns) + u(1,ns)
                cxs(2) = x(2,ns) + u(2,ns)
                if(ndm.eq.3) then
                  cxs(3) = x(3,ns) + u(3,ns)
                else
                  cxs(3) = 0.0d0
                endif

c               Find gap to rigid segment

                call criggap(cs02,ch2(1,kset),cn,cxs,ndm,resv,tanm, 1)

c               Compute # of active elements
c               (for Lagrangian Multipliers, special augmentation, ...)

                istgn = nint(ch2(p1(1),kset))
                if (istgn.gt.0) then   ! contact
                  nacte = nacte + 1
                endif
              endif
            end do ! kn
          end do ! ke

c       else

          do ke = 1,neps1
            nel1 = ke
            do kn = 1,nope1
              nod1 = kn

c             Skip if node just checked

              if (dnope1.eq.1 .or. nope1.eq.4 .or.
     &            (nod1.eq.1) .or. (ix1(dnope1,nel1).eq.0)) then

                kset = nel1
                if (nope1.eq.4) then
                  kset = 4*(nel1-1)+nod1
                elseif (nod1.eq.2) then
                  kset = neps1 + 1
                endif

                istgn = nint(ch2(p1(1),kset))
                if (istgn.gt.0) then
                  nacte = nacte + 1
                endif
              endif
            end do ! kn
          end do ! ke

        endif

c       Save total active elements

        cp0(11,-1) = nacte

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from PMACR5 to show element informations

      elseif (csw.eq.200) then

        write (*,2001)

c-----[--.----+----.----+----.-----------------------------------------]
c     Printout of contact status

      elseif(csw.eq.204) then

c      Get printout flag and range

       call setcprt (ifprt,fel,lel)

c       Print title

        if (ifprt) then
          write (iow,2000) npair

c         Loop over requested surface 1 elements

          do ke = fel,lel
            nel1 = ke
            do kn = 1,nope1
              nod1 = kn

c             Skip if node just checked

              if (dnope1.eq.1 .or. nope1.eq.4 .or.
     &            (nod1.eq.1) .or. (ix1(dnope1,nel1).eq.0)) then
                kset = nel1
                if (nope1.eq.4) then
                  kset = 4*(nel1-1)+nod1
                elseif (nod1.eq.2) then
                  kset = neps1 + 1
                endif

                call printc1 (npair,nel1,nod1,ix1,ix2,ch2(1,kset))
              endif
            end do ! kn
          end do ! ke
        endif

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from PTIMPL to compute tplot values                : 206

      elseif (csw.eq.206) then

c       Contact geometry
c       Loop over all surface 1 elements

        if(ncplts.gt.0) then
          if(max(cp0(4,0),abs(cp0(5,0))).eq.0.0d0) then
            cn(1) = cp0(3,2)
            cn(2) = 0.0d0
          else
            cn(1) = cp0(4,0)
            cn(2) = cp0(5,0)
          endif
          do ke = 1,neps1
            nel1 = ke
            do kn = 1,nope1
              nod1 = kn

c             Skip if node just checked

              if (dnope1.eq.1 .or. nope1.eq.4 .or.
     &            (nod1.eq.1) .or. (ix1(dnope1,nel1).eq.0)) then

c               Search for closest master node (global search

                kset = nel1
                if (nope1.eq.4) then
                  kset = 4*(nel1-1)+nod1
                elseif (nod1.eq.2) then
                  kset = neps1 + 1
                endif

c               Form residual

                istgn = nint(ch2(p1(1),kset))
                if (istgn .gt. 0) then
                  ns = ix1(nod1,nel1)
                  cxs(1) = x(1,ns) + u(1,ns)
                  cxs(2) = x(2,ns) + u(2,ns)
                  if(ndm.eq.3) then
                    cxs(3) = x(3,ns) + u(3,ns)
                  else
                    cxs(3) = 0.0d0
                  endif

                  call crigplt(cs02,ch2(1,kset),cn(1),cxs)

                endif
              endif
            end do ! kn
          end do ! ke
        endif ! ncplts > 0

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from PPLOTF for plot of contact geometry

      elseif (csw.eq.305) then
        if(nint(cs02(1,0)).eq.6) then
          call c2rigplt (ix1,cs02(2,0),2,8)
        else
          call c2geoplt (ix1,ix2,2,8)
        endif

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from PPLOTF to set profile and range to plot variable
c     Called from CONTACT for plot contours of a contact variable

      elseif ((csw.eq.308) .or.
     &        (csw.eq.408)    ) then

        call c2varplt (ix1,ch1,ch2,ch3,npair,csw)

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from PCONTR for initialization of history variables

      elseif (csw.eq.313) then

c       Activate needed history variables

        call acthvptr()

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from UPDATE to update lagrange multiplier

      elseif (csw.eq.314) then

c       Loop over all surface 1 elements

        if (ifsolm.eq.2) then
          do ke = 1,neps1
            nel1 = ke
            do kn = 1,nope1
              nod1 = kn

c             Skip if node just checked

              if (dnope1.eq.1 .or. nope1.eq.4 .or.
     &            (nod1.eq.1) .or. (ix1(dnope1,nel1).eq.0)) then
                kset = nel1
                if (nope1.eq.4) then
                  kset = 4*(nel1-1)+nod1
                elseif (nod1.eq.2) then
                  kset = neps1 + 1
                endif

c               Update for lagrange multipliers

                if (nint(ch2(p1(1),kset)).gt.0) then
                  ilm(1) = ix1(nod1,nel1)
                  n1     = 1
                  call getlagm(ilm,n1,n1,ch2(p1(21),kset))
                else
                  ch2(p1(21),kset) = 0.0d0
                endif

              endif
            end do ! kn
          end do ! ke
        endif

c-----[--.----+----.----+----.-----------------------------------------]
c     Reset profile for contacts

      elseif (csw.eq.403) then

c       Loop over all surface 1 elements

        do ke = 1,neps1
          nel1 = ke
          do kn = 1,nope1
            nod1 = kn

c           Skip if node just checked

            if (dnope1.eq.1 .or. nope1.eq.4 .or.
     &          (nod1.eq.1) .or. (ix1(dnope1,nel1).eq.0)) then
              kset = nel1
              if (nope1.eq.4) then
                kset = 4*(nel1-1)+nod1
              elseif (nod1.eq.2) then
                kset = neps1 + 1
              endif

c             Reset profile for active contacts

              istgn = nint(ch2(p1(1),kset))
              if (istgn .gt. 0) then

                ixl(1) = ix1(nod1,nel1)

c               Penalty -> modify profile

                if (ifsolm.eq.1) then
                  call modprof (ixl,ida,1,ndm)

c               Lagrangian Multipliers -> compute equation numbers

                elseif (ifsolm.eq.2) then
                  ilm(1) = ixl(1)
                  call modprofl(ixl,ida,1,ndm,ilm,1,1)
                endif
              endif
            endif
          end do ! kn
        end do ! ke

c-----[--.----+----.----+----.-----------------------------------------]
c     Determine extra equations for special augmentation

      elseif (csw.eq.410) then

        lcnt = 0

c-----[--.----+----.----+----.-----------------------------------------]

      endif

2000  format (/'     C o n t a c t   O u t p u t   f o r   P a i r ',i5)

2001  format(//10x,'2-D: Node to Rigid Surface Contact')

3000  format (/
     &'  *** WARNING *** Unsymmetric solver required for friction'/)

4000  format( 5x,a/(1p,7e11.3))

      end
