c$Id:$
      subroutine cnts2d (ndm,ndf,x,u,csw,npair,cs02,cp0,
     &                   ix1,ix2,cm1,ch1,ch2,ch3,w1,w3)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Remove cp0 from call to cstif1, geopar1,         21/04/2007
c          geopar2 and penares
c       2. Add axisymmetric modification for penalty terms  09/07/2010
c       3. Remove 'resfl' from call to 'penares'            31/01/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: Contact DRIVER for 2D NTS

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
c         cm1(*)  - Contact materials data for surface 1

c      Outputs:
c         ch1(*)  - Contact history variables (old)
c         ch2(*)  - Contact history variables (current)
c         ch3(*)  - Contact history variables (static)
c         w1(*)   - Dictionary of variables for CH1 & CH2
c         w3(*)   - Dictionary of variables for CH3
c                 - Data exchange with main program subroutine calls
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_comnd.h'
      include  'c_contac.h'
      include  'c_geom.h'
      include  'c_keyh.h'
      include  'c_mate.h'
      include  'c_pair.h'
      include  'c_tole.h'
      include  'cdata.h'
      include  'compac.h'
      include  'compas.h'
      include  'counts.h'
      include  'eqsym.h'
      include  'iofile.h'
      include  'pointer.h'
      include  'print.h'
      include  'tdata.h'
      include  'comblk.h'

      logical   ifprt,once, compg
      character w1(*)*(*),w3(*)*(*)
      integer   ndm,ndf, csw,npair,ix1(dnope1,*),ix2(dnope2,*)
      integer   ke,kn,nel1,nod1,ns,kset,istgn,fel,lel,masts
      integer   iffricor,ixl(3),ida(2),ilm(2)
      integer   istgt,istfr,n1,n2,i,j
      integer   dopt, lcnt,liter
      real*8    cs02(nr0,n0c1:*)
      real*8    cp0(nr0,n0c3:*),cm1(*), cn,alpha
      real*8    ch1(lh1,*),ch2(lh1,*),ch3(lh3,*), x(ndm,*),u(ndf,*)
      real*8    cxs(2),tanm(8,8),resv(8),gn,gt,dtd, rnorm,nnorm

      real*8    fn, fntot, fttot, fratio

      save

c     Set active dof and dof order (idl length = ndf)

      data      ida     /1,2/
      data      ilm     /0,0/

      call cdebug0 ('    cnts2d',csw)

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from PCONTR for activation of requested history variables

      if(csw.eq.0) then
        once = .true.

      elseif (csw.eq.1) then

c       Once actions

        if (once) then
          once = .false.

c         Load dictionary of history variables

          call defhvar1 (w1,w3)

c         Print a warning if unsymmetric solver is needed

          if (iffric.eq.1) then
            write (*,3000)
          endif
        endif

c-----[--.----+----.----+----.-----------------------------------------]
c       Initial check automatically carried out with csw = 14

      elseif (csw.eq.2) then
        continue

c-----[--.----+----.----+----.-----------------------------------------]
c       Called from FORMFE to compute stiffness and residual

      elseif (csw.eq.3) then

        if (ifdb) then
          write (*,4000)
        endif

c       Loop over all surface 1 elements

        do ke = 1,neps1
          nel1 = ke
          do kn = 1,nope1
            nod1 = kn

c           Skip if node just checked

            if (dnope1.eq.1 .or.
     &          (nod1.eq.1) .or. (ix1(dnope1,nel1).eq.0)) then
              kset = nel1
              if (nod1.eq.2) then
                kset = neps1 + 1
              endif

c             Form stiffness and residual

              istgn = nint(ch2(p1(4),kset))

              if (istgn .ge. 0) then

c               Call material law

                if (ifdb) then
                  write(99,*) ' IFMTY1 = ',ifmty1
                endif
                if    (ifmty1.le.1) then
                  call cmatl1 (cp0,cm1,ch1(1,kset),ch2(1,kset),
     &                         ch3(1,kset))
                elseif(ifmty1.eq.2) then
                  call cmatl2 (cp0,cm1,ch1(1,kset),ch2(1,kset),
     &                         ch3(1,kset))
                elseif(ifmty1.eq.3) then
                  call cumatl (cp0,cm1,ch1(1,kset),ch2(1,kset),
     &                         ch3(1,kset))
                endif

c               Form stiffness and residual from material law

                call cstif1 (ch2(1,kset),ch3(1,kset),8,tanm,resv)

c               Get involved dof

                ixl(1) = ix1(nod1,nel1)
                masts  = nint(ch2(p1(1),kset))
                ixl(2) = ix2(1,masts)
                ixl(3) = ix2(2,masts)

c               Check for axisymmetric contact

                if(cp0(3,9).gt.0.0d0) then
                  cxs(1) = x(1,ixl(1))      ! radius of slave node
                  if(cxs(1).eq.0.0d0) then
                    cxs(1) = 1.0d-02*(x(1,ixl(2)) + x(1,ixl(3)))
                  endif
                  do i = 1,8
                    resv(i) = resv(i)*cxs(1)
                    do j = 1,8
                      tanm(j,i) = tanm(j,i)*cxs(1)
                    end do ! j
                  end do ! i
                endif

c               Assemble  stiffness and residual

                if(ifsolm.eq.1) then
                  ilm(1) = 0
                  n1     = 0
                elseif(ifsolm.eq.2) then
                  ilm(1) = ixl(1)
                  n1     = 1
                endif
                call constass (ixl,ida,3,2,ilm,n1,n1,8,tanm,resv)

c               Dump on file of stiffness matrix

                if (ifdb) then
                  write (99,*) 'ixl',ixl
                  write (99,4002) 'RESV-Ntos', (resv(i),i=1,7)
                  write (99,4002) 'TANM-Ntos',((tanm(j,i),i=1,7),j=1,7)
                endif
              endif

c             Screen printout for debugging

              if (ifdb) then
                istgt   = nint(ch2(p1(3),kset))
                gn      = ch2(p1(9),kset)
                gt      = ch2(p1(10),kset)
                if (iffric.eq.1) then
                  dtd   = ch2(p1(12),kset)
                  istfr = nint(ch2(p1(58),kset))
                else
                  dtd   = 0.d0
                  istfr = 0
                endif
                masts   = nint(ch2(p1(1),kset))
                ns      = ix1(nod1,nel1)
                n1      = ix2(1,masts)
                n2      = ix2(2,masts)
                write (iow,4001) gn,istgn,istgt,istfr,dtd,ns,n1,n2
                write (  *,4001) gn,istgn,istgt,istfr,dtd,ns,n1,n2
              endif
            endif
          end do
        end do

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from FORMFE to compute residual                    :   6
c     Called from PTIMPL to compute residual                    : 206

      elseif (csw.eq.6 .or. csw.eq.206) then

c       Contact geometry
c       Loop over all surface 1 elements

        do ke = 1,neps1
          nel1 = ke
          do kn = 1,nope1
            nod1 = kn

c           Skip if node just checked

            if (dnope1.eq.1 .or.
     &          (nod1.eq.1) .or. (ix1(dnope1,nel1).eq.0)) then
              ns     = ix1(nod1,nel1)
              cxs(1) = x(1,ns) + u(1,ns)
              cxs(2) = x(2,ns) + u(2,ns)

c             Search for closest master node (global search

              kset = nel1
              if (nod1.eq.2) then
                kset = neps1 + 1
              endif

              call gloscln (ndm,ndf,x,u,ix2,ns,cxs,ch2(1,kset))

c             Find master segment

              call mastseg (ndm,ndf,x,u,ix2,ns,cxs,ch2(1,kset))

c             Compute geometrical parameters

              call geopar1 (ndm,ndf,x,u,nel1,nod1,ix1,ix2,cxs,
     &                      ns,mr(np(31)+ndf*numnp),
     &                      ch1(1,kset),ch2(1,kset),ch3(1,kset))

              istgn = nint(ch2(p1(4),kset))

              if (istgn .ge. 0) then

c               Call material law

                if    (ifmty1.le.1) then
                  call cmatl1 (cp0,cm1,ch1(1,kset),ch2(1,kset),
     &                         ch3(1,kset))
                elseif(ifmty1.eq.2) then
                  call cmatl2 (cp0,cm1,ch1(1,kset),ch2(1,kset),
     &                         ch3(1,kset))
                elseif(ifmty1.eq.3) then
                  call cumatl (cp0,cm1,ch1(1,kset),ch2(1,kset),
     &                         ch3(1,kset))
                endif

c               Form stiffness and residual

                call penares (ch2(1,kset),ch3(1,kset),resv)

c               Get involved dof

                ixl(1) = ix1(nod1,nel1)
                masts  = nint(ch2(p1(1),kset))
                ixl(2) = ix2(1,masts)
                ixl(3) = ix2(2,masts)

c               Clean for security stiffness

                do j = 1,8
                  do i = 1,8
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
                call constass (ixl,ida,3,2,ilm,n1,n1,8,tanm,resv)

              endif
            endif
          end do
        end do

c-----[--.----+----.----+----.-----------------------------------------]
c     Augmentation

      elseif (csw.eq.10) then

c       Check if augmentation is active

        if (ifaugm.ne.1) then

c         Perform update of variables after the last iteration
c         Loop to update geometrical variables and contact forces

          include 'c_loop01.h'
            ns     = ix1(nod1,nel1)
            cxs(1) = x(1,ns) + u(1,ns)
            cxs(2) = x(2,ns) + u(2,ns)

            call gloscln (ndm,ndf,x,u,ix2,ns,cxs,ch2(1,kset))

            call mastseg (ndm,ndf,x,u,ix2,ns,cxs,ch2(1,kset))

            call geopar1 (ndm,ndf,x,u,nel1,nod1,ix1,ix2,cxs,
     &                    ns,mr(np(31)+ndf*numnp),
     &                    ch1(1,kset),ch2(1,kset),ch3(1,kset))

            if    (ifmty1.le.1) then
              call cmatl1 (cp0,cm1,ch1(1,kset),ch2(1,kset),
     &                     ch3(1,kset))
            elseif(ifmty1.eq.2) then
              call cmatl2 (cp0,cm1,ch1(1,kset),ch2(1,kset),
     &                     ch3(1,kset))
            elseif(ifmty1.eq.3) then
              call cumatl (cp0,cm1,ch1(1,kset),ch2(1,kset),
     &                     ch3(1,kset))
            endif

          include 'c_end01.h'

c         Set to zero the residual norm

          rnorm = 0.d0

c         Perform augmentation with different schemes
c         Simple linear augmentation

          if (ifaugm.eq.2) then
            cn    = cp0(3,2)
            alpha = cp0(3,5)
            if (alpha.eq.0.d0) then
              alpha = 1.d0
            endif

            include 'c_loop01.h'
              call aug2 (alpha,ch2(1,kset),rnorm)
            include 'c_end01.h'
          endif

c         Residual norms

          rnorm = sqrt(rnorm)
          if (naugm.eq.0) then
            if (rnorm.ne.0.d0) then
              nnorm = rnorm
            else
              nnorm = 1.d0
            endif
          else if (nnorm.eq.1.d0) then
            if (rnorm.ne.0.d0) then
              nnorm = rnorm
            endif
          endif

        endif

c-----[--.----+----.----+----.-----------------------------------------]
c     History variables initialization

      elseif (csw.eq.14) then

c       Loop over all surface 1 elements

        do ke = 1,neps1
          nel1 = ke
          do kn = 1,nope1
            nod1 = kn

c           Skip if node just checked

            if (dnope1.eq.1 .or.
     &          (nod1.eq.1) .or. (ix1(dnope1,nel1).eq.0)) then

c             Coordinates of master node

              ns     = ix1(nod1,nel1)
              cxs(1) = x(1,ns) + u(1,ns)
              cxs(2) = x(2,ns) + u(2,ns)

c             Search for closest master node (global search)

              kset = nel1
              if (nod1.eq.2) then
                kset = neps1 + 1
              endif

              call gloscln (ndm,ndf,x,u,ix2,ns,cxs,ch2(1,kset))

c             Find master segment

              call mastseg (ndm,ndf,x,u,ix2,ns,cxs,ch2(1,kset))

c             Check for initial penetrations automatically performed

              if(ifpck) then
                call ipenck1 (ndm,x,ns,ix2,cxs,ch2(1,kset))
              endif

c             Determine geometry variables for initialization at t=0
c             skipping tangential disp variables computation

              iffricor = iffric
              if (iffric.eq.1) then
                iffric = 2
              endif
              call geopar1 (ndm,ndf,x,u,nel1,nod1,ix1,ix2,cxs,
     &                      ns,mr(np(31)+ndf*numnp),
     &                      ch1(1,kset),ch2(1,kset),ch3(1,kset))

c             Reset original tg flag

              iffric = iffricor

c             Set history variables at time T=0

              call inchv1 (ch1(1,kset),ch2(1,kset))

            endif
          end do
        end do

c-----[--.----+----.----+----.-----------------------------------------]
c     CSW = 103: Called from PMACR1 to check geometry of contacts
c     CSW = 304: Called from PMACR3 to check geometry of contacts

      elseif (csw.eq.103 .or. csw.eq.304) then

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

c       Loop over all surface 1 elements

        if(compg) then
          do ke = 1,neps1
            nel1 = ke
            do kn = 1,nope1
              nod1 = kn

c             Skip if node just checked

              if (dnope1.eq.1 .or.
     &            (nod1.eq.1) .or. (ix1(dnope1,nel1).eq.0)) then
                ns     = ix1(nod1,nel1)
                cxs(1) = x(1,ns) + u(1,ns)
                cxs(2) = x(2,ns) + u(2,ns)

c               Search for closest master node (global search)

                kset = nel1
                if (nod1.eq.2) then
                  kset = neps1 + 1
                endif

c               Find master segment

                call gloscln (ndm,ndf,x,u,ix2,ns,cxs,ch2(1,kset))

                call mastseg (ndm,ndf,x,u,ix2,ns,cxs,ch2(1,kset))

c               Compute geometrical parameters

                call geopar1 (ndm,ndf,x,u,nel1,nod1,ix1,ix2,cxs,
     &                        ns,mr(np(31)+ndf*numnp),
     &                        ch1(1,kset),ch2(1,kset),ch3(1,kset))
              endif
            end do
          end do

        else

          do ke = 1,neps1
            nel1 = ke
            do kn = 1,nope1
              nod1 = kn

c             Skip if node just checked

              if (dnope1.eq.1 .or.
     &            (nod1.eq.1) .or. (ix1(dnope1,nel1).eq.0)) then
                ns     = ix1(nod1,nel1)
                cxs(1) = x(1,ns) + u(1,ns)
                cxs(2) = x(2,ns) + u(2,ns)

c               Search for closest master node (global search

                kset = nel1
                if (nod1.eq.2) then
                  kset = neps1 + 1
                endif

                istgn = nint(ch2(p1(4),kset))
                if (istgn.ge.0) then

                  call geopar2 (ndm,ndf,x,u,nel1,nod1,ix1,ix2,cxs,
     &                          ns,mr(np(31)+ndf*numnp),
     &                          ch1(1,kset),ch2(1,kset),ch3(1,kset))
                endif
              endif
            end do
          end do

        endif

c-----[--.----+----.----+----.-----------------------------------------]
c     Called from PMACR5 to show element informations

      elseif (csw.eq.200) then

        write (*,2001)

c-----[--.----+----.----+----.-----------------------------------------]
c     Printout of contact status

      elseif(csw.eq.204) then

c       Get printout flag and range

        call setcprt (ifprt,fel,lel)

c       Print title

        if (ifprt) then
          write (iow,2000) npair
          do ke = fel,lel
            nel1 = ke
            do kn = 1,nope1
              nod1 = kn

c             Skip if node just checked

              if (dnope1.eq.1 .or.
     &            (nod1.eq.1) .or. (ix1(dnope1,nel1).eq.0)) then
                kset = nel1
                if (nod1.eq.2) then
                  kset = neps1 + 1
                endif

                call printc1 (npair,nel1,nod1,ix1,ix2,ch2(1,kset))

              endif
            end do ! kn
          end do ! ke

        endif ! ifprt

c       Loop over requested surface 1 elements

        fntot = 0.0d0
        fttot = 0.0d0
        do ke = 1,neps1
          nel1 = ke
          do kn = 1,nope1
            nod1 = kn

c           Skip if node just checked

            if (dnope1.eq.1 .or.
     &          (nod1.eq.1) .or. (ix1(dnope1,nel1).eq.0)) then
              kset = nel1
              if (nod1.eq.2) then
                kset = neps1 + 1
              endif

c             Extract and accumlate normal/tangential forces

              fn = 0.0d0
              if(ifsolm.eq.2) then
                fn = ch2(p1(21),kset)
              elseif(ifaugm.ne.1) then
                fn = ch2(p1(51),kset) + ch3(p1(6),kset)
              else
                fn = ch2(p1(51),kset)
              endif
              if(ch2(p1(9),kset).gt.0.0001d0) then
                fn = 0.0d0
              endif
              if(fn.ne.0.0d0) then
                fntot = fntot + fn
                if(p1(53).ne.0) then
                  fttot = fttot + ch2(p1(53),kset)
                endif
              endif
            endif
          end do ! kn
        end do ! ke

c       Output total normal/tangential forces

        if(fntot.ne.0.0d0) then
          fratio = fttot/fntot
        else
          fratio = 0.0d0
        endif
        write(iow,2002) ' FNTOT: ',ttim,fntot,fttot,fratio

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

        call acthvar1

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

              if (dnope1.eq.1 .or.
     &            (nod1.eq.1) .or. (ix1(dnope1,nel1).eq.0)) then
                kset = nel1
                if (nod1.eq.2) then
                  kset = neps1 + 1
                endif

c               Update for Lagrange multipliers

                if (nint(ch2(p1(4),kset)).ge.0) then
                  ilm(1) = ix1(nod1,nel1)
                  n1     = 1
                  call getlagm(ilm,n1,n1,ch2(p1(21),kset))
                else
                  ch2(p1(21),kset) = 0.0d0
                endif

              endif
            end do
          end do
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

            if (dnope1.eq.1 .or.
     &          (nod1.eq.1) .or. (ix1(dnope1,nel1).eq.0)) then
              kset = nel1
              if (nod1.eq.2) then
                kset = neps1 + 1
              endif

c             Reset profile for active contacts

              istgn = nint(ch2(p1(4),kset))
              if (istgn .ge. 0) then

                ixl(1) = ix1(nod1,nel1)
                masts  = nint(ch2(p1(1),kset))
                ixl(2) = ix2(1,masts)
                ixl(3) = ix2(2,masts)

c               Penalty -> modify profile

                if (ifsolm.eq.1) then
                  call modprof (ixl,ida,3,2)

c               Lagrangian Multipliers -> compute equation numbers

                elseif (ifsolm.eq.2) then
                  ilm(1) = ixl(1)
                  call modprofl(ixl,ida,3,2,ilm,1,1)
                endif
              endif
            endif
          end do
        end do

c-----[--.----+----.----+----.-----------------------------------------]
c     Determine extra equations for special augmentation

      elseif (csw.eq.410) then

        lcnt = 0

c-----[--.----+----.----+----.-----------------------------------------]

      endif

2000  format (/'     C o n t a c t   O u t p u t   f o r   P a i r ',i5)

2001  format(//10x,'2-D: Node to Surface Contact')

2002  format(a,1p,4e13.5)
3000  format (/
     &'  *** WARNING *** Unsymmetric solver required for friction'/)

4000  format(/22x,'gn istgn istgt istfr',21x,'dtd  ns  n1  n2')
4001  format (e24.15,3i6,e24.15,3i4)
4002  format( 5x,a/(1p,7e11.3))

      end
