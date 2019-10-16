c$Id:$
      subroutine cnts3d (ndm,ndf,x,u,csw,npair,cp0,
     &                   ix1,ix2,cm1,ch1,ch2,ch3,w1,w3)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Anna Haraldsson             February 1998            1.0

c      Acronym: Contact DRIVER for 3D NTS

c      Purpose: Driver 'nts' for 3d contact (node-to-segment-contact)
c               sucht nachbar zu kante entsprechend Knoten
c               surf1 = quad = slave
c               surf2 = quad = master

c      Inputs:
c         ndm     - Space dimension of mesh
c         ndf     - Number dof/node
c         x(*)    - Nodal coordinates
c         u(*)    - Current nodal solution vectors
c         csw     - Contact switch
c         npair   - # of current pair
c         cp0(*)  - Contact pair control data
c         ix1(*)  - Element nodal connection list for surface 1
c         ix2(*)  - Element nodal connection list for surface 2
c         cm1(*)  - Contact materials data storage for surface 1

c      Outputs:
c         ch1(*)  - Contact history variables (old)
c         ch2(*)  - Contact history variables (current)
c         ch3(*)  - Contact history variables (static)
c         w1(*)   - Dictionary of variables for CH1 & CH2
c                 - Data exchange with main program subroutine calls
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'c_0.h'
      include   'c_comnd.h'
      include   'c_contac.h'
      include   'c_geom.h'
      include   'c_keyh.h'
      include   'c_mate.h'
      include   'c_pair.h'
      include   'c_tole.h'
      include   'augdat.h'
      include   'compas.h'
      include   'iofile.h'
      include   'ndata.h'
      include   'print.h'
      include   'pointer.h'
      include   'p_int.h'
      include   'comblk.h'

      logical    ifprt,errck,active,flag,test, setval,palloc
      logical    change(22),chngto(22),newflag(22),reibg(22)
      character  w1(*)*(*),w3(*)*(*)
      integer    ndm,ndf, csw,npair,ix1(dnope1,*),ix2(dnope2,*)
      integer    i,j,ke,nod2,ns,kset,istgn,fel,lel,masts,maxseg
      integer    ixl(5),ida(3),ilm(3),istgi(0:4),im(4)
      real*8     x(ndm,*),u(ndf,*), cp0(nr0,n0c3:*),cm1(*)
      real*8     ch1(lh1,*),ch2(lh1,*),ch3(lh3,*),xs(3),xm(3,4)
      real*8     mue ,mreibg,mreibg2,tanm(18,18),resv(18)

      save

c     Set active dof and dof order (idl length = ndf)

      data      ida     /1,2,3/
      data      ilm     /0,0,0/

      newflag(npair)=.false.

      if(npair.ge.22) then
        write(*,*)' *ERROR* Increase size of array reibg to ',npair
        call plstop()
      endif

      call cdebug0 ('    cnts3d',csw)

      if(csw.ne.1      .and.    csw.ne.12 .and.
     &   change(npair) .and. reibg(npair)) then

        if(ifdb .and. indb.ge.2) then
          if(iffric.eq.1) then
            write(*,*) 'friction was on for contact pair',npair,
     &                 ', mue=',cm1(1)
          else
            write(*,*)'friction was off for contact pair',npair
          endif
        endif

        if (chngto(npair)) then
          if(iffric.eq.1) then
            cm1(1) =  abs(cm1(1))
            if(ifdb .and. indb.ge.2) then
              write(*,*) 'and changed to on for contact pair',npair,
     &                   ', mue=',cm1(1)
            endif
          endif
          newflag(npair) = .true.
        else
          if(iffric.eq.1) then
            cm1(1) = -abs(cm1(1))
            if(ifdb .and. indb.ge.2) then
              write(*,*) 'and changed to off for contact pair',npair
            endif
          endif
        endif

        change(npair) = .false.

      endif

c-----[--.----+----.----+----.-----------------------------------------]
cc--> csw <--> 001 <--> csw <--> 001 <--> csw <--> 001 <--> csw <-->
c  Called from PCONTR for activation of requested history variables
c-----[--.----+----.----+----.-----------------------------------------]

      if (csw.eq.0) then

      elseif (csw.eq.1) then

        mue = cm1(1)

c       Load dictionary of history variables

c       masts   - MASTter Segment
c       istgt   - Index of STatus for GT
c       istgn   - Index of STatus for GN
c       gn      - Normal projection Gap  (+ if open)
c       nvec    - normal vector to master segment
c       xi      - ...
c       a1      - vector tangential to surface
c       a2      - vector tangential to surface
c       a11,a12,a22 - metric tensor
c       diffxi1 - difference of xi(1) in timestep to previous one
c       diffxi2 - difference of xi(2) ...

c       var. and lin:    slave
c                        1.master
c                        2.master
c                        3.master
c                        4.master

c       CH1 & CH2 VARIABLES (CH2 copied in CH1 at every time step)

        w1( 1) = 'masts'
        w1( 3) = 'istgt'
        w1( 4) = 'istgn'
        w1(10) = 'diffxi1'
        w1(11) = 'diffxi2'
        w1(19) = 'a1'
        w1(21) = 'a11'
        w1(22) = 'a12'
        w1(23) = 'a22'
        w1(24) = 'xi'
        w1(25) = 'nachbar'
        w1(26) = 'gt'

        w1(31) = 'lagmu'
        w1(32) = 'lagmv'

c       CH3 VARIABLES (CH3 never copied)

        w3( 2) = 'knflg'
        w3( 5) = 'ft3'
        w3( 9) = 'gn'
        w3(13) = 'faug'
        w3(14) = 'fn'
        w3(15) = 'ft'
        w3(16) = 'flaeche'
        w3(18) = 'nvec'
        w3(20) = 'a2'

        w3(31) = 'lagmn'
        w3(32) = 'lagmt'

        change(npair) = .false.
        reibg(npair)  = .false.
        if(iffric.eq.1) then
          reibg(npair) = .true.
          if(ifdb .and. indb.ge.2) then
            write(*,*) 'Contact pair',npair,' has friction.'
          endif
        endif

c-----[--.----+----.----+----.-----------------------------------------]
cc--> csw <--> 002 <--> csw <--> 002 <--> csw <--> 002 <--> csw <-->
c       Initial check automatically carried out with csw = 14
c-----[--.----+----.----+----.-----------------------------------------]

      elseif (csw.eq.2) then
        continue

c-----[--.----+----.----+----.-----------------------------------------]
cc--> csw <--> 003 <--> csw <--> 003 <--> csw <--> 003 <--> csw <-->
c       Called from FORMFE to compute stiffness and residual
c-----[--.----+----.----+----.-----------------------------------------]

      elseif (csw.eq.3) then

c       Penalty for normal stiffness

        fp(1)  = np(191) + mr(np(192)+nsurf1-1)

c       Loop over all surface 1 nodes = Loop over all slave nodes

        do ke = 1,mr(fp(1))
          kset   = ke
          ns     = mr(fp(1)+kset)
          xs(1) = x(1,ns) + u(1,ns)
          xs(2) = x(2,ns) + u(2,ns)
          xs(3) = x(3,ns) + u(3,ns)

c         Compute and store all geometrical parameters

          if(iffric.eq.1) then
            mue  = cm1(1)
          else
            mue  = 0.0d0
          endif
          test = .false.
          test = .true.
          call gnqtq (kset,x,u,ix1,ix2,xs,ch1(1,kset),ch2(1,kset),
     &                ch3(1,kset),ida,im,mr(np(192)),mue,test)

          istgn = nint(ch2(p1(4),kset))

          if (test) then

c           Initialize tangent and residual

            do i = 1,18
              resv(i) = 0.0d0
              do j = 1,18
                tanm(j,i) = 0.0d0
              end do ! j
            end do ! i

c           No contact case

            if(istgn.eq.0) then

c           Contact with facet

            elseif(istgn.eq.1) then

c             Get involved dof

              ixl(1) = ns
              masts  = nint(ch2(p1(1),kset))
              ixl(2) = ix2(1,masts)
              ixl(3) = ix2(2,masts)
              ixl(4) = ix2(3,masts)
              ixl(5) = ix2(4,masts)

              do nod2 = 1,4
                do j=1,3
                  xm(j,nod2) = x(j,ixl(nod2+1)) + u(j,ixl(nod2+1))
                end do ! j
              end do ! nod2

c             Form stiffness and residual

              call stfqtq1 (cp0,cm1,ch2(1,kset),ch1(1,kset),
     &                      ch3(1,kset),tanm,resv,xm)

c             Assemble stiffness and residual

              if(ifsolm.eq.1) then
                ilm(1) = 0
                call constass (ixl,ida,5,3,ilm,0,0,18,tanm,resv)
              elseif(ifsolm.eq.2) then
                ilm(1) = ixl(1)
                call constass (ixl,ida,5,3,ilm,1,1,18,tanm,resv)
              endif

            elseif((istgn.eq.2).or.(istgn.eq.3)) then

c             Get involved dof

              ixl(1) = ns
              ixl(2) = nint(ch2(p1( 1),kset))
              ixl(3) = nint(ch2(p1(25),kset))

              xs(1)  = x(1,ixl(1)) + u(1,ixl(1))
              xs(2)  = x(2,ixl(1)) + u(2,ixl(1))
              xs(3)  = x(3,ixl(1)) + u(3,ixl(1))
              do nod2 = 1,2
                do j = 1,3
                  xm(j,nod2) = x(j,ixl(nod2+1)) + u(j,ixl(nod2+1))
                end do ! j
              end do ! nod2

c             Form stiffness and residual

c             if(iffric.eq.1) then
c               mue  = cm1(1)
c             else
c               mue  = 0.0d0
c             endif
              call stfqtq2 (cp0,cm1,ch2(1,kset),ch1(1,kset),
     &                      ch3(1,kset),tanm,resv,xm)

c             Assemble stiffness and residual

              if(ifsolm.eq.1) then
                ilm(1) = 0
                call constass (ixl,ida,3,3,ilm,0,0,18,tanm,resv)
              elseif(ifsolm.eq.2) then
                ilm(1) = ixl(1)
                call constass (ixl,ida,3,3,ilm,1,1,18,tanm,resv)
              endif

            elseif(istgn.eq.4) then

c             Get involved dof

              ixl(1) = ns
              ixl(2) = nint(ch2(p1(1),kset))

              xs(1)  = x(1,ixl(1)) + u(1,ixl(1))
              xs(2)  = x(2,ixl(1)) + u(2,ixl(1))
              xs(3)  = x(3,ixl(1)) + u(3,ixl(1))
              do j = 1,3
                 xm(j,1) = x(j,ixl(2)) + u(j,ixl(2))
              end do ! j

c             Form stiffness and residual

c             if(iffric.eq.1) then
c               mue  = cm1(1)
c             else
c               mue  = 0.0d0
c             endif
              call stfqtq4 (cp0,cm1,ch2(1,kset),ch1(1,kset),
     &                      ch3(1,kset),tanm,resv)

c             Assemble  stiffness and residual

              if(ifsolm.eq.1) then
                ilm(1) = 0
                call constass (ixl,ida,2,3,ilm,0,0,18,tanm,resv)
              elseif(ifsolm.eq.2) then
                ilm(1) = ixl(1)
                call constass (ixl,ida,2,3,ilm,1,1,18,tanm,resv)
              endif

            endif

            else
          endif
        end do

c-----[--.----+----.----+----.-----------------------------------------]
cc--> csw <-->   6 <--> csw <-->   6 <--> csw <-->   6 <--> csw <-->
c     Called from FORMFE to compute residual
cc--> csw <--> 206 <--> csw <--> 206 <--> csw <--> 206 <--> csw <-->
c     Called from PTIMPL to compute residual
c-----[--.----+----.----+----.-----------------------------------------]

      elseif (csw.eq.6 .or. csw.eq.206) then

c       Contact geometry
c       Loop over all surface 1 elements = Loop over all slave nodes

        fp(1)  = np(191) + mr(np(192)+nsurf1-1)
        do ke = 1,mr(fp(1))

          kset   = ke
          ns     = mr(fp(1)+kset)
          xs(1) = x(1,ns) + u(1,ns)
          xs(2) = x(2,ns) + u(2,ns)
          xs(3) = x(3,ns) + u(3,ns)

c         Compute and store all geometrical parameters

          if(iffric.eq.1) then
            mue  = cm1(1)
          else
            mue  = 0.0d0
          endif

          test = .false.
          test = .true.
          call gnqtq (kset,x,u,ix1,ix2,xs,ch1(1,kset),ch2(1,kset),
     &                ch3(1,kset),ida,im,mr(np(192)),mue,test)

          istgn = nint(ch2(p1(4),kset))

          if (test) then

c           Initialize tangent and residual

            do i = 1,18
              resv(i) = 0.0d0
              do j = 1,18
                tanm(j,i) = 0.0d0
              end do ! j
            end do ! i

c           No contact case

            if(istgn.eq.0) then

c           Contact with facet

            elseif(istgn.eq.1) then

c             Get involved dof

              ixl(1) = mr(fp(1)+kset)
              masts  = nint(ch2(p1(1),kset))
              ixl(2) = ix2(1,masts)
              ixl(3) = ix2(2,masts)
              ixl(4) = ix2(3,masts)
              ixl(5) = ix2(4,masts)

              do nod2 = 1,4
                do j = 1,3
                  xm(j,nod2) = x(j,ixl(nod2+1)) + u(j,ixl(nod2+1))
                end do ! j
              end do ! nod2

c             Form stiffness and residual

              call stfqtq1 (cp0,cm1,ch2(1,kset),ch1(1,kset),
     &                      ch3(1,kset),tanm,resv,xm)

c             Clean for security stiffness

              do i = 1,18
                do j = 1,18
                  tanm(j,i) = 0.0d0
                end do ! j
              end do ! i

c             Assemble  residual

              if(ifsolm.eq.1) then
                ilm(1) = 0
                call constass (ixl,ida,5,3,ilm,0,0,18,tanm,resv)
              elseif(ifsolm.eq.2) then
                ilm(1) = ixl(1)
                call constass (ixl,ida,5,3,ilm,1,1,18,tanm,resv)
              endif

            elseif((istgn.eq.2).or.(istgn.eq.3)) then

c             Get involved dof

              ixl(1) = mr(np(191) + mr(np(192)+nsurf1-1)+kset)
              ixl(2) = nint(ch2(p1(1),kset))
              ixl(3) = nint(ch2(p1(25),kset))

              xs(1)  = x(1,ixl(1)) + u(1,ixl(1))
              xs(2)  = x(2,ixl(1)) + u(2,ixl(1))
              xs(3)  = x(3,ixl(1)) + u(3,ixl(1))
              do nod2 = 1,2
                do j = 1,3
                  xm(j,nod2) = x(j,ixl(nod2+1)) + u(j,ixl(nod2+1))
                end do ! j
              end do ! nod2

c             Form stiffness and residual

c             if(iffric.eq.1) then
c               mue  = cm1(1)
c             else
c               mue  = 0.0d0
c             endif
              call stfqtq2 (cp0,cm1,ch2(1,kset),ch1(1,kset),
     &                      ch3(1,kset),tanm,resv,xm)

c             Clean for security stiffness

              do i = 1,18
                do j = 1,18
                  tanm(j,i) = 0.0d0
                end do ! j
              end do ! i

c             Assemble  stiffness and residual

              if(ifsolm.eq.1) then
                ilm(1) = 0
                call constass (ixl,ida,3,3,ilm,0,0,18,tanm,resv)
              elseif(ifsolm.eq.2) then
                ilm(1) = ixl(1)
                call constass (ixl,ida,3,3,ilm,1,1,18,tanm,resv)
              endif

            elseif(istgn.eq.4) then

c             Get involved dof

              ixl(1) = mr(np(191) + mr(np(192)+nsurf1-1)+kset)
              ixl(2) = nint(ch2(p1(1),kset))

              xs(1)  = x(1,ixl(1)) + u(1,ixl(1))
              xs(2)  = x(2,ixl(1)) + u(2,ixl(1))
              xs(3)  = x(3,ixl(1)) + u(3,ixl(1))
              do j = 1,3
                xm(j,1) = x(j,ixl(2)) + u(j,ixl(2))
              end do ! j

c             Form stiffness and residual

c             if(iffric.eq.1) then
c               mue  = cm1(1)
c             else
c               mue  = 0.0d0
c             endif
              call stfqtq4 (cp0,cm1,ch2(1,kset),ch1(1,kset),
     &                      ch3(1,kset),tanm,resv)

c             Clean for security stiffness

              do i = 1,18
                do j = 1,18
                  tanm(j,i) = 0.0d0
                end do ! j
              end do ! i

c             Assemble  stiffness and residual

              if(ifsolm.eq.1) then
                ilm(1) = 0
                call constass (ixl,ida,2,3,ilm,0,0,18,tanm,resv)
              elseif(ifsolm.eq.2) then
                ilm(1) = ixl(1)
                call constass (ixl,ida,2,3,ilm,1,1,18,tanm,resv)
              endif

            endif
          endif
        end do

c-----[--.----+----.----+----.-----------------------------------------]
cc--> csw <--> 010 <--> csw <--> 010 <--> csw <--> 010 <--> csw <-->
c     Augmentation step
c-----[--.----+----.----+----.-----------------------------------------]

      elseif (csw.eq.10) then

        if(ifaugm.ge.2) then

c         Loop over all surface 1 elements = Loop over all slave nodes

          fp(1) = np(191)+mr(np(192)+nsurf1-1)
          do ke = 1,mr(fp(1))
            kset = ke
            if(nint(ch2(p1(4),kset)).gt.0) then
              ch3(p3(13),kset) = ch3(p3(13),kset)
     &                + cp0(3,2)*ch3(p3(16),kset)*ch3(p3(9),kset)
              augg = max(augg,abs(ch3(p3(9),kset)))
c             if(ch3(p3(13),kset).lt.0.0d0) then
c               ch3(p3(13),kset) = 0.0d0
c             endif
            endif
          end do ! ke

        endif


c-----[--.----+----.----+----.-----------------------------------------]
cc--> csw <--> 014 <--> csw <--> 014 <--> csw <--> 014 <--> csw <-->
c     History variables initialization
c-----[--.----+----.----+----.-----------------------------------------]

      elseif (csw.eq.14) then

        if (nope2.ne.4) then
          write(*,3000) nope2
        endif

c       Loop over all surface 1 elements = Loop over all slave nodes

        fp(1) = np(191)+mr(np(192)+nsurf1-1)
        fp(2) = np(191)+mr(np(192)+nsurf2-1)
        maxseg = 0
        do ke = 1,mr(fp(2))
          maxseg = max(maxseg,mr(np(191)+mr(fp(2)+mr(fp(2))+ke)))
        end do ! ke
        setval = palloc( 193,'INSEG', maxseg   ,1)
        setval = palloc( 194,'CNSEG', maxseg   ,1)
        setval = palloc( 195,'PNSEG', maxseg   ,1)
        setval = palloc( 196,'XISEG', maxseg*2 ,2)
        do ke = 1,mr(fp(1))
          kset = ke
          ns   = mr(fp(1)+kset)

c         Coordinates of slave node

          xs(1) = x(1,ns) + u(1,ns)
          xs(2) = x(2,ns) + u(2,ns)
          xs(3) = x(3,ns) + u(3,ns)

c         Compute and store all geometrical parameters
c         skipping tangential disp variables computation

          mue = -1.0d0
          call geoqtq (kset,x,u,ix1,ix2,ns,xs,ch1(1,kset),ch2(1,kset),
     &                 ch3(1,kset),   ida ,mr(np(192)),mr(np(193)),
     &                 mr(np(194)),mr(np(195)),hr(np(196)),mue,csw)

c         Set history variables at time T=0

          ch1(p1( 1)  ,kset) = ch2(p1( 1)  ,kset)     !masts
          ch1(p1( 4)  ,kset) = ch2(p1( 4)  ,kset)     !istgn
          ch1(p1(19)  ,kset) = ch2(p1(19)  ,kset)     !a1
          ch1(p1(19)+1,kset) = ch2(p1(19)+1,kset)
          ch1(p1(19)+2,kset) = ch2(p1(19)+2,kset)
          ch1(p1(24)  ,kset) = ch2(p1(24)  ,kset)     !xi
          ch1(p1(24)+1,kset) = ch2(p1(24)+1,kset)     !

          if(iffric.eq.1) then
            ch1(p1( 3)  ,kset) = 0.d0     !istgt
            ch1(p1(10)  ,kset) = 0.d0     !diffxi1
            ch1(p1(11)  ,kset) = 0.d0     !diffxi2
            ch1(p1(26)  ,kset) = 0.d0     !gt
            ch1(p1(26)+1,kset) = 0.d0
            ch1(p1(26)+2,kset) = 0.d0
          endif
        end do
        setval = palloc( 196,'XISEG', 0,2)
        setval = palloc( 195,'PNSEG', 0,1)
        setval = palloc( 194,'CNSEG', 0,1)
        setval = palloc( 193,'INSEG', 0,1)

c-----[--.----+----.----+----.-----------------------------------------]
cc--> csw <--> 021 <--> csw <--> 021 <--> csw <--> 021 <--> csw <-->

c     Output of special data for umacr RAUS

c-----[--.----+----.----+----.-----------------------------------------]

      elseif (csw.eq.21) then

c       Loop over all surface 1 elements = Loop over all slave nodes

        mreibg  = 0.0d0
        mreibg2 = 0.0d0
        fp(1)     = np(191)+mr(np(192)+nsurf1-1)
        do ke = 1,mr(fp(1))
          kset    = ke
          ns      = mr(fp(1)+kset)
          xs(1)  = x(1,ns) + u(1,ns)
          xs(2)  = x(2,ns) + u(2,ns)
          xs(3)  = x(3,ns) + u(3,ns)

c         Compute and store all geometrical parameters

          if(iffric.eq.1) then
            mue  = cm1(1)
          else
            mue  = 0.0d0
          endif

          test = .false.
          test = .true.
          call gnqtq (kset,x,u,ix1,ix2,xs,ch1(1,kset),ch2(1,kset),
     &                ch3(1,kset),ida,im,mr(np(192)),mue,test)
          istgn = nint(ch2(p1(4),kset))

          if (test) then
            if(istgn.eq.0) then
            elseif(istgn.eq.1) then

c             Get involved dof

              ixl(1) = ns
              masts  = nint(ch2(p1(1),kset))
              ixl(2) = ix2(1,masts)
              ixl(3) = ix2(2,masts)
              ixl(4) = ix2(3,masts)
              ixl(5) = ix2(4,masts)

              do nod2 = 1,4
                do j = 1,3
                  xm(j,nod2) = x(j,ixl(nod2+1)) + u(j,ixl(nod2+1))
                end do ! j
              end do ! nod2

c             Form stiffness and residual

              call stfqtq1 (cp0,cm1,ch2(1,kset),ch1(1,kset),
     &                      ch3(1,kset),tanm,resv,xm)
              if(iffric.eq.1) then
                mreibg2 =mreibg2 + ch3(p3( 5),kset)
                mreibg  =mreibg  + ch3(p3(15),kset)*ch3(p3(16),kset)
              endif

            elseif((istgn.eq.2).or.(istgn.eq.3)) then

c             Get involved dof

              ixl(1) = ns
              ixl(2) = nint(ch2(p1(1),kset))
              ixl(3) = nint(ch2(p1(25),kset))

              xs(1)  = x(1,ixl(1)) + u(1,ixl(1))
              xs(2)  = x(2,ixl(1)) + u(2,ixl(1))
              xs(3)  = x(3,ixl(1)) + u(3,ixl(1))
              do nod2 = 1,2
                do j = 1,3
                  xm(j,nod2) = x(j,ixl(nod2+1)) + u(j,ixl(nod2+1))
                end do ! j
              end do ! nod2

c             Form stiffness and residual

              call stfqtq2 (cp0,cm1,ch2(1,kset),ch1(1,kset),
     &                      ch3(1,kset),tanm,resv,xm)
              if(iffric.eq.1) then
                mreibg = mreibg + ch3(p3(15),kset)*ch3(p3(16),kset)
              endif

            elseif(istgn.eq.4) then

c             Get involved dof

              ixl(1) = ns
              ixl(2) = nint(ch2(p1(1),kset))

              xs(1) = x(1,ixl(1)) + u(1,ixl(1))
              xs(2) = x(2,ixl(1)) + u(2,ixl(1))
              xs(3) = x(3,ixl(1)) + u(3,ixl(1))
              do j=1,3
                xm(j,1)=x(j,ixl(2))+u(j,ixl(2))
              end do ! j

c             Form stiffness and residual

              call stfqtq4 (cp0,cm1,ch2(1,kset),ch1(1,kset),
     &                      ch3(1,kset),tanm,resv)
              if(iffric.eq.1) then
                mreibg = mreibg + ch3(p3(15),kset)*ch3(p3(16),kset)
              endif

            endif
          endif
        end do

        if(iffric.eq.1) then
          j = mr(fp(1))
          write(  *,321) 'ftrial/ix:',(ch3(p3(5),i),mr(fp(1)+i),i=1,j)
          write( 99,321) 'ftrial/ix:',(ch3(p3(5),i),mr(fp(1)+i),i=1,j)
          write(iow,321) 'ftrial/ix:',(ch3(p3(5),i),mr(fp(1)+i),i=1,j)

          write(  *,333) 'Mantelreibung:',mreibg*4,'Setzung:',u(3,1)*100
          write( 99,333) 'Mantelreibung:',mreibg*4,'Setzung:',u(3,1)*100
          write(iow,333) 'Mantelreibung:',mreibg*4,'Setzung:',u(3,1)*100
        endif

 321    format(a10,20(f10.2,i3))
 333    format(a20,f10.2,a20,f20.8)

c-----[--.----+----.----+----.-----------------------------------------]
cc--> csw <--> 103 <--> csw <--> 103 <--> csw <--> 103 <--> csw <-->
c     Contact geometry

c     calculate gap etc.

c-----[--.----+----.----+----.-----------------------------------------]

      elseif (csw.eq.103) then

c       Loop over all surface 1 elements = Loop over all slave nodes

        fp(1)  = np(191)+mr(np(192)+nsurf1-1)
        do ke = 1,mr(fp(1))

          kset  = ke
          ns    = mr(fp(1)+kset)
          xs(1) = x(1,ns) + u(1,ns)
          xs(2) = x(2,ns) + u(2,ns)
          xs(3) = x(3,ns) + u(3,ns)

c         Compute and store all geometrical parameter

          if(iffric.eq.1) then
            mue  = cm1(1)
          else
            mue  = 0.0d0
          endif
          test = .true.
          call gnqtq (kset,x,u,ix1,ix2,xs,ch1(1,kset),ch2(1,kset),
     &                ch3(1,kset),ida,im,mr(np(192)),mue,test)

        end do

c-----[--.----+----.----+----.-----------------------------------------]
cc--> csw <--> 200 <--> csw <--> 200 <--> csw <--> 200 <--> csw 200
c     Called from PMACR5 to show element informations
c-----[--.----+----.----+----.-----------------------------------------]

      elseif (csw.eq.200) then
         write (*,2000)

c-----[--.----+----.----+----.-----------------------------------------]
cc--> csw <--> 204 <--> csw <--> 204 <--> csw <--> 204 <--> csw 204
c     Printout of contact status
c-----[--.----+----.----+----.-----------------------------------------]

      elseif(csw.eq.204) then

c       Get printout flag and range

        call setcprt (ifprt,fel,lel)

c       Print title
        if (ifprt) then
          write (iow,2001) npair

        endif

c-----[--.----+----.----+----.-----------------------------------------]
cc--> csw <--> 103 <--> csw <--> 103 <--> csw <--> 103 <--> csw 103
cc--> csw <--> 304 <--> csw <--> 304 <--> csw <--> 304 <--> csw 304
c-----[--.----+----.----+----.-----------------------------------------]
c     Called from PMACR1 to reset profile                   csw = 103
c     Called from PMACR3 to reset profile                   csw = 304

      elseif ((csw.eq.103 .and.      ifistgn) .or.
     &        (csw.eq.304 .and. .not.ifistgn) ) then

c       Loop over all surface 1 elements = Loop over all slave nodes

        fp(1) = np(191)+mr(np(192)+nsurf1-1)
        fp(2) = np(191)+mr(np(192)+nsurf2-1)
        maxseg = 0
        do ke = 1,mr(fp(2))
          maxseg = max(maxseg,mr(np(191)+mr(fp(2)+mr(fp(2))+ke)))
        end do ! ke
        setval = palloc( 193,'INSEG', maxseg   ,1)
        setval = palloc( 194,'CNSEG', maxseg   ,1)
        setval = palloc( 195,'PNSEG', maxseg   ,1)
        setval = palloc( 196,'XISEG', maxseg*2 ,2)

        do i = 0,4
          istgi(i) = 0
        end do ! i

        do ke = 1,mr(fp(1))

          kset   = ke
          ns     = mr(fp(1)+kset)
          xs(1) = x(1,ns) + u(1,ns)
          xs(2) = x(2,ns) + u(2,ns)
          xs(3) = x(3,ns) + u(3,ns)

c         Compute and store all geometrical parameter

          if(iffric.eq.1) then
            mue  = cm1(1)
          else
            mue  = 0.0d0
          endif
          call geoqtq (kset,x,u,ix1,ix2,ns,xs,ch1(1,kset),ch2(1,kset),
     &                 ch3(1,kset),   ida ,mr(np(192)),mr(np(193)),
     &                 mr(np(194)),mr(np(195)),hr(np(196)),mue,csw)

          istgn        = nint(ch2(p1(4),kset))
          istgi(istgn) = istgi(istgn) + 1
          if(newflag(npair)) then

            if(ifdb .and. indb.ge.2) write(*,*)'frict. hist.var. init'

c           History initialization for normal variables

            ch1(p1( 1)  ,kset) = ch2(p1( 1)  ,kset)     !masts
            ch1(p1( 4)  ,kset) = ch2(p1( 4)  ,kset)     !istgn
            ch1(p1(19)  ,kset) = ch2(p1(19)  ,kset)     !a1
            ch1(p1(19)+1,kset) = ch2(p1(19)+1,kset)
            ch1(p1(19)+2,kset) = ch2(p1(19)+2,kset)
            ch1(p1(24)  ,kset) = ch2(p1(24)  ,kset)     !xi
            ch1(p1(24)+1,kset) = ch2(p1(24)+1,kset)     !
            ch1(p1(25)  ,kset) = ch2(p1(25)  ,kset)     !nachbar

c           History initialization for friction variables

            if (iffric.eq.1) then
              ch1(p1( 3)  ,kset) = 0.d0                 !istgt
              ch1(p1(10)  ,kset) = 0.d0                 !diffxi1
              ch1(p1(11)  ,kset) = 0.d0                 !diffxi2
              ch1(p1(26)  ,kset) = 0.d0                 !gt(1)
              ch1(p1(26)+1,kset) = 0.d0                 !gt(2)
              ch1(p1(26)+2,kset) = 0.d0                 !gt(3)
            endif

          endif
        end do
        setval = palloc( 196,'XISEG', 0,2)
        setval = palloc( 195,'PNSEG', 0,1)
        setval = palloc( 194,'CNSEG', 0,1)
        setval = palloc( 193,'INSEG', 0,1)

        newflag(npair)=.false.

        if(prt) then
          write(iow,2002) npair,mr(fp(1)),istgi
          if(ior.lt.0) then
            write(*,2002) npair,mr(fp(1)),istgi
          endif
        endif

c-----[--.----+----.----+----.-----------------------------------------]
cc--> csw <--> 305 <--> csw <--> 305 <--> csw <--> 305 <--> csw 305
c     Called from PPLOTF for plot of contact geometry
c-----[--.----+----.----+----.-----------------------------------------]

      elseif (csw.eq.305) then

        call c3geoplt (ix1,ix2,2,8)

c-----[--.----+----.----+----.-----------------------------------------]
cc--> csw <--> 308 <--> csw <--> 308 <--> csw <--> 308 <--> csw 308
cc--> csw <--> 408 <--> csw <--> 408 <--> csw <--> 408 <--> csw 408
c     Called from PPLOTF to set profile and range to plot variable
c     Called from CONTACT for plot contours of a contact variable
c-----[--.----+----.----+----.-----------------------------------------]

      elseif ((csw.eq.308) .or.
     &        (csw.eq.408)    ) then

        flag=.true.

        call c3varplt (ix1,ch1,ch2,ch3,npair,csw)

c       call c3varpltnew (ix1,ch1,ch2,ch3,npair,csw,flag)

c-----[--.----+----.----+----.-----------------------------------------]
cc--> csw <--> 310 <--> csw <--> 310 <--> csw <--> 310 <--> csw 310
c     Called from PMACR3 to reset penalty flag
c-----[--.----+----.----+----.-----------------------------------------]

      elseif (csw.eq.310) then

        if(iffron) then
          change(npair) = .true.
          chngto(npair) = .true.
        else
          change(npair) = .true.
          chngto(npair) = .false.
        endif

c-----[--.----+----.----+----.-----------------------------------------]
cc--> csw <--> 313 <--> csw <--> 313 <--> csw <--> 313 <--> csw 313
c     Called from PCONTR to initialize history variables
c-----[--.----+----.----+----.-----------------------------------------]

      elseif (csw.eq.313) then

c       Activate needed history variables

        errck = active ('masts',1)
        errck = active ('flaeche',1)
        errck = active ('nachbar',1)
        errck = active ('knflg',1)
        errck = active ('istgn',1)
        errck = active ('gn',1)
        errck = active ('nvec',3)
        errck = active ('a1',3)
        errck = active ('a2',3)
        errck = active ('fn',1)
        errck = active ('xi',2)
        if (iffric.eq.1) then
          errck = active ('istgt',1)
          errck = active ('diffxi1',1)
          errck = active ('diffxi2',1)
          errck = active ('gt',3)
          errck = active ('ft',1)
          errck = active ('ft3',1)
        endif

        if(ifaugm.ge.2) then
          errck = active ('faug',1)
        endif

        if(ifsolm.eq.2) then
c         errck = active ('lagmn',2)
          errck = active ('lagmu',1)
        endif

c       Stop variable activation defining # of data set

        fp(1) = np(191)+mr(np(192)+nsurf1-1)
        errck = active('stop',mr(fp(1)))

c-----[--.----+----.----+----.-----------------------------------------]
cc--> csw <--> 314 <--> csw <--> 314 <--> csw <--> 314 <--> csw 314
c     Called from UPDATE to update Lagrange multiplier values
c-----[--.----+----.----+----.-----------------------------------------]

      elseif (csw.eq.314) then

c       Loop over all surface 1 elements = Loop over all slave nodes

        if( ifsolm.eq.2 ) then
          fp(1)  = np(191)+mr(np(192)+nsurf1-1)
          do ke = 1,mr(fp(1))
            kset   = ke
            if(nint(ch2(p1(4),kset)).gt.0) then
              ilm(1) = mr(fp(1)+kset)
              call getlagm(ilm,1,1,ch2(p1(31),kset))
            else
              ch2(p1(31),kset) = 0.0d0
            endif
          end do ! ke
        endif ! ifsolm = 2

c-----[--.----+----.----+----.-----------------------------------------]
cc--> csw 403 <--> csw <--> 403 <--> csw <--> 403 <--> csw <--> 403    |
c     Reset profile for contacts                                       |
c-----[--.----+----.----+----.-----------------------------------------]

      elseif (csw.eq.403) then

c       Loop over all surface 1 elements = Loop over all slave nodes

        fp(1)  = np(191)+mr(np(192)+nsurf1-1)
        do ke = 1,mr(fp(1))
          kset = ke

c         Reset profile for active contacts

          istgn = nint(ch2(p1(4),kset))
          if(nint(ch3(p3(2),kset)).eq.0) then

            if (istgn.eq.1) then

              ixl(1) = mr(fp(1)+kset)
              masts  = nint(ch2(p1(1),kset))
              ixl(2) = ix2(1,masts)
              ixl(3) = ix2(2,masts)
              ixl(4) = ix2(3,masts)
              ixl(5) = ix2(4,masts)

c             Penalty -> modify profile
c             Modify profile for active contacts

              if(ifsolm.eq.1) then
                call modprof (ixl,ida,5,3)
              else
                ilm(1) = ixl(1)
                call modprofl(ixl,ida,5,3,ilm,1,1)
              endif

            elseif((istgn.eq.2).or.(istgn.eq.3))then

              ixl(1) = mr(fp(1)+kset)
              ixl(2) = nint(ch2(p1 (1) ,kset ))
              ixl(3) = nint(ch2(p1(25) ,kset ))

              if(ifsolm.eq.1) then
                call modprof (ixl,ida,3,3)
              else
                ilm(1) = ixl(1)
                call modprofl(ixl,ida,3,3,ilm,1,1)
              endif

            elseif(istgn.eq.4)then

              ixl(1) = mr(fp(1)+kset)
              ixl(2) = nint(ch2(p1 (1) ,kset ))

              if(ifsolm.eq.1) then
                call modprof (ixl,ida,2,3)
              else
                ilm(1) = ixl(1)
                call modprofl(ixl,ida,2,3,ilm,1,1)
              endif

            endif

          endif
        end do

c-----[--.----+----.----+----.-----------------------------------------]

      endif

2000  format(10x,'3-D Point to Quadrilateral Facet Contact Element',
     &       ' (QTQ)')

2001  format (/
     &'     C o n t a c t   O u t p u t   f o r   P a i r ',i5)

2002  format(/'   Pair Number =',i4,',  Total Slave Points =',i6//
     &        '    0 - Not in Contact           =',i6/
     &        '    1 - Number to a Facet        =',i6/
     &        '    2 - Number to a Concave Edge =',i6/
     &        '    3 - Number to a Convex  Edge =',i6/
     &        '    4 - Number to a Point        =',i6/ )

3000  format(3x,'*ERROR* CNTS3D must have master surface = QUAD.'/
     &          '        Input was:',i3,'-node segment surfaces.')

      end
