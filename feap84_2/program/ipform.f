c$Id:$
      subroutine ipform(ul1,ul2,xl1,xl2,tl1,tl2,ld,p,s,
     &                  ie,d,id,x,ix,intel,f,t,jp,
     &                  u,ud,b,a,al,aufl,bfl,alfl,dfl,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Separate id and eq on call to plocal             29/04/2009
c       2. Add set of nrvn for velocity at t_n              21/05/2011
c       3. Dimension un(20), dun(20)                        15/06/2011
c       4. Increase to nrot1(3) & nrot2(3)                  11/03/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute interface element arrays and assemble into
c               global arrays

c      Inputs:
c         ie(nie,*)   - Assembly information for material set
c         d(ndd,*)    - Material set parameters
c         id(ndf,*)   - Equation numbers for each active dof
c         x(ndm,*)    - Nodal coordinates of mesh
c         ix(nen1,*)  - Element nodal connections of mesh
c         rben(*)     - Rigid/modal element indicator
c         f(ndf,*,2)  - Nodal force and displacement values
c         t(*)        - Nodal temperature values
c         jp(*)       - Pointer array for row/columns of tangent
c         u(*)        - Nodal solution values
c         ud(*)       - Nodal rate values
c         ndd         - Dimension for d array
c         nie         - Dimension for ie array
c         ndf         - Number dof/node
c         ndm         - Spatial dimension of mesh
c         nen1        - Dimension for ix array
c         nst         - Dimension for element array
c         aufl        - Flag, assemble coefficient array if true
c         bfl         - Flag, assemble vector if true
c         alfl        - Flag, coefficient array unsymmetric if true
c         dfl         - Flag, assemble reactions if true
c         isw         - Switch to control quantity computed

c      Local element arrays:
c         ul1(ndf,*)  - Element solution and rate values
c         ul2(ndf,*)  - Element solution and rate values
c         xl1(ndm,*)  - Element nodal coordinates
c         xl2(ndm,*)  - Element nodal coordinates
c         ld(*)       - Element local/global equation numbers
c         p(*)        - Element vector
c         s(nst,*)    - Element array

c      Outputs:
c         b(*)        - Global vector
c         a(*)        - Global matrix, diagonal and upper part
c         al(*)       - Global matrix, lower part
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'arcler.h'
      include  'cdata.h'
      include  'cdat1.h'
      include  'complx.h'
      include  'ddata.h'
      include  'elcount.h'
      include  'eldata.h'
      include  'hdatam.h'
      include  'ieldat.h'
      include  'lmdata.h'
      include  'pointer.h'
      include  'prld1.h'
      include  'prlod.h'
      include  'sdata.h'
      include  'comblk.h'

      include  'p_int.h'

      logical   aufl,bfl,alfl,dfl,efl,gfl,rel,rfl,arotflg
      logical   actrg,nobnd
      integer   isw, jsw, ksw
      integer   i, n1,n2,nint, nov
      integer   nrot1(3), nrot2(3), nub
      integer   ld(*),ie(nie,*),id(ndf,*),ix(nen1,*),intel(8,*),jp(*)
      integer   intnod(5)
      real*8    ul1(*),ul2(*),xl1(*),xl2(*),tl1(*),tl2(*),p(*),s(*)
      real*8    d(ndd,*),x(ndm,*) ,f(ndf,numnp),u(ndf,*),ud(*),t(*)
      real*8    un(20), dun(20), prope, b(*), a(*), al(*)

      save

c     Set element proportional loading value

      prope = (theta(3)*(prop - propo) + propo)
      if(isw.ne.23) then
        prope = prope*rlnew
      endif

c     Recover nh1, nh2, nh3 pointers for both elements

      n1h1 = np(50)
      n1h2 = np(51)
      n1h3 = np(52)

      n2h1 = np(50) + nhmax
      n2h2 = np(51) + nhmax
      n2h3 = np(52) + nh3max

      nih1 = np(215)
      nih2 = np(216)
      nih3 = np(217)

c     Set parameters for solution

      nsts = nst*2
      nub  = nsts*2 + 1
      iel1 = 0
      iel2 = 0
      efl  = .false.
      if(.not.dfl.and.isw.eq.6) efl = .true.
      if(bfl.and.isw.eq.3)      efl = .true.
      arotflg = aufl .or. bfl

      if(isw.eq.19) then
        if(bfl) efl = .true.
        jsw = 5
        ksw = 5
        gfl = .false.
      else
        jsw = isw
        ksw = 3
        gfl = .true.
      endif

      nneq   = numnp*ndf
      nrkn   = nrk*nneq - nneq
      nrcn   = nrc*nneq - nneq
      nrmn   = nrm*nneq - nneq
      nrvn   = nrt*nneq - nneq - nneq

c     Set for no rigid body

      rfl = .false.

c     Loop over interface

      nint = 1
      do while(intel(1,nint).ne.0)

        n1        = intel(1,nint)
        n2        = intel(2,nint)
        intnod(1) = intel(3,nint)
        intnod(2) = intel(4,nint)
        intnod(3) = intel(5,nint)
        intnod(4) = n1
        intnod(5) = n2
        nint      = nint + 1
        ma1       = ix(nen1,n1)
        iel1      = ie(nie-1,ma1)

c       Set interface history variables

        if(intel(7,nint)-intel(6,nint).gt.0) then
          nit1 = np(214) + intel(6,nint)
          nit2 = np(214) + intel(7,nint)
          do i = 0,intel(7,nint)-intel(6,nint)-1
            hr(nih1+i) = hr(nit1+i)
            hr(nih2+i) = hr(nit2+i)
          end do ! i
        endif

        if(intel(8,nint)-intel(6,nint+1).gt.0) then
          nit3 = np(214) + intel(8,nint)
          do i = 0,intel(8,nint)-intel(6,nint+1)-1
            hr(nih3+i) = hr(nit3+i)
          end do ! i
        endif

c       Check for boundary and active regions

        if( n2.gt.0 ) then
          actrg = ix(nen1-1,n1).ge.0 .and. ix(nen1-1,n2).ge.0
          nobnd = .true.
          ma2   = ix(nen1,n2)
          iel2  = ie(nie-1,ma2)
        else
          actrg = ix(nen1-1,n1).ge.0
          nobnd = .false.
          n2    = n1
          ma2   = ma1
          iel2  = 0
          call pzeroi(ld(nst+1),nst)
        endif

        if(actrg) then

c         Set local and history variables for element 1

          if(ie(nie-2,ma1).eq.ix(nen1,n1)) then

            if(ie(nie,ma1).gt.0) then
              n1t1 = np(49) + ix(nen+1,n1) + ie(nie-3,ma1)
              n1t2 = np(49) + ix(nen+2,n1) + ie(nie-3,ma1)
              do i = 0,ie(nie,ma1)-1
                hr(n1h1+i) = hr(n1t1+i)
                hr(n1h2+i) = hr(n1t2+i)
              end do ! i
            endif

            if(ie(nie-5,ma1).gt.0) then
              n1t3 = np(49) + ix(nen+3,n1) + ie(nie-4,ma1)
              do i = 0,ie(nie-5,ma1)-1
                hr(n1h3+i) = hr(n1t3+i)
              end do ! i
            endif
          endif

c         Set local arrays and rotate to global frame

          n     = n1
          fp(1) = ndf*nen*(ma1-1) + np(240)              ! iedof
          call plocal(ld,id,mr(np(31)+nneq),ix(1,n1),ie(1,ma1),
     &                mr(fp(1)),xl1,ul1,ule1,tl1,p(nub),x,f,u,ud,t,
     &                un,dun,nrot1, rfl,rel, jsw)
          if(nrot1(1)+nrot1(2)+nrot1(3).gt.0) then
            call ptlocal(ul1,p,s, .false., ndf,nen,nel,nst,nrot1, 1)
          endif

c         Set local and history variables for element 2

          if(nobnd) then

            if(ie(nie-2,ma2).eq.ix(nen1,n2)) then

              if(ie(nie,ma2).gt.0) then
                n2t1 = np(49) + ix(nen+1,n2) + ie(nie-3,ma2)
                n2t2 = np(49) + ix(nen+2,n2) + ie(nie-3,ma2)
                do i = 0,ie(nie,ma1)-1
                  hr(n1h1+i) = hr(n1t1+i)
                  hr(n1h2+i) = hr(n1t2+i)
                end do ! i
              endif

              if(ie(nie-5,ma1).gt.0) then
                n2t3 = np(49) + ix(nen+3,n2) + ie(nie-4,ma2)
                do i = 0,ie(nie-5,ma1)-1
                  hr(n1h3+i) = hr(n1t3+i)
                end do ! i
              endif
            endif

c           Set local arrays for element 2

            n     = n2
            fp(1) = ndf*nen*(ma1-1) + np(240)              ! iedof
            call plocal(ld(nst+1),id,mr(np(31)+nneq),ix(1,n2),ie(1,ma2),
     &                  mr(fp(1)),xl2,ul2,ule2,tl2,p(nst+nub),x,f,u,ud,
     &                  t,un,dun,nrot2,rfl,rel, jsw)
            if(nrot2(1)+nrot2(2)+nrot2(3).gt.0) then
              call ptlocal(ul2,p,s, .false., ndf,nen,nel,nst,nrot2, 1)
            endif

c         Boundary arrays:

          else

            call pzero(xl2 , ndm*nen)
            call pzero(ul2 , ndf*nen*7)
            call pzero(tl2 ,     nen)
            call pzero(ule2, 100)

          endif

c         Call element library

          dm = prope
          call ielmlib(d(1,ma1),ul1,xl1,ix(1,n1),tl1,
     &                 d(1,ma2),ul2,xl2,ix(1,n2),tl2,
     &                 intnod,s,p,jsw)

c         Modify for rotated dof's : PROBABLY DOES NOT WORK

          if(nrot1(1)+nrot1(2)+nrot1(3).gt.0 .and. arotflg) then
            call ptlocal(ul1,p,s, cplxfl, ndf,nen,nel,nst,nrot1, 2)
          endif

c         Assemble arrays as necessary

          if(jsw.eq.3 .or. jsw.eq.5 .or. jsw.eq.6 .or.jsw.eq.9) then

            call passble(s,p,ld,ix(1,n1), jp,a,al,b,
     &                   alfl,aufl,bfl,dfl,gfl,rel, nsts,nov, jsw)
          endif

c         Position update terms 'nt1,nt2' from 'nh1,nh2' to save

          if(hflgu) then

c           Interface history

            if(intel(7,nint)-intel(6,nint).gt.0) then
              do i = 0,intel(7,nint)-intel(6,nint)-1
                hr(nit1+i) = hr(nih1+i)
                hr(nit2+i) = hr(nih2+i)
              end do ! i
            endif

            if(intel(8,nint)-intel(6,nint+1).gt.0) then
              do i = 0,intel(8,nint)-intel(6,nint+1)-1
                hr(nit3+i) = hr(nih3+i)
              end do ! i
            endif

c           Element 1 history

            do i = 0,ie(nie,ma1)-1
              hr(n1t1+i) = hr(n1h1+i)
              hr(n1t2+i) = hr(n1h2+i)
            end do ! i

            do i = 0,ie(nie-5,ma1)-1
              hr(n1t3+i) = hr(n1h3+i)
            end do ! i

c           Element 2 history

            if(nobnd) then
              do i = 0,ie(nie,ma2)-1
                hr(n2t1+i) = hr(n2h1+i)
                hr(n2t2+i) = hr(n2h2+i)
              end do ! i

              do i = 0,ie(nie-5,ma2)-1
                hr(n2t3+i) = hr(n2h3+i)
              end do ! i
            endif
          endif

        end if ! regions

      end do ! while int(1,nint) > 0

      end
