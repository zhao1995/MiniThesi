c$Id:$
      subroutine constass (ixl,ida,nnod,ndof,ilm,lnod,nlag,ntan,
     &                     tanm,resv)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Increase storage of idl to 200                   11/10/2011
c       2. Change du to du(1)                               01/05/2012
c       3. Add triad rotation check                         09/03/2013
c          Check for existence of parallel 'eq' (np(245))
c       4. Recode modification for displacement dofs        15/03/2013
c       5. Correct dimension of un(1) to un(10)             30/05/2013
c       6. Set 'k1' based on 'ndof' for 'du' index          25/06/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0
c               Robert L Taylor          October 28, 1996            1.1

c      Acronym: CONtact STiffness ASSembling

c      Purpose: Call assembling subroutine in blind way for the user

c      Inputs :
c         ixl(*)  - IX Local vector with involved nodes
c         ida(*)  - ID Active vector with active contact dof
c         nnod    - # of nodes in ixl
c         ndof    - # of dof in ida
c         ilm(*)  - Nodes for Lagrange Multipliers
c         lnod    - Number of multiplier nodes
c         nlag    - Number LAGrange multipliers/node
c         ntan    - SIZE of tangent matrix (dimension of array)
c         tanm(*) - TANgent Matrix
c         resv(*) - RESidual Vector

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_contac.h'
      include  'c_tanfl.h'
      include  'cdata.h'
      include  'compas.h'
      include  'corset.h'
      include  'counts.h'
      include  'eqsym.h'
      include  'idptr.h'
      include  'ndata.h'
      include  'part0.h'
      include  'pointer.h'
      include  'ptdat4.h'
      include  'rigid1.h'
      include  'sdata.h'
      include  'tdatb.h'
      include  'comblk.h'

      logical   rel
      integer   ixl(*),ida(*),ilm(*),nnod,ndof,lnod,nlag,ntan,nsiz
      integer   idl(1400),nrot(3),i,i1,i2,iu,ix,j,k1, idlag
      real*8    tanm(*),resv(*), un(10),dun(10),du(200),tol

      save

      data      tol / 1.0d-10 /

      call cdebug0 ('      constass ->',-1)

c     Use defined idl vector if big enough

      if (ntan.le.200) then

c       Get involved DOF

        do i = 1,ntan
          idl(i) = 0
        end do ! i
        call defdof (ixl,ida,idl,nnod,ndof,mr(id31))

c       Add Lagrange multipler equations to idl(*)

        i1 = nnod*ndof
        if(lnod.gt.0) then
          do i2 = 1,lnod
            if(ilm(i2).gt.0) then
              iu = idlag(mr(np(224)),ilm(i2)) - nlag
              do i = 1,nlag
                idl(i1+i) = iu + i
                du(i1+i)  = 0.0d0
              end do ! i
              i1 = i1 + nlag
            endif
          end do ! i2
        endif
        nsiz = i1

c       Transform stiffness/residual for rotated boundary conditions

        if(anglefl) then
          call pangl(ixl,nnod,hr(np(46)),hr(np(45)),nrot(1))
        else
          nrot(1) = 0
        endif
        if(eulerfl) then
          call peule(ixl,nnod,hr(np(243)),hr(np(242)),nrot(2))
        else
          nrot(2) = 0
        endif
        if(triadfl) then
          call pltriad(ixl,nnod,hr(np(275)),hr(np(274)),nrot(3))
        else
          nrot(3) = 0
        endif

c       Modify for rotated boundary conditions

        if(nrot(1)+nrot(2)+nrot(3).gt.0) then
          call ptlocal(du,resv,tanm,.false.,ndof,1,nnod,ntan,nrot,2)
        endif

c       Modify for specified boundary displacements

        do i2 = 1,ndof
          un(i2)  = 0.0d0
          dun(i2) = 0.0d0
        end do ! i2

        do i1 = 1,nnod
          k1 = ndof*(i1-1)
          do i2 = 1,ndof
            if(idl(k1+i2).le.0) then
              iu             = ndf*(ixl(i1) - 1) + ida(i2) - 1
              du(k1+ida(i2)) = (hr(np(30)+iu) - hr(np(40)+iu))*cc3
              un(i2)         = max(un(i2) ,abs(hr(np(30)+iu)),
     &                                     abs(hr(np(40)+iu)))
              dun(i2)        = max(dun(i2),abs(du(k1+ida(i2))))
            else
              du(k1+ida(i2)) = 0.0d0
            endif
          end do ! i2
        end do ! i1

        do i2 = 1,ndof
          if(dun(i2).gt.tol*un(i2)) then
            call modify(resv,tanm,du,nsiz,ntan)
            exit
          endif
        end do ! i2

c       Check for rigid body interface nodes

        rel = .false.
        if(rbody .and.nrbprt.eq.npart) then
          iu = 0
          ix = 0
          call pzero(hr(np(44)), ndm*nnod)
          call pzero(hr(np(41)),ndof*nnod)
          do i = 1,nnod
            if(mr(np(100)+ixl(i)-1).ne.0) then
              rel  = .true.
              do j = 1,ndm
                hr(np(44)+ix+j) = hr(np(43)+(ixl(i)-1)*ndm+j)
              end do
              do j = 1,ndof
                hr(np(41)+iu+ida(j)) = hr(np(40)+(ixl(i)-1)*ndf+ida(j))
              end do
            endif
            iu = iu + ndof
            ix = ix + ndm
          end do
        endif

c       Adjust assembly for parallel solutions

        if(np(245).ne.0) then
          call pparlo(idl,mr(np(245)),ixl,mr(np(31)+ndf*numnp),
     &                ntan,nnod)
          call pmodify(idl,resv,tanm,du,nsiz,ntan)
        endif

c       Transform and assemble rigid body part

        if(rel) then

          call rasbly(tanm,resv,hr(np(44)),hr(np(41)),idl,
     &                mr(np(20+npart)),ixl,mr(np(100)),mr(np(96)),
     &                mr(np(99)),hr(np(95)),ndm,ndof,nnod,nnod,ntan,
     &                lafl,uafl,dbfl,ddfl,hr(np(26)),
     &                hr(nal),hr(nau),hr(na))

c       Assemble for real arithmetic

        else

c         Assemble tangent and residual into global equations

          call dasble (tanm,resv,idl,mr(np(20+npart)),ntan,neqs,
     &                 uafl,dbfl,hr(np(26)),hr(nal),hr(nau),hr(na))
        endif

c       Store time history contact plot data

        if (ncplts.gt.0) then

          do i1 = 1,nnod
            do j = 1, ncplts
              if (icpl(1,j).eq.ixl(i1)) then
                do i2 = 1,ndof
                  if(icpl(2,j).eq.ida(i2)) then
                    cpl(j) = cpl(j) - resv(ndof*(i1-1) + i2)
                  endif
                end do ! i2
              end if
            end do ! j
          end do ! i1

        end if ! ncplts

c     Allocate scratch idl vector
c     WARNING still to be done

      else
        write (*,*) 'CONSTASS - idl vector too small'
        call plstop()
      endif

      end

      integer function idlag(ida,i2)

      include  'cdata.h'
      integer   ida(numnp,3), i2, ie

      ie    = ida(i2,2)
      if(ie.ne.0) then
        if(ida(ie,1).ne.0) then
          idlag = ida(ie,1) + ida(i2,3)
        else
          idlag = -1
        endif
      else
        idlag = -2
      endif

      end
