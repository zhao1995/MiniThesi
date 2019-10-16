c$Id:$
      subroutine plocal(ld,eq,id,ix,ie,iedof,xl,ul,ule,tl,ub, x,f,u,ud,
     &                  t,un,dun, nrot, rfl,rel, jsw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add global equation numbers to ld array          27/03/2008
c       2. Add 'eltyp' to define each element topology type 03/12/2008
c       3. Separate id and eq in argument of plocal         29/04/2009
c       4. Set ul(*,*,6) for velocity at t_n                21/05/2011
c       5. Set local nurbs weight for all 'nurbfl' true     13/07/2011
c       6. Set 'hpltfl' to false                            05/01/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set local arrays for each element

c      Inputs:
c        eq(*)    - Global equation numbers
c        id(*)    - Boundary restraints
c        ie(*)    - Element descriptor parameters
c        iedof(*) - Element descriptor parameters
c        x(*)     - Global nodal coordinates
c        f(*)     - Global nodal forces/displacements
c        u(*)     - Global nodal solution parameters
c        ud(*)    - Global nodal rate parameters
c        t(*)     - Global temp variables
c        rfl      - Rigid body indicator
c        jsw      - Switching parameter

c      Scratch
c        ubl(*)   - Local array for boundary displacements

c      Outputs:
c        ld(*,*)  - Element global equation numbers
c        xl(*)    - Element nodal coordinates
c        ul(*)    - Element nodal solution parameters
c        ule(*)   - Element internal solution parameters
c        tl(*)    - Element temp values
c        ub(*)    - Element boundary displacement modify values
c        un,dun   - Boundary modification indicators
c        nrot(3)  - Number dof's with rotated directions
c        rel      - Rigid element indicator
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'c_tanfl.h'
      include   'cdata.h'
      include   'cdat1.h'
      include   'complx.h'
      include   'corset.h'
      include   'crotas.h'
      include   'ddata.h'
      include   'eldata.h'
      include   'eldatp.h'
      include   'mdata.h'
      include   'pglob1.h'
      include   'qudshp.h'
      include   'p_point.h'
      include   'rdata.h'
      include   'rdat0.h'
      include   'sdata.h'
      include   'setups.h'
      include   'part0.h'
      include   'part7.h'
      include   'pointer.h'
      include   'comblk.h'

      include   'p_int.h'

      logical    rfl,rel, rfln
      integer    nrot(3), jsw, i,j,k, iid,ild
      integer    ld(nst,*),eq(ndf,*),id(ndf,*),ix(*),ie(*),iedof(ndf,*)
      real*8     un(*),dun(*), eul(3), ang
      real*8     xl(ndm,*),ul(ndf,nen,*),ule(*),tl(*), ub(*), ubl(20)
      real*8     x(ndm,*),f(ndf,*),u(ndf,*),ud(*),t(*)

      save

c     Zero array used to store local displ, veloc, and accel

      do i = 1,nst
        ld(i,1) = 0
        ub(i)   = 0.0d0
      end do ! i

      do k = 1,7*ipc
        do j = 1,nen
          do i = 1,ndf
            ul(i,j,k) = 0.0d0
          end do ! i
        end do ! j
      end do ! k

      do i = 1,ndl
        ule(i) = 0.0d0
      end do ! i

c     Zero array used to store local tl and coordinates

      do i = 1,nen
        tl(i) = 0.0d0
        do j = 1,ndm
          xl(j,i) = 0.0d0
        end do ! j
      end do ! i

      do j = 1,ndf
        un(j)  =  0.0d0
        dun(j) =  0.0d0
      end do ! j

      rel    = .false.
      hpltfl = .false.
      plhmax =  0

c     Set up local nodal rotation array for inclined b.c.

      if(anglefl) then
        call pangl(ix,nen,hr(np(46)),hr(np(45)),nrot(1))
      else
        nrot(1) = 0
      endif
      if(eulerfl) then
        call peule(ix,nen,hr(np(243)),hr(np(242)),nrot(2))
      else
        nrot(2) = 0
      endif
      if(triadfl) then
        call pltriad(ix,nen,hr(np(275)),hr(np(274)),nrot(3))
      else
        nrot(3) = 0
      endif

c     Set element type

      eltyp = ix(nen+7)  ! N.B. FE elements have negative type
      elty2 = ix(nen+8)  ! Used for NURBS 2-d & 3-d
      elty3 = ix(nen+9)  ! Used for NURBS 3-d

c     Set individual nodal values

      do i = 1,nen

        if(ix(i).gt.0) then

c         Check for rigid body interface nodes

          rfln = .true.
          if(rfl) then
            if(mr(np(100)+ix(i)-1).ne.0) then
              rel  = .true.
              rfln = .false.
            endif
          endif

c         Set up localized solution parameters

          iid = ix(i)*ndf - ndf
          ild =     i*ndf - ndf
          nel = i
          tl(i) = t(ix(i))
          do j = 1,ndm
            xl(j,i) = x(j,ix(i))
          end do ! j
          if(eltyp.gt.0 .or. nurbfl) then
            hr(np(264)+i-1) = hr(np(263)+ix(i)-1)  ! NURB weight
          endif
          do j = 1,ndf
            ubl(j) = u(j,ix(i))
          end do ! j
          if(anglefl) then
            ang = hr(np(46)+i-1)
            if(ang.ne.0.0d0) then
              call upang(dal,ang,ubl,ndf,1)
              if(ral(1).gt.0) then
                call upang(ral,ang,ubl,ndf,1)
              endif
            endif
          endif
          if(eulerfl) then
            do j = 1,3
              eul(j) = hr(np(243)+3*i+j-4)
            end do ! j
            if(eul(1).ne.0.0d0 .or.
     &         eul(2).ne.0.0d0 .or.
     &         eul(3).ne.0.0d0) then
              call upeul(dal,eul,ubl,ndf,1)
              if(ral(1).gt.0) then
                call upeul(ral,eul,ubl,ndf,1)
              endif
            endif
          endif
          if(triadfl) then
            point = np(275) + 9*i - 9
            if(hr(point).ne.-100.0d0) then
              call uptriad(dal,hr(point),ubl,ndf,1)
              if(ral(1).gt.0) then
                call uptriad(ral,hr(point),ubl,ndf,1)
              endif
            endif
          endif
          do j = 1,ndf
            if(iedof(j,i).gt.0) then

c             Set solution, total increment, last increment

              ul(j,i,1) = u(iedof(j,i),ix(i))
              ul(j,i,2) = u(iedof(j,i),ix(i)+numnp)
              ul(j,i,3) = u(iedof(j,i),ix(i)+numnp*2)

c             Set dynamics solutions

              if(flp(9,ndfp(iedof(j,i)))) then
                k = iid+iedof(j,i)
                if(nrk.gt.0) then
                  ul(j,i,1) = ud(nrkn+k)
                endif

                if(jsw.eq.13) then
                  ul(j,i,1) = u(iedof(j,i),ix(i))
                  ul(j,i,4) = ud(k)
                else
                  if(nrc.gt.0) ul(j,i,4) = ud(nrcn+k)
                  if(nrm.gt.0) ul(j,i,5) = ud(nrmn+k)
                endif

c               Set velocity at t_n

                ul(j,i,6) = ud(nrvn+k)

c             Set acceleration for specified shift

              elseif(shflg) then
                ul(j,i,5) = -shift*ul(j,i,1)
              endif

              un(j) = max(un(j),abs(u(iedof(j,i),ix(i))))

c             Check solution method

              if(nmeth.eq.2) then
                ul(j,i,1) = ul(j,i,1) - ul(j,i,2)
              endif

c             Check for active partition

              if( ndfp(iedof(j,i)).eq.npart .or. jsw.eq.15 ) then

c               Set increment for specified boundary values

                if( id(iedof(j,i),ix(i)).gt.0 .and. rfln ) then
                  if(expflg) then
                    ul(j,i,1) = f(iedof(j,i),ix(i))
                  else
                    ub(j+ild) = f(iedof(j,i),ix(i)) - ubl(iedof(j,i))
                    if(nmeth.eq.2) then
                      ub(j+ild) = u(j+iid,2) - ul(j,i,2)
                    endif
                    dun(j) = max(dun(j),abs(ub(j+ild)))
                  endif
                endif

c               Add periodic increments for boundary correction

                k         = (ix(i) - 1)*ndf + iedof(j,i) - 1
                if(perflg) then
                  ub(j+ild) = ub(j+ild) + u(iedof(j,i),ix(i)+numnp*3)
                endif
                dun(j)    = max(dun(j),abs(ub(j+ild)))

c               Set local/global map for assembly step

                if(ddfl) then

c                 Set k for reactions

                  k = iid + iedof(j,i)
                else

c                 Set k for assembly

                  k = eq(iedof(j,i),ix(i))
                endif
              else

c               Reset k for active partitions only

                k = 0
              endif

c             Form assembly array

              ld(j+ild,1) = k

            endif
          end do ! j

c         Localize nodal rotation matrices

          if(frotas) then
            call rotloc(hr(np(82)),hr(np(83)),mr(np(81)),numnp,i,ix(i))
          endif

        endif
      end do ! i

c     Lagrange multipliers

      if(ie(nie-8).gt.0 .and. np(211).ne.0) then

c       Set equation numbers

        ild = ndf*nen
        iid = mr(np(211) + ix(nen+4) - 1) + ix(nen+5) - 1
        do j = 1,ie(nie-8)
          ld(j+ild,1) = iid + j
        end do ! j

c       Set element solution parameters

        fp(1) = np(213) + ndl*(n-1) - 1
        do j = 1,ie(nie-8)
          ule(j) = hr(fp(1)+j)
        end do ! j

      endif

c     Global equation set

      if(ie(nie-10).gt.0 .and. (gpart.eq.0 .or. gpart.eq.npart)) then
        ild = ndf*nen
        fp(1) = np(258) - 1
        do j = 1,ie(nie-10)
          ld(j+ild,1) = gneq + j
          ule(j)      = hr(fp(1)+j)
        end do ! j
      endif

c     Copy to column assembly

      do i = 1,nst
        ld(i,2) = ld(i,1)
      end do ! i

c     For complex problems set up remaining parts

      if(cplxfl) then

        call cpform(ul,ul(1,1,8),hr(np(35)),hr(np(36)),hr(np(46)),
     &              iedof,ix,un,u(1,3*numnp+1),ud(nrt*nneq+1),
     &              ndf,nst,nrot,nneq,jsw)
      endif

      end
