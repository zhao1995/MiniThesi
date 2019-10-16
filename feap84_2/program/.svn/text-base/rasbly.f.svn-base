c$Id:$
      subroutine rasbly(s,p,xl,ul,ld,jp,ix,ixt,eqrb,irb,rcg,ndm,ndf,
     &                  nel,nen,nst,alfl,aufl,bfl,dfl,b,al,au,ad)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Transform flexible connections to rigid bodies

c               This version includes options for:
c                  3: Cayley transform (THETA_n+a)
c                 -3: Cayley transform (THETA_n+1)
c                 -4: Exponential map
c                 -6: Explicit central difference
c                 -7: Linear

c      Inputs:
c         s(*,*)    - Element tangent matrix
c         p(*)      - Element residual vector
c         xl(ndm,*) - Element nodal coordinates
c         ul(ndf,*) - Element solution/rate vector
c         ld(*)     - Element global equation numbers
c         jp(*)     - Pointer to end of row/columns in profile
c         ix(*)     - Element nodal connections
c         ixt(*)    - Nodal rigid body numbers (0 = flexible node)
c         eqrb(*)   - Rigid body update method
c         irb(*)    - Rigid body equation numbers
c         rcg(3,*,*)- Rigid body translational solution/rates
c         ndm       - Spatial dimension of mesh
c         ndf       - Number dof/node
c         nel       - Active number of nodes/element
c         nen       - Maximum number nodes/element
c         nst       - Dimension of s-array
c         alfl      - Flag, assemble unsymmetric matrix if true.
c         aufl      - Flag, assemble tangent if true
c         bfl       - Flag, assemble compressed residual if true
c         dfl       - Flag, form uncompressed reaction if true

c      Outputs:
c         b(*)      - Residual with rigid body coupling terms
c         al(*)     - Lower tangent with rigid body coupling terms
c         au(*)     - Upper tangent with rigid body coupling terms
c         ad(*)     - Diagonal tangent with rigid body coupling terms
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'ddata.h'
      include  'eqsym.h'
      include  'iofile.h'
      include  'rigid1.h'

      logical   aufl,alfl,bfl,dfl, nonlin
      integer   i,j,k,l, i1,i2,i3,ii,il,jf,jr, j1, ndm,ndf,nel,nen,nst
      integer   efl(27),erl(27), lr(48),lrb(48)
      integer   ix(nel),ixt(*),eqrb(*),irb(nrbdof,*),ld(ndf,nel),jp(*)
      real*8    pdotr,factr, facth
      real*8    p(nst),s(nst,nst), xl(ndm,nel), ul(ndf,nen,2)
      real*8    b(*),ad(*),au(*),al(*), pr(48), sr(48,48)
      real*8    rcg(3,11,*), r(3,27), ra(3,27), rh(3,27), h1(3,6)

      save

c     Storage/address information:

c       jr : # rigid    nodes on element
c       jf : # flexible nodes on element
c       il = last flexible dof on local element.

c       pr : rigid/flexible local residual
c       sr : rigid/flexible local tangent matrix
c       lr : local/global numbers for assembly

c       eqr : global equation number -1 for rigid node 'ii'
c       erl : local position in sr for 'jr' rigid''' node.
c       efl = local position in sr for 'jf' flexible node
c       lrb : temporary storage for local-global dof map.

      jr = 0
      jf = 0
      il = 0

      call pzeroi(lr, 48   )

      i1 = 0
      i2 = 0

      do i = 1,nel

        ii = ix(i)

        if( ii.gt.0 ) then

c         Check if node rigid [ ixt(ii) > 0 ]

          if( ixt(ii).gt.0 ) then

c           Set tangent factors for transformations

            nonlin = .true.
            if(eqrb(ixt(ii)).eq.3) then
              facth = 1.d0/theta(3)
              factr = facth - 1.0d0
            elseif(eqrb(ixt(ii)).eq.-3) then
              facth = theta(3)
              factr = 1.d0/theta(3) - 1.0d0
            elseif(eqrb(ixt(ii)).eq.-6) then
              facth = 0.0d0
              factr = 0.0d0
            elseif(eqrb(ixt(ii)).eq.-7) then
              facth = 0.0d0
              factr = 0.0d0
              nonlin= .false.
            else
              facth = 1.0d0
              factr = 0.0d0
            endif

c           Rigid computations:

            jr      = jr + 1
            erl(jr) = i1
            if(dfl) then
              do j = 1,nrbdof
                lrb(i2+j) = 0
              end do ! j
              do j = 1,ndf
                lrb(i2+j) = ld(j,i)
              end do ! j
            else
              do j = 1,nrbdof
                lrb(i2+j) = irb(j,ixt(ii))
              end do ! j
            end if
            i2     = i2 + nrbdof

c           Set up relative position vector for rigid node 'jr'

            r(3,jr)  = 0.0d0
            ra(3,jr) = 0.0d0
            rh(3,jr) = 0.0d0
            if(nonlin) then
              do j = 1,ndm
                r(j,jr)  = xl(j,i) + ul(j,i,1) + factr*ul(j,i,2)
     &                                         - rcg(j,2,ixt(ii))
                ra(j,jr) = xl(j,i) + ul(j,i,1) - rcg(j,7,ixt(ii))
                rh(j,jr) = facth*r(j,jr)
              end do ! j
            else
              do j = 1,ndm
                r(j,jr)  = xl(j,i) - rcg(j,1,ixt(ii))
                ra(j,jr) = xl(j,i) - rcg(j,1,ixt(ii))
                rh(j,jr) = facth*r(j,jr)
              end do ! j
            endif

          else

c           Flexible computations:

            do j = 1,ndf
              lr(il+j) = ld(j,i)
            end do ! j
            jf      = jf + 1
            efl(jf) = i1
            il      = il + ndf

          end if
        end if

        i1 = i1 + ndf

      end do ! i

c     Adjust lr array to include rigid equations

      do i = 1,i2
        lr(il+i) = lrb(i)
      end do ! i

c     Assembly for rigid dof: ndm = 2 or 3

      if( bfl ) then

        i1 = 0
        do i = 1,jf

c         Translational dof: R_f = P_a

          do j = 1,ndf
            pr(i1+j) = p(efl(i)+j)
          end do ! j

          i1 = i1 + ndf
        end do ! i

        i1 = il
        do i = 1,jr

c         Translational dof: R_rb = P_a

          if(ndm.eq.2) then
            pr(i1+3) = 0.0d0
          elseif(ndm.eq.3) then
            pr(i1+4) = 0.0d0
            pr(i1+5) = 0.0d0
            pr(i1+6) = 0.0d0
          endif

          do j = 1,ndf
            pr(i1+j) = p(erl(i)+j)
          end do ! j

c         Rotational dof: M_rb = r_a x P_a

          if(ndm.eq.2) then
            pr(i1+3) = pr(i1+3)
     +               + ra(1,i)*p(erl(i)+2) - ra(2,i)*p(erl(i)+1)
          elseif(ndm.eq.3) then

            pr(i1+4) = pr(i1+4)
     +               + ra(2,i)*p(erl(i)+3) - ra(3,i)*p(erl(i)+2)
            pr(i1+5) = pr(i1+5)
     +               + ra(3,i)*p(erl(i)+1) - ra(1,i)*p(erl(i)+3)
            pr(i1+6) = pr(i1+6)
     +               + ra(1,i)*p(erl(i)+2) - ra(2,i)*p(erl(i)+1)
          endif

          i1 = i1 + nrbdof
        end do ! i

      end if

c     Assemble Stiffness dof: ndm = 2 or 3

      if( aufl ) then

c       Set rigid-flexible stiffness to zero

        do i = 1,48
          do j = 1,48
            sr(i,j) = 0.0d0
          end do ! j
        end do ! i

c       Linearization of flexible residual: coupling to rigid body

        i1 = 0
        do i = 1,jf

          j1 = 0
          do j = 1,jf

c           Translational flexible to flexible nodes

            do l = 1,ndf
              do k = 1,ndf
                sr(i1+k,j1+l) = s(efl(i)+k,efl(j)+l)
              end do ! k
            end do ! l

            j1 = j1 + ndf
          end do ! j

          j1 = il
          do j = 1,jr

c           Translational dof coupling to flexible nodes

            do l = 1,ndf
              do k = 1,ndf
                sr(i1+k,j1+l) = s(efl(i)+k,erl(j)+l)
              end do ! k
            end do ! l

c           Rotational dof coupling to flexible nodes

            if(ndm.eq.2) then
              do k = 1,ndf
                sr(i1+k,j1+3) = sr(i1+k,j1+3)
     +                        + s(efl(i)+k,erl(j)+2)*r(1,j)
     +                        - s(efl(i)+k,erl(j)+1)*r(2,j)
              end do ! k
            elseif(ndm.eq.3) then
              do k = 1,ndf
                sr(i1+k,j1+4) = sr(i1+k,j1+4)
     +                        + s(efl(i)+k,erl(j)+3)*r(2,j)
     +                        - s(efl(i)+k,erl(j)+2)*r(3,j)
                sr(i1+k,j1+5) = sr(i1+k,j1+5)
     +                        + s(efl(i)+k,erl(j)+1)*r(3,j)
     +                        - s(efl(i)+k,erl(j)+3)*r(1,j)
                sr(i1+k,j1+6) = sr(i1+k,j1+6)
     +                        + s(efl(i)+k,erl(j)+2)*r(1,j)
     +                        - s(efl(i)+k,erl(j)+1)*r(2,j)
              end do ! k
            endif

            j1 = j1 + nrbdof
          end do ! j
          i1 = i1 + ndf
        end do ! i

c       Linearization of rigid body residual: coupling to flexible body

        i1 = il
        do i = 1,jr

          j1 = 0
          do j = 1,jf

c           Translational dof coupling to flexible nodes
            do l = 1,ndf
              do k = 1,ndf

                sr(i1+k,j1+l) = s(erl(i)+k,efl(j)+l)

              end do ! k
            end do ! l

c           Rotational dof coupling to flexible nodes

            if(ndm.eq.2) then
              do k = 1,ndf
                sr(i1+3,j1+k) = sr(i1+3,j1+k)
     +                        + ra(1,i)*s(erl(i)+2,efl(j)+k)
     +                        - ra(2,i)*s(erl(i)+1,efl(j)+k)
              end do ! k
            elseif(ndm.eq.3) then
              do k = 1,ndf

                sr(i1+4,j1+k) = sr(i1+4,j1+k)
     +                        + ra(2,i)*s(erl(i)+3,efl(j)+k)
     +                        - ra(3,i)*s(erl(i)+2,efl(j)+k)
                sr(i1+5,j1+k) = sr(i1+5,j1+k)
     +                        + ra(3,i)*s(erl(i)+1,efl(j)+k)
     +                        - ra(1,i)*s(erl(i)+3,efl(j)+k)
                sr(i1+6,j1+k) = sr(i1+6,j1+k)
     +                        + ra(1,i)*s(erl(i)+2,efl(j)+k)
     +                        - ra(2,i)*s(erl(i)+1,efl(j)+k)
              end do ! k
            end if

            j1 = j1 + ndf
          end do ! j

c         Linearization of rigid body residual: coupling to rigid body

          i3 = i1 + 3
          j1 = il
          do j = 1,jr

c           Translational dof coupling to translational rigid nodes

            do k = 1,ndf
              do l = 1,ndf
                sr(i1+k,j1+l) = s(erl(i)+k,erl(j)+l)
              end do ! l
            end do ! k

c           Rotational dof coupling to translational rigid nodes

            if(ndm.eq.2) then
              do l = 1,ndf
                sr(i1+l,j1+3) = sr(i1+l,j1+3)
     +                        + s(erl(i)+l,erl(j)+2)*r(1,j)
     +                        - s(erl(i)+l,erl(j)+1)*r(2,j)

                h1(3,l)       = ra(1,i)*s(erl(i)+2,erl(j)+l)
     +                        - ra(2,i)*s(erl(i)+1,erl(j)+l)

                sr(i1+3,j1+l) = sr(i1+3,j1+l) + h1(3,l)

              end do ! l

              sr(i1+3,j1+3) = sr(i1+3,j1+3)
     +                      + h1(3,2)*r(1,j) - h1(3,1)*r(2,j)

            elseif(ndm.eq.3) then
              do l = 1,ndf
                sr(i1+l,j1+4) = sr(i1+l,j1+4)
     +                        + s(erl(i)+l,erl(j)+3)*r(2,j)
     +                        - s(erl(i)+l,erl(j)+2)*r(3,j)
                sr(i1+l,j1+5) = sr(i1+l,j1+5)
     +                        + s(erl(i)+l,erl(j)+1)*r(3,j)
     +                        - s(erl(i)+l,erl(j)+3)*r(1,j)
                sr(i1+l,j1+6) = sr(i1+l,j1+6)
     +                        + s(erl(i)+l,erl(j)+2)*r(1,j)
     +                        - s(erl(i)+l,erl(j)+1)*r(2,j)

                h1(1,l)       = ra(2,i)*s(erl(i)+3,erl(j)+l)
     +                        - ra(3,i)*s(erl(i)+2,erl(j)+l)
                h1(2,l)       = ra(3,i)*s(erl(i)+1,erl(j)+l)
     +                        - ra(1,i)*s(erl(i)+3,erl(j)+l)
                h1(3,l)       = ra(1,i)*s(erl(i)+2,erl(j)+l)
     +                        - ra(2,i)*s(erl(i)+1,erl(j)+l)

                sr(i1+4,j1+l) = sr(i1+4,j1+l) + h1(1,l)
                sr(i1+5,j1+l) = sr(i1+5,j1+l) + h1(2,l)
                sr(i1+6,j1+l) = sr(i1+6,j1+l) + h1(3,l)

              end do ! l

c             Rotational dof coupling to rotational rigid nodes

              do l = 1,3
                sr(i3+l,j1+4) = sr(i3+l,j1+4)
     +                        + h1(l,3)*r(2,j) - h1(l,2)*r(3,j)
                sr(i3+l,j1+5) = sr(i3+l,j1+5)
     +                        + h1(l,1)*r(3,j) - h1(l,3)*r(1,j)
                sr(i3+l,j1+6) = sr(i3+l,j1+6)
     +                        + h1(l,2)*r(1,j) - h1(l,1)*r(2,j)
              end do ! l
            end if

            j1 = j1 + nrbdof
          end do ! j

c         Add geometric stiffness term

          if(ndm.eq.2) then

            sr(i3,i3) = sr(i3,i3)
     +                + pr(i1+1)*rh(1,i) + pr(i1+2)*rh(2,i)

          elseif(ndm.eq.3) then

            pdotr = pr(i1+1)*rh(1,i)
     +            + pr(i1+2)*rh(2,i)
     +            + pr(i1+3)*rh(3,i)

            do k = 1,3

              sr(i3+k,i3+k) = sr(i3+k,i3+k) + pdotr

              if(alfl) then

                do l = 1,3
                  sr(i3+k,i3+l) = sr(i3+k,i3+l) - rh(k,i)*pr(i1+l)
                end do ! l

              else

                do l = 1,3
                  sr(i3+k,i3+l) = sr(i3+k,i3+l)
     &              - (rh(k,i)*pr(i1+l) + rh(l,i)*pr(i1+k))*0.5d0
                end do ! l

              endif

            end do ! k

          endif

          i1 = i1 + nrbdof

        end do ! i

      endif

c     Perform assembly into global arrays

      call dasble(sr,pr,lr,jp,48,neqs,aufl,bfl,b,al,au,ad)

      end
