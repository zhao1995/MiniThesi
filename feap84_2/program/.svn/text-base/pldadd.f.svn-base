c$Id$
      subroutine pldadd(ldtab,ldnod,ldval, f,id, ns,nf, prop)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/01/2009
c       1. Add ns,nf to argument list                       05/03/2009
c       2. Increase ldtab to store spin number of displ.    09/03/2009
c       3. Add id to argument and use to set f.             07/04/2009
c          Add id on call to protnd.
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Add loads and nodes to tables

c      Inputs:
c        ldtab(4,2,*)   - Pointer, length prop no. table
c        ldnod(*)       - List of nodes
c        ldval(ndf,*)   - List of values
c        id(ndf,*)      - Boundary code values
c        ns             - Start for loop (1 or 2) -- 1 = Force
c        nf             - End   for loop (1 or 2) -- 2 = Displ
c        prop           - Total proportional load value

c      Outputs:
c        f(ndf,*)       - Force/displacement values for solution step
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'pload1.h'
      include   'prld1.h'
      include   'sdata.h'

      include   'pointer.h'
      include   'comblk.h'

      integer    ns,nf, i,j,n,nn, nod, nd,nl,npl, nprop
      real*8     prop, prp,prv, edge
      integer    ldtab(4,2,*),ldnod(*), id(ndf,*)
      real*8     ldval(ndf,*), f(ndf,*)
      real*8     theta,nv(3),xc(3),v0(3)

c     Loop over groups

      do n = 1,ldtot

c       Do forces (1) and displacements (2)

        do i = ns,nf
          nl = ldtab(2,i,n)
          if(nl.gt.0) then

c           Set proportional load value for group

            npl = ldtab(3,i,n)
            if(npl.le.0) then
              prp = prop
            else
              prp = prldv(npl)
            endif

c           Check for spin condition

            if(i.eq.2 .and. ldtab(4,2,n).gt.0) then

c             Set conditions for spin nodes

              spnum = ldtab(4,2,n)
              nd    = ldtab(1,i,n) + 1
              call pspinset(hr(np(268)),spnum,
     &                      nprop,theta,nv,xc,v0,edge)
              if(nprop.le.0) then
                prv = prop
              else
                prv = prldv(nprop)
              endif

              call protnd(nl,prp,prv,theta,nv,xc,v0,
     &                      hr(np(43)),f,id,ldnod(nd))

c           Normal Force/Displacement set

            else

c             Loop over nodes in group

              nd = ldtab(1,i,n)
              do nn = 1,nl
                nod = ldnod(nd+nn)
                do j = 1,ndf
                  if(    i.eq.1 .and. id(j,nod).eq.0) then ! Force value
                    f(j,nod) = f(j,nod) + ldval(j,nd+nn)*prp
                  elseif(i.eq.2 .and. id(j,nod).ne.0) then ! Displ value
                    f(j,nod) = f(j,nod) + ldval(j,nd+nn)*prp
                  endif
                end do ! j
              end do ! nn

            endif

          endif
        end do ! i
      end do ! n

      end
