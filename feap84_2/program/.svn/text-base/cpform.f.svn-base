c$Id:$
      subroutine cpform(ulr,uli,p,s,angl,iedof,ix,un,
     &                  ui,udi,ndf,nst,nrot,nneq,jsw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add set of nrvn for velocity at t_n              21/05/2011
c       2. Add triad rotation                               11/02/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set element variables for complex imaginary parts

c      Inputs:
c         angl(*)      - Value of angle for sloping boundary conditions
c         iedof(ndf,*) - Assembly information for material sets
c         ix(*)        - Node numbers connected to element
c         ui(ndf,*)    - Imaginary part of solutions
c         udi(*)       - Imaginary part of rate terms
c         ndf          - Number dof/node
c         nst          - Size of elemnt arrays
c         nrot         - Number nodes requiring rotation transformation
c         nneq         - Total number of terms in solution arrays
c         jsw          - Solution option

c      Outputs:
c         ulr(*)       - reali     parts of local solution variables
c         uli(*)       - Imaginary parts of local solution variables
c         s(nst,nst)   - Zeroed array
c         p(nst)       - Zeroed array
c         un           - Maximum entry in solution
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'ddata.h'
      include  'eldata.h'
      include  'eltran.h'
      include  'mdata.h'
      include  'part0.h'
      include  'part7.h'
      include  'rdata.h'
      include  'pointer.h'
      include  'comblk.h'

      integer   i, j, k, nst, jsw, numnp2, nneq,ndf, nrot(3)
      integer   iid, iedof(ndf,*), ix(*)
      real*8    ulr(ndf,nen,*), uli(ndf,nen,*), ui(ndf,*), udi(*)
      real*8    p(nst,*), s(nst,nst,*), angl(*), un(*), omega

      save

c     Set indexing parameters

      omega = sqrt(abs(shift))  ! shift = omega*omega
      numnp2 = numnp*2
      nneq   = numnp*ndf
      nrkn   = nrk*nneq - nneq
      nrcn   = nrc*nneq - nneq
      nrmn   = nrm*nneq - nneq
      nrvn   = nrt*nneq - nneq - nneq

c     Initialize residual and array

      do i = 1,nst
        do j = 1,nst
          s(j,i,1) = 0.0d0 ! real      part
          s(j,i,2) = 0.0d0 ! imaginary part
        end do ! j
        p(i,1) = 0.0d0 ! real      part
        p(i,2) = 0.0d0 ! imaginary part
      end do ! i

c     Set imaginary parts of element arrays

      do i = 1,nen

        if(ix(i).gt.0) then
          iid = ix(i)*ndf - ndf
          do j = 1,ndf
            if(iedof(j,i).gt.0) then
              uli(j,i,1) = ui(iedof(j,i),ix(i))
              uli(j,i,2) = ui(iedof(j,i),ix(i)+numnp)
              uli(j,i,3) = ui(iedof(j,i),ix(i)+numnp2)
              if(flp(9,ndfp(iedof(j,i)))) then
                k = iid + iedof(j,i)
                if(nrk.gt.0) uli(j,i,1) = udi(nrkn+k)

                if(jsw.eq.13) then
                   uli(j,i,1) = ui(iedof(j,i),ix(i))
                   uli(j,i,5) = udi(k)
                elseif (noi.eq.5) then
                   uli(j,i,4) = udi(k)
                   uli(j,i,5) = udi(nneq+k)
                   uli(j,i,6) = udi(nrmn+k)
                else
                   if(nrc.gt.0) uli(j,i,4) = udi(nrcn+k)
                   if(nrm.gt.0) uli(j,i,5) = udi(nrmn+k)
                end if

c             Set velocity and acceleration for specified shift

              elseif(shflg) then
                uli(j,i,4) =  omega*ulr(j,i,1)
                ulr(j,i,4) = -omega*uli(j,i,1)
                uli(j,i,5) = -shift*uli(j,i,1)
                ctan(2)      =  0.0d0
              endif

              un(j) = max(un(j),abs(uli(j,i,1)))

            end if

          end do ! j
        end if

      end do ! i

c     Rotate angle parameters if necessary

      if(nrot(1).gt.0) then
        call ptrans(dal,angl,uli,p,s,nel,ndf,nst,1)
        if(ral(1).gt.0) then
          call ptrans(ral,angl,uli,p,s,nel,ndf,nst,1)
        end if
      end if

c     Rotate Euler parameters if necessary

      if(nrot(2).gt.0) then
        call petrans(dal,hr(np(243)),uli,p,s,nel,ndf,nst,1)
        if(ral(1).ne.0) then
          call petrans(dal,hr(np(243)),uli,p,s,nel,ndf,nst,1)
        endif
      endif

c     Rotate triad parameters if necessary

      if(nrot(3).gt.0) then
        call pttrans(dal,hr(np(275)),uli,p,s,nel,ndf,nst,1)
        if(ral(1).ne.0) then
          call pttrans(dal,hr(np(275)),uli,p,s,nel,ndf,nst,1)
        endif
      endif

      end
