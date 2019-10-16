c$Id:$
      subroutine pndiff(ix,ie,ul,st,k1,du,flg)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Change 'j' to 'fp(1)' at line 271                17/04/2011
c       2. Set pointer to 'fp(3)' for property array        21/06/2011
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute tangent by numerical differentiation using
c               central difference

c      Inputs:
c        ix(*)          - Node connection list & material number
c        ie(nie,*)      - Element information
c        ul(ndf,nen,*)  - Array for local solutions
c        st(nst,nst,2)  - Stiffness storage
c        k1             - Element number to compute
c        du             - Increment to use for computation
c        flg            - Compute analytical tangent if true

c      Outputs:
c        st(nst,nst,2)  - Stiffness: 1 = numerical tangent
c                                    2 = analytical tangent
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'cdat1.h'
      include  'ddata.h'
      include  'debugs.h'
      include  'eldata.h'
      include  'eltran.h'
      include  'idptr.h'
      include  'iofile.h'
      include  'lmdata.h'
      include  'part0.h'
      include  'pointer.h'
      include  'sdata.h'
      include  'tdata.h'
      include  'tdatb.h'
      include  'comblk.h'

      include  'p_int.h'

      logical   flg, elflg
      integer   i,j,k,m, mel, k1, ix(nen1), ie(nie,*)
      integer   nrw, ncl
      real*8    du,sqrtdig,usize(30),usizel, denom(30), denomel
      real*8    temp1,temp2,temp3,temp4,temp5
      real*8    ul(ndf,nen,*),st(nst,nst,2)

      save

c     Tolerance for numerical diffentiation

      data      sqrtdig /1.d-8/

c     Numerical differentiation to compute element tangent

      call pzero(st,nst*nst)

c     Compute current residual

      if(flg) then
        call pform(ul,hr(np(44)),hr(np(39)),mr(np(34)),
     &             hr(np(35)),hr(np(36)),mr(np(32)),hr(np(25)),
     &             mr(id31),hr(np(43)),mr(np(33)),mr(np(181)),
     &             hr(np(30)),hr(np(38)),mr(np(20+npart)),hr(np(40)),
     &             hr(np(42)),hr(np(26)),hr(np(26)),hr(np(26)),
     &             ndd,nie,ndf,ndm,nen1,nst,
     &             .false.,.false.,.false.,.true.,3,k1,k1,1)

c       Save element tangent

        do m = 1,nst
          fp(1) = np(36) + (m-1)*nst - 1
          do i = 1,nst
            st(i,m,2) = hr(fp(1)+i)
          end do ! i
        end do ! m
      endif

      do j = 1,ndf
        usize(j) = 0.0d0
        do i = 1,nen
          usize(j) = max(usize(j),abs(ul(j,i,1)))
        end do ! i
      end do ! j
      if(debug) then
        write(iow,*) ' Computed max U:',cc1,cc2,c3,c4,c5,du
        call mprint(usize,1,ndf,1,'USIZE')
      endif

      do j = 1,ndf
        if(du.ne.0.0d0) then
          usize(j) = du
        elseif(usize(j).ne.0.0d0) then
          usize(j) = sqrtdig*usize(j)
        else
          usize(j) = sqrtdig
        endif
        denom(j) = 0.5d0/usize(j)
      end do ! j

c     Check for element variables

      elflg = .false.
      ma    = ix(nen1)
      fp(3) = np(25) + ndd*(ma - 1)
      if(ie(nie-8,ma).gt.0) then
        fp(2)  = np(213) + ndl*(k1-1) - 1
        elflg  = .true.
        usizel = 0.0d0
        do j = 1,ie(nie-8,ma)
          ule(j) = hr(fp(2)+j)
          usizel = max(usizel,abs(ule(j)))
        end do ! j
        if(usizel.ne.0.0d0) then
          usizel = sqrtdig*usizel
        else
          usizel = sqrtdig
        endif
        if(debug) then
          call mprint(ule,1,ie(nie-8,ma),1,'ULE')
        endif
      endif ! elflg = true

c     For each dof compute perturbed residual

      do i = 1,nen
        do j = 1,ndf
          temp1     = ul(j,i,1)
          temp2     = ul(j,i,2)
          temp3     = ul(j,i,3)
          temp4     = ul(j,i,4)
          temp5     = ul(j,i,5)
          ul(j,i,1) = temp1 + cc1*usize(j)
          ul(j,i,2) = temp2 + cc2*usize(j)
          if(np(42).ne.0) then
            ul(j,i,1) = temp1 + c3*usize(j)
            if(noi.eq.5) then
              ul(j,i,4) = temp4 + c4*usize(j)
              ul(j,i,5) = temp5 + c5*usize(j)
            endif
          endif

c         Compute perturbed residual

          call elmlib(hr(fp(3)),ul,hr(np(44)),ix,hr(np(39)),
     &                hr(np(36)),hr(np(35)),ndf,ndm,nst,iel,6)

c         Compute difference for tangent

          do k=1,nst
            st(k,ndf*(i-1)+j,1) = -hr(np(35)+k-1)
          end do

c         Set displacements back to original values

          ul(j,i,1)  = temp1
          ul(j,i,2)  = temp2
          ul(j,i,3)  = temp3
          ul(j,i,4)  = temp4
          ul(j,i,5)  = temp5

c         Compute perturbed displacements

          ul(j,i,1) = temp1 - cc1*usize(j)
          ul(j,i,2) = temp2 - cc2*usize(j)
          if(np(42).ne.0) then
            ul(j,i,1) = temp1 - c3*usize(j)
            if(noi.eq.5) then
              ul(j,i,4) = temp4 - c4*usize(j)
              ul(j,i,5) = temp5 - c5*usize(j)
            endif
          endif

c         Compute perturbed residual

          call elmlib(hr(fp(3)),ul,hr(np(44)),ix,hr(np(39)),
     &                hr(np(36)),hr(np(35)),ndf,ndm,nst,iel,6)

c         Compute difference for tangent

          do k=1,nst
            st(k,ndf*(i-1)+j,1) = (st(k,ndf*(i-1)+j,1)
     &                          +  hr(np(35)+k-1))*denom(j)
          end do

c         Set displacements back to original values

          ul(j,i,1)  = temp1
          ul(j,i,2)  = temp2
          ul(j,i,3)  = temp3
          ul(j,i,4)  = temp4
          ul(j,i,5)  = temp5
        end do ! j
      end do ! i

c     If element variables do perturbations

      if(elflg) then

c       Compute element residual

        denomel = 0.5d0/usizel
        mel     = ndf*nen
        do j = 1,ie(nie-8,ma)

c         Set value and perturbed state

          temp1  = ule(j)
          ule(j) = temp1 + cc1*usizel

c         Compute perturbed residual

          call elmlib(hr(fp(3)),ul,hr(np(44)),ix,hr(np(39)),
     &                hr(np(36)),hr(np(35)),ndf,ndm,nst,iel,6)

c         Set value for finite difference derivative

          do k=1,nst
            st(k,mel+j,1) = -hr(np(35)+k-1)
          end do ! k

c         Compute perturbed displacement

          ule(j) = temp1 - cc1*usizel

c         Compute perturbed residual

          call elmlib(hr(fp(3)),ul,hr(np(44)),ix,hr(np(39)),
     &                hr(np(36)),hr(np(35)),ndf,ndm,nst,iel,6)

c         Compute difference for tangent

          do k=1,nst
            st(k,mel+j,1) = (st(k,mel+j,1) +  hr(np(35)+k-1))*denomel
          end do ! k

c         Set value back to current state

          ule(j) = temp1

        end do ! j
      endif

c     Output numerical tangent

      if(flg) then

c       Compute size to print

        nrw = 0
        ncl = 0
        do n = 1,nst
          do i = 1,nst
            if(st(i,n,1).ne.0.0d0 .or. st(i,n,2).ne.0.0d0) then
              nrw = max(i,nrw)
              ncl = max(n,ncl)
            endif
          end do ! i
        end do ! n

c       If nearly full print all

        if(nrw+ndf.gt.nst) nrw = nst
        if(ncl+ndf.gt.nst) ncl = nst

        call mprint(st(1,1,2),nrw,ncl,nst,'S_tangent')
        call mprint(st(1,1,1),nrw,ncl,nst,'S_numerical')

c       Compute and output difference in tangent

        do n = 1,nst
          do i = 1,nst
            st(i,n,1) = st(i,n,2) - st(i,n,1)
          end do ! i
        end do ! n

        call mprint(st,nrw,ncl,nst,'Difference: S_t - S_num')


c     Replace tangent by numerical tangent

      else
        do m = 1,nst
          fp(1) = np(36) + (m-1)*nst - 1
          do i = 1,nst
            hr(fp(1)+i) = st(i,m,1)
          end do ! i
        end do ! m
      endif

      end
