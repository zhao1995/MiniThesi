c$Id:$
      subroutine pnodal(id,cstif,cdamp,cmass,pmass,ad,dr,u,ud,nnel,
     &                  consist,lumped,react,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add 'DTANG' for use in sparse matrix assembly    21/03/2013
c          Assemble nodal values into correct diagonals.
c       2. Correct set of mass value for parallel solution  31/03/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Assemble nodal stiffness, damping, and mass effects

c      Inputs:
c         id(*)     - Equation numbers of active dof
c         cstif(*)  - Stiffness coefficient for each dof
c         cdamp(*)  - Damping   coefficient for each dof
c         cmass(*)  - Lump mass coefficient for each dof
c         pmass(*)  - Lump mass prop loads  for each dof
c         u(*)      - Solution values
c         ud(*)     - Rate values
c         nnel      - Number terms in id
c         consist   - Flag, assemble consistent mass diagonal if true
c         lumped    - Flag, assemble lumped mass diagonal if true
c         react     - Flag, assemble uncompressed residual if true
c         isw       - Switch on action to perform

c      Outputs:
c         ad(*)     - Diagonals of tangent array
c         dr(*)     - Residual array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'compac.h'
      include  'compas.h'
      include  'pointer.h'
      include  'comblk.h'

      include  'ddata.h'
      include  'eltran.h'
      include  'fdata.h'
      include  'pfeapb.h'
      include  'prld1.h'
      include  'ptdat6.h'
      include  'rdat1.h'
      include  'refng.h'
      include  'sdata.h'

      include  'iofile.h'

      logical   consist,lumped,react, setvar,palloc
      integer   e,i,n,nnel,nrhs, isw
      integer   id(nnel),pmass(nnel)
      real*8    cstif(nnel),cdamp(nnel),cmass(nnel)
      real*8    ad(*),dr(*),u(*),ud(nnel,*), mom, dres

      save

c     Check for parallel or serial solution to set number of equations

      if(pfeap_on) then
        nrhs = ndf*numpn
      else
        nrhs = nnel
      endif

c     Modify tangent: diagonal only

      if(isw.eq.3) then

c       Parallel assembly

        if(pfeap_on) then
          setvar = palloc(320,'DTANG',nrhs,2)
          do n = 1,nrhs
            hr(np(320)+n-1) = cstif(n)*(ctan(1)+gray(2)*ctan(2))
     &                      + cdamp(n)* ctan(2)
     &                      + cmass(n)*(ctan(3)+gray(1)*ctan(2))
          end do ! n
          call uassdk(hr(np(320)),mr(np(245)))
          setvar = palloc(320,'DTANG',   0,2)

c       Serial assembly

        else
          if(compfl) then
            setvar = palloc(320,'DTANG',nrhs,2)
            do n = 1,nrhs
              hr(np(320)+n-1) = cstif(n)*(ctan(1)+gray(2)*ctan(2))
     &                        + cdamp(n)* ctan(2)
     &                        + cmass(n)*(ctan(3)+gray(1)*ctan(2))
            end do ! n
            call cassdi(ad,hr(np(320)),mr(np(94)),mr(np(93)),id,nrhs,
     &                  kbycol,kdiag,kall)
            setvar = palloc(320,'DTANG',   0,2)
          else
            do n = 1,nrhs
              if(id(n).gt.0) then
                ad(id(n)) = ad(id(n))
     &                    + cstif(n)*(ctan(1)+gray(2)*ctan(2))
     &                    + cdamp(n)* ctan(2)
     &                    + cmass(n)*(ctan(3)+gray(1)*ctan(2))
              endif
            end do ! n
          endif

        endif

c     Add point mass terms

      elseif(isw.eq.5) then

c       Parallel assembly

        if(pfeap_on) then
          setvar = palloc(320,'DTANG',nrhs,2)
          do n = 1,nrhs
            hr(np(320)+n-1) = cmass(n)
          end do ! n
          call uassdm(hr(np(320)),mr(np(245)),consist,lumped)
          setvar = palloc(320,'DTANG',   0,2)

c       Serial assembly

        else

          if(compfl) then
            setvar = palloc(320,'DTANG',nrhs,2)
            do n = 1,nrhs
              hr(np(320)+n-1) = cmass(n)
            end do ! n
            call cassdi(ad,hr(np(320)),mr(np(91)),mr(np(90)),id,nrhs,
     &                  mbycol,mdiag,mall)
            setvar = palloc(320,'DTANG',   0,2)
          else
            do n = 1,nrhs
              if(id(n).gt.0) then

c               Consistent mass

                if(consist) ad(id(n)) = ad(id(n)) + cmass(n)

c               Lumped mass

                if(lumped ) dr(id(n)) = dr(id(n)) + cmass(n)

              endif
            end do ! n
          endif
        endif

c     Add diagonal damper terms

      elseif(isw.eq.9) then

        if(pfeap_on) then
          if(ior.lt.0) then
            write(*,*) ' * * ERROR * * No Parallel Point Dampers'
          endif
          write(iow,*) ' * * ERROR * * No Parallel Point Dampers'
          call plstop()
        elseif(compfl) then
          setvar = palloc(320,'DTANG',nrhs,2)
          do n = 1,nrhs
            hr(np(320)+n-1) = cdamp(n)
          end do ! n
          call cassdi(ad,hr(np(320)),mr(np(204)),mr(np(203)),id,nrhs,
     &                dbycol,ddiag,dall)
          setvar = palloc(320,'DTANG',   0,2)
        else
          do n = 1,nrhs
            if(id(n).gt.0) then
              ad(id(n)) = ad(id(n)) + cdamp(n)
            endif
        end do ! n
        endif

c     Compute momentum and energy

      elseif(isw.eq.13) then

        do n = 1,nrhs,ndf
          do i = 0,ndm-1
            mom    = cmass(n+i)*ud(n+i,1)
            epl(i) = epl(i) + mom
            epl(7) = epl(7) + 0.5d0*mom*ud(n+i,1)
          end do ! i
        end do ! n

      endif

c     Correct residual for stiffness (damping/mass of transient problem)

      if(isw.eq.3 .or. isw.eq.6) then

        do n = 1,nrhs

c         Lumped stiffness: linear spring

          if(react) then           ! Reaction outputs
            e = n
          else                     ! Residual for solutions
            e = id(n)
          endif

c         Compute contribution for residual or reaction

          if(e.gt.0 .and. cstif(n).ne.0.0d0) then
            if(nrk.gt.0) then
              dres = cstif(n)*ud(n,nrk)
            else
              dres = cstif(n)*u(n)
            endif

c           Assemble residual or reaction

            dr(e) = dr(e) - dres

c           Compute residual norm for convergence check

            if(compre) then
              rnorm1 = rnorm1 + abs(dres)
              rnormn = rnormn + 1.d0
            endif

          endif

        end do ! n

        if(fl(9)) then

          do n = 1,nrhs

            if(react) then           ! Reaction outputs
              e = n
            else                     ! Residual for solutions
              e = id(n)
            endif

c           Lumped damping: linear damper

            if(nrc.gt.0) then
              dr(e) = dr(e) - cdamp(n)*ud(n,nrc)
            endif

c           Lumped inertia: point mass

            if(nrm.gt.0) then
              dr(e) = dr(e) - cmass(n)*ud(n,nrm)
            endif

          end do ! n

c         Modify residual for proportional acceleration terms

          do i = 1,ndf
            if(gfac(i).ne.0.0d0) then
              do n = i,nrhs,ndf
                if(pmass(n).gt.0) then
                  if(react) then
                    dr(n) = dr(n) - cmass(n)*prldv(pmass(n))*gfac(i)
                  elseif(id(n).gt.0) then
                    dr(id(n)) = dr(id(n))
     &                        - cmass(n)*prldv(pmass(n))*gfac(i)
                  endif
                endif
              end do ! n
            endif ! gfac
          end do ! i

c         Rayleigh damping residual modifications

          if(max(abs(gray(1)),abs(gray(2))).gt.0.0d0) then
            do n = 1,nrhs
              if(nrc.gt.0) then
                if(react) then
                  dr(n) = dr(n) - (gray(1)*cmass(n)
     &                  + gray(2)*cstif(n))*ud(n,nrc)
                elseif(id(n).gt.0) then
                  dr(id(n)) = dr(id(n)) - (gray(1)*cmass(n)
     &                      + gray(2)*cstif(n))*ud(n,nrc)
                endif
              endif
            end do ! n
          endif

        endif

      end if

      end
