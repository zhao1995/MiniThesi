c----------------------------------------------------------------------+
c     Interface  to use SuperLU solver in FEAP                         |
c----------------------------------------------------------------------+
c
      subroutine mkptr3(isymcsr)
c-----------------------------------------------------------------------
c
c     Purpose: Define arrays for SuperLU solver
c              pointer with CSR storage in mkptr_csr
c
c     Inputs: -
c
c     Output: isymcsr  ! 1= symmetric 2= unsymmetric matrix
c             idrslu   array necessary for b and u
c
c-----------------------------------------------------------------------
      USE cdata
      USE pdata2
      USE slu
      USE doalloc
       
      implicit double precision (a-h,o-z)
      logical ldummy

      if(idev.ne.4) stop 'solver only possible for idev=4 (SALFORD)'

      isymcsr  = 2
      ifactors = 0

      call ralloc(  drslu,neq,'SUPERLU-dr  ',ldummy)
      call ralloc(diagslu,neq,'SUPERLU-diag',ldummy)

      return
      end
c
c      subroutine datri3(a,ia,neq,nneg) in SR PLOT_SAL_CW

      subroutine datri3nd(ad,neq,nneg)
c----------------------------------------------------------------------
c
c      Purpose: Calculate number of negative diagonals
c               Print/Plot data
c      Inputs:
c         ad(*)  - Factored diagonal terms
c
c      Outputs:
c         nneg   - number of negativ diagonal entries
c
c----------------------------------------------------------------------
      USE dii
      implicit double precision (a-h,o-z)
      dimension ad(*)
c...  initial values
      zero = 0.0d0
      nneg = 0
      call pzeroi(ndii,150)
      call pzeroi(ii,3)

c.... loop through factorized diagonal entries
      do i = 1,neq
        adi=ad(i)
        if(abs(adi).gt.1.d-30) then
          if(ad(i).lt.zero) then
            nneg=nneg+1       ! count
            call prtdii1(i,1) ! plot/print
          end if
        end if
      end do

      return
      end
c
c     subroutine dasol3 (a,b,ia,neq,energy) in SR PLOT_SAL_CW.for
c
      subroutine detkt3(neq,nneg,det1)
c----------------------------------------------------------------------
c     Purpose: calculate determinant of stiffness matrix
c              for SuperLU solver: value det1
c      Inputs:
c         m(idrslu)  - Factored diagonal terms of A-array
c         neq        - Length of A array
c
c      Outputs:
c        det1 - value of exp(determinant)
c
c----------------------------------------------------------------------
      USE slu 
      implicit double precision(a-h,o-z)

      call detkt31(diagslu,neq,nneg,det1)
      return
      end
c
      subroutine detkt31(ad,neq,nneg,det1)
c----------------------------------------------------------------------
c
c     Purpose: calculate value det1 of determinant of stiffness matrix
c              called from SuperLU solver
c
c      Inputs:
c         ad(neq)  - Factored diagonal terms of A-array
c         neq      - Length of A array
c
c     Outputs:
c        det1 - value of exp(determinant)
c
c----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension ad(*)
      zero = 0.0d0
c.... calculate determinant
      do i = 1,neq
        adi=ad(i)
        if(abs(adi).gt.1.d-30) then
          if(adi.lt.zero) then
            xx = -adi
            nneg = nneg + 1
          else
            xx = adi
          end if
          det1 = det1 + log(xx)
        end if
      end do

      return
      end
c


