c$Id:$
      subroutine pfiber(d, a0, f, j, sig, dd, ftype)

c-----[--+---------+---------+---------+---------+---------+---------+-]
c      Purpose: Fiber models for finite deformation constitution
c               N.B. - Use I_4 not \bar{I_4}
c               1.) Holzapfel-Gasser model

c      Psi(I_4) = 0.5*h1/h2*[exp(h2*(I_4 - 1)**2) - 1)]
c      I_4      = A_I*C_IJ*A_J = (F_iI*A_I)*(F_iJ*A_J) = a_i*a_i
c      a_i      = F_iI*A_I

c      Inputs:
c        d(*)    - Material parameters
c        a0(3)   - Structure vector
c        f(3,3)- Deformation and displacement gradient
c        j       - Jacobian determinant
c        ftype   = 1 uses F; else uses G = displacement gradeint

c      Outputs:
c        sig(6)  - Cauchy stress
c        dd(6,6) - Moduli (current configuration)
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      real*8     d(*), a0(3), f(3,3), j, sig(6), dd(6,6)

      integer    ftype, i
      real*8     a(3), tau(6), i4, ee, h1,h2, dpsi1, dpsi2

c     Compute current structure tensor

      if(ftype.eq.1) then             ! Use deformation gradient
        a(:) = 0.0d0
        do i = 1,3
          a(:) = a(:) + f(:,i)*a0(i)
        end do ! i
      else                           ! Use displacement gradient
        a(:) = a0(:)
        do i = 1,3
          a(:) = a(:) + f(:,i)*a0(i)
        end do ! i
      endif

c     Compute I4

      i4 = a(1)**2 + a(2)**2 + a(3)**2

c     Stress transformation to Voigt notation
      tau(1) = a(1)**2
      tau(2) = a(2)**2
      tau(3) = a(3)**2
      tau(4) = a(1)*a(2)*2.d0
      tau(5) = a(2)*a(3)*2.d0
      tau(6) = a(3)*a(1)*2.d0

c     Compute remaining terms

      if(i4.ge.1.d0) then

        if(nint(d(1)).eq.1) then ! Holzapfel-Gasser model

c         Factors in stored energy function
          i4 = i4 - 1.0d0
          h1 = d(2)
          h2 = d(3)
          ee = h1*exp(h2*i4**2)

c         Compute first derivative of stored energy function w/r I_4
          dpsi1 = ee*i4/j

c         Compute second derivative of stored energy function w/r I_4
          dpsi2 = ee*(1.d0 + 2.d0*h2*i4**2)/j

c       Weiss model

        elseif(nint(d(1)).eq.2) then

c         Factors in stored energy function
          h1 = d(2)
          h2 = d(3)
          ee = h1*exp(h2*(i4-1.d0))

c         Compute first derivative of stored energy function w/r I_4
          dpsi1 = (ee - h1*i4**(h2-1.0d0))/j

c         Compute second derivative of stored energy function w/r I_4
          dpsi2 = (ee*h2 - h1*(h2-1.0d0)*i4**(h2-2.0d0))/j

c       No model
        else

          dpsi1 = 0.0d0
          dpsi2 = 0.0d0

        endif

c       Multiply by 2 (S_IJ = 2 * dPsi/C_IJ)
        tau    = 2.d0*tau

c       Compute tangent moduli
        do i = 1,6
          dd(:,i) = dd(:,i) + tau(:)*tau(i)*dpsi2
        end do ! i

c       Compute stress
        sig = sig + tau*dpsi1

      endif

      end
