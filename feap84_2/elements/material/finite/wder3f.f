c$Id:$
      subroutine wder3f(d,bpr, tautil,atilp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use constants from 'pconstant.h'                   14/11/2006
c     2. Correct comment on definition of w_i               07/06/2007
c-----[--.----+----.----+----.-----------------------------------------]

c     Stored Energy Function in Principal Stretches
c          W      = w(lam_1) + w(lam_2) + w(lam_3)

c          lam _i = (Je)^(-1/3)*lambda_i (deviatoric stretches)

c     Ogden strain energy functions in terms of Seth strains
c          w_i  = c-alpha/(m-alpha)*(lam_i**m-alpha - 1.) -> wengy

c     Inputs:
c          d(22) = first  ci
c          d(23) = first  ni
c          d(24) = second ci
c          d(25) = second ni
c          d(26) = third  ci
c          d(27) = third  ni
c          d(28) = no. terms in series expansion of W(lam) = 3, max

c          bpr(3) principal values of left Cauchy-Green tensor (-1)

c     Outputs:
c          tautil(3)   Principal deviatoric Kirchoff stresses
c          atilp(6,6)  Deviatoric tangent matrix in principal basis
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdamag.h'
      include  'pconstant.h'

      integer   nn, i, j, k, n
      real*8    jthrd, c1, c2, c3, c4, tol, trtautil, w
      real*8    d(*),bpr(3),lamt(3),lam1(3),lam2(3),tautil(3)
      real*8    atilp(6,6),dwtil(4,3),dtdl(3,3)

      real*8    a1,a2,a3,a4,a5,a6

      save

      data      tol    /  1.0d-06 /

c     Initialize arrays

      do i = 1,3
        tautil(i) = 0.0d0
        do j = 1,3
          dtdl(j,i) = 0.0d0
        end do ! j
        do j = 1,4
          dwtil(j,i) = 0.0d0
        end do ! j
      end do ! i
      do i = 1,6
        do j = 1,6
          atilp(j,i) = 0.0d0
        end do ! j
      end do ! i

c     Compute J**(-1/6): w = J - 1 (Use series for small)

      w = bpr(1) + bpr(2) + bpr(3)
     &  + bpr(1) * bpr(2) + bpr(2) * bpr(3) + bpr(3) * bpr(1)
     &  + bpr(1) * bpr(2) * bpr(3)
      if(w.gt.0.001d0) then
        jthrd = (1.d0 + w)**(-one6) - 1.0d0
      else
        a1 =-one6
        a2 = one2*(a1-1.0d0)
        a3 = a2*one3*(a1-2.0d0)
        a4 = a3*one4*(a1-3.0d0)
        a5 = a4*one5*(a1-4.0d0)
        a6 = a5*one6*(a1-5.0d0)
        jthrd = a1*(w + a2*w**2 + a3*w**3 + a4*w**4
     &                + a5*w**5 + a6*w**6)
      endif

c     Compute modified (Flory) stretches - 1.0 (Use series for small)

      do i = 1,3
        if(abs(bpr(i)).gt.0.001d0) then
          lamt(i) = (1.d0+jthrd)*sqrt(1.d0+bpr(i)) - 1.d0
        else
          w  = bpr(i)
          a1 = one2
          a2 = one2*(a1-1.0d0)
          a3 = a2*one3*(a1-2.0d0)
          a4 = a3*one4*(a1-3.0d0)
          a5 = a4*one5*(a1-4.0d0)
          a6 = a5*one6*(a1-5.0d0)
          lam1(i) = a1*(w + a2*w**2 + a3*w**3 + a4*w**4
     &                    + a5*w**5 + a6*w**6)
          lamt(i) = jthrd + lam1(i) + jthrd * lam1(i)
        endif
        lam1(i) = 1.d0 + lamt(i)
        lam2(i) = lam1(i)**2
      end do ! i

c     Compute strain energy function derivatives

      wengy = 0.0d0
      nn    = nint(d(28))

      do n = 1,nn
        c1   = d(20+2*n)
        c2   = d(20+2*n+1)
        do i = 1,3

c         Compute Kirchhoff principal stresses (Use series for small)

          if(abs(lamt(i)).gt.0.001d0) then
            tautil(i) = tautil(i) + c1*(lam1(i)**c2 - 1.d0)
          else
            w  = lamt(i)
            a1 = c2
            a2 = one2*(a1-1.0d0)
            a3 = a2*one3*(a1-2.0d0)
            a4 = a3*one4*(a1-3.0d0)
            a5 = a4*one5*(a1-4.0d0)
            a6 = a5*one6*(a1-5.0d0)
            tautil(i)  = tautil(i)
     &                 + c1*a1*(w + a2*w**2 + a3*w**3 + a4*w**4
     &                            + a5*w**5 + a6*w**6)
          endif

c         Compute energy and higher derivatives for tangent

          wengy      = wengy      + (c1/c2)*(lam1(i)**c2 - 1.d0)
          dwtil(1,i) = dwtil(1,i) + c1*lam1(i)**(c2-1.d0)
          dwtil(2,i) = dwtil(2,i) + (c2-1.d0)*c1*lam1(i)**(c2-2.d0)
          dwtil(3,i) = dwtil(3,i)
     &                  + (c2-2.d0)*(c2-1.d0)*c1*lam1(i)**(c2-3.d0)
          dwtil(4,i) = dwtil(4,i)
     &        + (c2-3.d0)*(c2-2.d0)*(c2-1.d0)*c1*lam1(i)**(c2-4.d0)
        end do ! i
      end do ! n

c     Compute trace stress term

      trtautil = (tautil(1) + tautil(2) + tautil(3))*one3

c     Compute deviatoric principal Kirchoff stresses

      do i=1,3
        tautil(i) = tautil(i) - trtautil
      end do ! i

c     Compute deviatoric tangent matrix in principal basis

      do i=1,3
        dtdl(i,i) = dwtil(2,i)*lam1(i) + dwtil(1,i)
        do j = 1,3
          dtdl(i,j) = dtdl(i,j)
     &              - ( dwtil(2,j)*lam1(j) + dwtil(1,j) )*one3
        end do ! j
      end do ! i

c     Compute deviatoric tangent matrix in principal basis

      do i = 1,3
        k = 1+mod(i,3)

        if( abs(lamt(i)-lamt(k)).gt.tol )then
          atilp(i+3,i+3) =
     &         ( lam2(i)*tautil(k) - lam2(k)*tautil(i) )/
     &         ( lam2(k) - lam2(i) )
        else
          c3=1.d0/(2.d0 + lamt(i) + lamt(k))
          c4=lamt(k) - lamt(i)
          atilp(i+3,i+3) = c3*lam1(i)*lam1(k) * ( - dwtil(1,i)
     &                   + dwtil(2,i)*lam1(i)
     &                   + dwtil(3,i)*lam1(i)*c4*0.5d0
     &                   + dwtil(4,i)*lam1(i)*c4**2*one6 )
     &                 + ( dwtil(1,1)*lam1(1)
     &                   + dwtil(1,2)*lam1(2)
     &                   + dwtil(1,3)*lam1(3) )*one3
        endif

        atilp(i,i) = -2.d0*tautil(i)
        do j = 1,3
          atilp(i,j) = atilp(i,j) + lam1(j)*dtdl(i,j)
     &                          - ( lam1(1)*dtdl(i,1)
     &                            + lam1(2)*dtdl(i,2)
     &                            + lam1(3)*dtdl(i,3) )*one3
        end do ! j
      end do ! i

      end
