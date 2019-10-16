c$Id:$
      subroutine mass3s(s,r,cfac,d,le,nst,ndm,ndf)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Frame mass matrix for small rotation angles

c     Inputs:
c       cfac       - Mass factor: 0 = lump; 1 = consistent.
c       d(*)       - Parameter set
c       le         - Element length
c       nst        - Array s dimension
c       ndm        - Mesh coordinate dimension (should be 3)
c       ndf        - Array r dimension

c     Outputs:
c       s(nst,nst) - Consistent mass array
c       r(ndf,*)   - Lumped mass array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   nst,ndm,ndf,i,ii,i1,j,jj,j1, ll
      real*8    cfac,lfac,le,lr,dv,s1,s2,s3
      real*8    d(*), r(ndf,*),s(nst,nst),sg(2,4),bb(2,2),db(2,2)
      real*8    nt(6,6,2),ntin(6,6),in(6,6)

      save

c     Initialize

      do i = 1,6
        do j = 1,6
          in(j,i)   = 0.0d0
          nt(j,i,1) = 0.0d0
          nt(j,i,2) = 0.0d0
        end do ! j
      end do ! i

c     Cross section given

      if(nint(d(100)).eq.0) then

c       Set inertial properties

        in(1,1)  = d(4)*d(32)
        in(2,2)  = in(1,1)
        in(3,3)  = in(1,1)
        in(4,4)  = d(4)*d(36)
        in(4,5)  = 0.0d0
        in(4,6)  = 0.0d0
        in(5,4)  = 0.0d0
        in(5,5)  = d(4)*d(33)*d(8)
        in(5,6)  = d(4)*d(35)*d(8)
        in(6,4)  = 0.0d0
        in(6,5)  = d(4)*d(35)*d(8)
        in(6,6)  = d(4)*d(34)*d(8)

c     Tubular cross section

      elseif(nint(d(100)).eq.1) then

        call b3tubm(d,in)

c     Rectangular cross section

      elseif(nint(d(100)).eq.2) then

        call b3rctm(d,in)

c     Section cross sections: wide(3);

      elseif(nint(d(100)).ge.3) then

        call b3secm(d,nint(d(100))-2, in)

      endif

c     Lumped mass matrix

      lr       = 1.d0/le

      do i = 1,ndm
        r(i,1) = 0.5d0*in(i,i)*le
        r(i,2) = r(i,1)
      end do ! i

c     Consistent mass matrix - Linear interpolation

      if(nint(d(79)).eq.0) then

        bb(1,1) = 0.0d0
        bb(1,2) = 0.0d0
        bb(2,1) = 0.0d0
        bb(2,2) = 0.0d0

        if(nint(d(182)).gt.0) then
          call int1dn(2, sg)
        else
          call int1d(2, sg)
        endif
        do ll = 1,2
          dv       = 0.5d0*le*sg(2,ll)
          s1      = 0.5d0 - 0.5d0*sg(1,ll)
          s2      = 0.5d0 + 0.5d0*sg(1,ll)
          bb(1,1) = bb(1,1) + s1*s1*dv
          bb(1,2) = bb(1,2) + s1*s2*dv
          bb(2,1) = bb(2,1) + s2*s1*dv
          bb(2,2) = bb(2,2) + s2*s2*dv
        end do ! ll

c       Assemble mass matrix

        i1 = 0
        do ii = 1,2
          j1 = 0
          do jj = 1,2
            do i = 1,6
              do j = 1,6
                s(i+i1,j+j1) = in(i,j)*bb(ii,jj)
              end do ! j
            end do ! i
            j1 = j1 + ndf
          end do ! jj
          i1 = i1 + ndf
        end do ! ii

c     Consistent mass matrix - Cubic interpolation

      else

        call int1d(4,sg)
        do ll = 1,4

          dv       = 0.5d0*le*sg(2,ll)

c         Axial and bending displacement functions

          s1       = 0.5d0 + 0.5d0*sg(1,ll)
          s2       = s1*s1
          s3       = s1*s2

          bb(1,2)  = 3.d0*s2 - s3 - s3
          bb(2,2)  = le*(s3 - s2)
          bb(1,1)  = 1.d0 - bb(1,2)
          bb(2,1)  = le*(s1 - s2) + bb(2,2)

          db(1,2)  = 6.d0*(s1 - s2)*lr
          db(2,2)  = 3.d0*s2 - 2.d0*s1
          db(1,1)  = -db(1,2)
          db(2,1)  = 1.d0 -2.d0*s1 + db(2,2)

          nt(1,1,1) = 1.d0 - s1
          nt(1,1,2) = s1
          nt(4,4,1) = 1.d0 - s1
          nt(4,4,2) = s1
          do ii = 1,2
            nt(2,2,ii) =  bb(1,ii)
            nt(2,6,ii) =  bb(2,ii)
            nt(3,3,ii) =  bb(1,ii)
            nt(3,5,ii) = -bb(2,ii)

            nt(5,3,ii) =  db(1,ii)
            nt(5,5,ii) = -db(2,ii)
            nt(6,2,ii) = -db(1,ii)
            nt(6,6,ii) = -db(2,ii)
          end do ! ii

c         Form consistent mass

          i1 = 0
          do ii = 1,2
            do i = 1,6
              do j = 1,6
                ntin(j,i) = dv*(nt(1,i,ii)*in(1,j)
     &                        + nt(2,i,ii)*in(2,j)
     &                        + nt(3,i,ii)*in(3,j)
     &                        + nt(4,i,ii)*in(4,j)
     &                        + nt(5,i,ii)*in(5,j)
     &                        + nt(6,i,ii)*in(6,j))
              end do ! j
            end do ! i

            j1 = 0
            do jj = 1,2

              do i = 1,6
                do j = 1,6
                  s(i+i1,j+j1) = s(i+i1,j+j1)
     &                         + ntin(1,i)*nt(1,j,jj)
     &                         + ntin(2,i)*nt(2,j,jj)
     &                         + ntin(3,i)*nt(3,j,jj)
     &                         + ntin(4,i)*nt(4,j,jj)
     &                         + ntin(5,i)*nt(5,j,jj)
     &                         + ntin(6,i)*nt(6,j,jj)
                end do ! j
              end do ! i
              j1 = j1 + ndf
            end do ! jj
            i1 = i1 + ndf
          end do ! ii
        end do ! ll

      endif

c     Interpolate mass between lumped and consistent
c     Consistent Mass: cfac = 1.0 ; Lumped Mass: cfac = 0.0

      lfac = 1.d0 - cfac
      do i = 1,nst
        do j = 1,nst
          s(i,j) = cfac*s(i,j)
        end do ! j
      end do ! i
      do i = 1,ndm
        s(i    ,i    ) = s(i    ,i    ) + lfac*r(i,1)
        s(i+ndf,i+ndf) = s(i+ndf,i+ndf) + lfac*r(i,2)
      end do ! i

      end
