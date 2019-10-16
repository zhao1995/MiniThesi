c$Id:$
      subroutine pqrblk(op,oq,or, np1,nq1,nr1, ne,ni, ma, ix,nen1)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    06/02/2009
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Generate a block of p x q x r Lagrangian elements

c      Inputs:
c         op         - Order of element in p-direction
c         oq         - Order of element in q-direction
c         or         - Order of element in r-direction
c         np1        - Number nodes in r-direction
c         nq1        - Number nodes in s-direction
c         nr1        - Number nodes in t-direction
c         ne         - Initial element number
c         ni         - Initial node number
c         ma         - Material number of block

c      Outputs:
c        ix(nen1,*)  - Element connection list
c        ne          - Final element number
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    op,oq,or, np,nq,nr, ne,ni, ma, nen1
      integer    ix(nen1,*)

      integer    i,j,k, ii,jj,kk, nii,njj,nkk, ij
      integer    ne1,ne2,ne3, np1,nq1,nr1

c     Generate elements pxqxr type

      np  = np1 - 1
      nq  = nq1 - 1
      nr  = nr1 - 1

      ne1 = np/op
      ne2 = nq/oq
      ne3 = nr/or

      nii = ni - 1
      do k = 1,ne3
        do j = 1,ne2
          do i = 1,ne1
            njj = np1*(nq1*or*(k-1) + oq*(j-1)) + op*(i-1) + nii
            ij  = 0
            do kk = 1,or+1
              nkk = njj
              do jj = 1,oq+1
                do ii = 1,op+1
                  ix(ii+ij,ne) = nkk + ii
                end do ! ii
                ij = ij + op + 1
                nkk = nkk + np1
              end do ! jj
              njj = njj + np1*nq1
            end do ! kk
            ix(nen1  ,ne) = ma   ! Material number
            ne = ne + 1
          end do ! i
        end do ! j
      end do ! k
      ne = ne - 1

      end
