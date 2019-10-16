c$Id:$
      subroutine eig3(v,d,rot)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute eigenvalues/vectors for 3x3 symmetric matrix

c      Inputs:
c         v(3,3) - matrix with initial values (only upper half used)

c      Outputs:
c         v(3,3) - matrix of eigenvectors (by column)
c         d(3)   - eigenvalues associated with columns of v
c         rot    - number of rotations to diagonalize
c-----[--.----+----.----+----.-----------------------------------------]
c     Storage done as follows:

c       | v(1,1) v(1,2) v(1,3) |     |  d(1)  a(1)  a(3)  |
c       | v(2,1) v(2,2) v(2,3) |  =  |  a(1)  d(2)  a(2)  |
c       | v(3,1) v(3,2) v(3,3) |     |  a(3)  a(2)  d(3)  |

c        Transformations performed on d(i) and a(i) and v(i,j) become
c        eigenvectors.  Thus, original array is destroyed.

c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   rot, its, i,j,k
      real*8    g,h, aij, sm,thresh, t, c,s,tau
      real*8    v(3,3), d(3), a(3), b(3), z(3)

      save

c     Move array into one-d arrays

      a(1) = v(1,2)
      a(2) = v(2,3)
      a(3) = v(1,3)

      do i = 1,3
        d(i) = v(i,i)
        b(i) = d(i)
        z(i) = 0.0d0
        do j = 1,3
          v(i,j) = 0.0d0
        end do ! j
        v(i,i) = 1.d0
      end do ! i

c     Check for diagonal case

      sm = abs(a(1)) + abs(a(2)) + abs(a(3))
      g  = abs(d(1)) + abs(d(2)) + abs(d(3))
      if (sm .lt. 1.d-13*g) return

      rot = 0
      do its = 1,50

c       Set convergence test and threshold

        sm = abs(a(1)) + abs(a(2)) + abs(a(3))
        if (sm.eq.0.0d0) return

        if(its.lt.4) then
          thresh = 0.011d0*sm
        else
          thresh = 0.0d0
        end if

c       Perform sweeps for rotations

        do i = 1,3
          j = mod(i,3) + 1
          k = mod(j,3) + 1

          aij  = a(i)
          g    = 100.d0*abs(aij)
          if(abs(d(i))+g.ne.abs(d(i)) .or.
     &       abs(d(j))+g.ne.abs(d(j))) then

            if(abs(aij).gt.thresh) then
              a(i) = 0.0d0
              h    = d(j) - d(i)
              if(abs(h)+g.eq.abs(h)) then
                t = aij/h
              else
                t = sign(2.d0,h/aij)/(abs(h/aij)+sqrt(4.d0+(h/aij)**2))
              endif

c             Set rotation parameters

              c    = 1.d0/sqrt(1.d0+t*t)
              s    = t*c
              tau  = s/(1.d0+c)

c             Rotate diagonal terms

              h    = t*aij
              z(i) = z(i) - h
              z(j) = z(j) + h
              d(i) = d(i) - h
              d(j) = d(j) + h

c             Rotate off-diagonal terms

              h    = a(j)
              g    = a(k)
              a(j) = h + s*(g - h*tau)
              a(k) = g - s*(h + g*tau)

c             Rotate eigenvectors

              do k = 1,3
                g      = v(k,i)
                h      = v(k,j)
                v(k,i) = g - s*(h + g*tau)
                v(k,j) = h + s*(g - h*tau)
              end do ! k

              rot = rot + 1

            end if
          else
            a(i) = 0.0d0
          end if
        end do ! i

c       Update diagonal terms

        do i = 1,3
          b(i) = b(i) + z(i)
          d(i) = b(i)
          z(i) = 0.0d0
        end do ! i

      end do ! its

      end
