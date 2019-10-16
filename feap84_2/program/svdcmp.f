c$Id:$
      subroutine svdcmp(u,w,v,t, mp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Singular valued decomposion.
c        Note: A = V'WU & U replace A for output. (N.B., A stored in U)
c        This form of the algorithm considers square arrays only using
c        full memory allocation for each array.

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i,its, j,jj, k, l, mp
      real*8    c, f, g, h, s, x, y, z, scale, anorm
      real*8    u(mp,*), w(*), v(mp,*), t(*)

      save

c     Householder tridiagonal form reduction

      g     = 0.d0
      scale = 0.d0
      anorm = 0.d0

      do i = 1,mp
        l     = i + 1
        t(i)  = scale*g
        g     = 0.d0
        s     = 0.d0

        scale = 0.d0
        do j = i,mp
          scale = scale + abs(u(j,i))
        end do ! j

        if(scale.ne.0.d0) then
          f = 1.d0/scale
          do j = i,mp
            u(j,i) = u(j,i)*f
            s      = s + u(j,i)*u(j,i)
          end do ! j

          f      = u(i,i)
          g      = -sign(sqrt(s),f)
          h      = f*g - s
          u(i,i) = f - g

          if(i.ne.mp) then
            do j = l,mp
              s = 0.d0
              do k = i,mp
                s = s + u(k,i)*u(k,j)
              end do ! k
              f = s/h
              do k = i,mp
                u(k,j) = u(k,j) + u(k,i)*f
              end do ! k
            end do ! j
          end if

          do j = i,mp
            u(j,i) = u(j,i)*scale
          end do ! j

        end if

        w(i)  = scale*g
        g     = 0.d0
        s     = 0.d0
        scale = 0.d0

        if(i.ne.mp) then

          do k = l,mp
            scale = scale + abs(u(i,k))
          end do ! k

          if(scale.ne.0.d0) then

            f = 1.d0/scale
            do j = l,mp
              u(i,j) = u(i,j)*f
              s      = s + u(i,j)*u(i,j)
            end do ! j

            f      =  u(i,l)
            g      = -sign(sqrt(s),f)
            h      =  f*g - s
            u(i,l) =  f - g

            f = 1.d0/h
            do j = l,mp
              t(j) = u(i,j)*f
            end do ! j

            do j = l,mp
              s = 0.d0
              do k = l,mp
                s = s + u(j,k)*u(i,k)
              end do ! k
              do k = l,mp
                u(j,k) = u(j,k) + t(k)*s
              end do ! k
            end do ! j

            do j = l,mp
              u(i,j) = u(i,j)*scale
            end do ! j

          end if

        end if

        anorm = max(anorm, (abs(w(i)) + abs(t(i))))

      end do ! i

c     Right side transformations

      do i = mp,1,-1

        if(g.ne.0.d0) then

          f =  1.d0/g
          do j = l,mp
            v(j,i) = (u(i,j)/u(i,l))*f
          end do ! j

          do j = l,mp

            s = 0.d0
            do k = l,mp
              s = s + u(i,k)*v(k,j)
            end do ! k

            do k = l,mp
              v(k,j) = v(k,j) + v(k,i)*s
            end do ! k

          end do ! j

        end if

        do j = l,mp
          v(i,j) = 0.d0
          v(j,i) = 0.d0
        end do ! j

        v(i,i) = 1.d0
        g      = t(i)
        l      = i
      end do ! i

c     Left side transformations

      do i = mp,1,-1
        l = i + 1
        g = w(i)

        do j = l,mp
          u(i,j) = 0.d0
        end do ! j

        if(g.ne.0.d0) then
          g = 1.d0/g
          do j = l,mp
            s = 0.d0
            do k = l,mp
              s = s + u(k,i)*u(k,j)
            end do ! k
            f = (s/u(i,i))*g
            do k = i,mp
              u(k,j) = u(k,j) + u(k,i)*f
            end do ! k
          end do ! j

          do j = i,mp
            u(j,i) = u(j,i)*g
          end do ! j

        else

          do j = i,mp
            u(j,i) = 0.d0
          end do ! j

        end if

        u(i,i) = u(i,i) + 1.d0

      end do ! i

c     Diagonalization step

      do k = mp,1,-1
        do its = 1,30

          do l = k,1,-1
            if((abs(t(l  ))+anorm).eq.anorm) go to 2  ! N.B. t(1) = 0.0
            if((abs(w(l-1))+anorm).eq.anorm) go to 1
          end do ! l

c         Reduce t(l) for l > 1

1         c = 0.d0
          s = 1.d0

          do i = l,k
            f = s*t(i)
            if((abs(f)+anorm).ne.anorm) then
              g    = w(i)
              h    = sqrt(f*f + g*g)
              w(i) = h
              h    = 1.d0/h
              c    =  g*h
              s    = -f*h
              do j = 1,mp
                z        =  (u(j,l-1)*c) + (u(j,i)*s)
                u(j,i  ) = -(u(j,l-1)*s) + (u(j,i)*c)
                u(j,l-1) =   z
              end do ! j
            end if
          end do ! i

2         z = w(k)

          if(l.eq.k) then

            if(z.lt.0.d0) then
              w(k) = -z
              do j = 1,mp
                v(j,k) = -v(j,k)
              end do ! j
            end if
            go to 3

          end if

          x  = w(l)
          y  = w(k-1)
          g  = t(k-1)
          h  = t(k)
          f  = ((y-z)*(y+z) + (g-h)*(g+h))/(2.d0*h*y)
          g  = sqrt(f*f +1.d0)
          f  = ((x-z)*(x+z) + h*((y/(f+sign(g,f))) - h))/x
          c  = 1.d0
          s  = 1.d0
          do j = l,k-1
            i    =  j + 1
            g    =  t(j+1)
            y    =  w(j+1)
            h    =  s*g
            g    =  c*g
            z    =  sqrt(f*f + h*h)
            t(j) =  z
            c    =  f/z
            s    =  h/z
            f    =  x*c + g*s
            g    = -x*s + g*c
            h    =  y*s
            y    =  y*c
            do jj = 1,mp
              z         =  (v(jj,j)*c) + (v(jj,j+1)*s)
              v(jj,j+1) = -(v(jj,j)*s) + (v(jj,j+1)*c)
              v(jj,j  ) =   z
            end do ! jj
            z    = sqrt(f*f + h*h)
            w(j) = z
            if(z.ne.0.d0) then
              z = 1.d0/z
              c = f*z
              s = h*z
            end if
            f =  c*g + s*y
            x = -s*g + c*y
            do jj = 1,mp
              z         =  (u(jj,j)*c) + (u(jj,j+1)*s)
              u(jj,j+1) = -(u(jj,j)*s) + (u(jj,j+1)*c)
              u(jj,j  ) =   z
            end do ! jj
          end do ! j
          t(l) = 0.d0
          t(k) = f
          w(k) = x
        end do ! its
        write(*,*) ' *WARNING* No convergence in 30 iterations.'

3       continue
      end do ! k

      end
