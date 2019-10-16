c$Id:$
      subroutine dsolmod(afd,afl,afu,wf,mf,nc,aflg,bflg,cflg,nmod)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set up modal solution for rigid body options

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'ddata.h'
      include  'modreg.h'
      include  'pointer.h'
      include  'tdata.h'
      include  'tdatb.h'
      include  'comblk.h'

      include  'p_int.h'

      logical   aflg, bflg, cflg
      integer   mf, nc, nmod, i,j,n
      real*8    cn3,cn4,cn5,cn6,ur1,ur2,drot(3)
      real*8    afd(*),afl(mf,*),afu(mf,*),wf(mf,*)

      save

c     Construct reduced tangent matrix

      if(aflg) then

c       Initialize arrays

        do i = 1,nc
          do j = 1,nc
            ar(j,i) = 0.0d0
          end do ! j
          br(i) = 0.0d0
        end do ! i

c       Reduce tangent vector

        do n = 1,mf
          afd(n) = 1.d0/afd(n)
          do j = 1,nc
            afu(n,j) = afd(n)*afu(n,j)
            do i = 1,nc
              ar(i,j) = ar(i,j) - afl(n,i)*afu(n,j)
            end do ! i
          end do ! j
        end do ! n

c       Reduce force vector

        do n = 1,mf
          wf(n,3) = afd(n)*wf(n,3)
          do i = 1,nc
            br(i) = br(i) - afl(n,i)*wf(n,3)
          end do ! i
        end do ! n

      endif

c     Do the initial solution for modal values
c     N.B., Energy-Momentum solution option only.

      if(bflg) then

        cn3 = 1.d0 - 0.5d0/theta(1)
        cn4 = 1.d0 - theta(2)/theta(1)
        cn5 = dt*(1.d0 - 0.5d0*theta(2)/theta(1))
        cn6 = 1.d0/(theta(1)*dt)

        do n = 1,mf
          ur1     =  cn4*wf(n,4) + cn5*wf(n,5)
          ur2     = -cn6*wf(n,4) + cn3*wf(n,5)
          wf(n,6) =  wf(n,1)
          wf(n,7) = (1.d0-theta(3))*wf(n,4) + theta(3)*ur1
          wf(n,8) = (ur1-wf(n,4))/dt
          wf(n,4) =  ur1
          wf(n,5) =  ur2
          wf(n,2) =  0.d0
        end do ! n
      endif

c     Do the final solution for modal values

      if(cflg) then

c       Compute the increment to modes

        fp(1) = np(108) + 6*(modbod(nmod) - 1) + 2
        do i = 1,3
          drot(i) = hr(fp(1)+i)
        end do ! i
        do n = 1,mf
          do i = 1,3
            wf(n,3) = wf(n,3) - afu(n,i)*drot(i)
          end do ! i
        end do ! n

c       Update the solution vectors

        do n = 1,mf
          wf(n,1) = wf(n,1) + cc1*wf(n,3) ! d_n+1
          wf(n,2) = wf(n,2) + cc1*wf(n,3) ! DELTA d_n+1

c         Update the mid-point values and n+1 rates

          wf(n,4) = wf(n,4) + c2*wf(n,3)  ! v_n+1
          wf(n,5) = wf(n,5) + c1*wf(n,3)  ! a_n+1
          wf(n,6) = wf(n,6) + c3*wf(n,3)  ! d_n+a
          wf(n,7) = wf(n,7) + c4*wf(n,3)  ! v_n+a
          wf(n,8) = wf(n,8) + c5*wf(n,3)  ! a_n+a
        end do ! n

      endif

      end
