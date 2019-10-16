c$Id:$
      subroutine modlbm(d,ta,eps,yy,zz,ww, hn,h1,nh,ii,istrt,
     &                  mhook, stress, isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    11/11/2008
c-----[--+---------+---------+---------+---------+---------+---------+-]
c      Purpose:  Compute local beam stress components from 3-d models

c      Inputs:
c        d(*)       - Material parameters
c        ta         - Temperature change
c        eps(3)     - Beam strains on cross section
c                     eps(1) = xx; eps(2) = xy; eps(3) = xz
c        hn(*)      - History variables at t_n
c        h1(*)      - History variables at t_n+1
c        nh         - Number history terms
c        ii         - Computation location
c        istrt      - Iteration start flag
c        isw        - FEAP element switch value

c      Outputs:
c        mhook(6,6) - Beam resultant moduli
c        stress(6)  - Beam resultant stresses
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      integer    nh,ii,istrt, isw
      real*8     ta, yy,zz,ww
      real*8     d(*),eps(3),hn(*),h1(*), mhook(6,6), stress(6)

c     Local variables

      logical    notconv
      integer    i,j, nit,nitmax, imap(6)
      real*8     dsnrm, snrm, tol, alam,ha, df
      real*8     strs(6), strn(6,2), ds(6,6,5), dm(6,6), tau(6)
      real*8     dstrn(6), tmp(3,3), dd(3,3), sig(3)

      real*8     dot

      data       nitmax / 10 /
      data       imap   / 1, 4, 6, 2, 3, 5 /
      data       tol    / 1.d-10 /

c     Map strains to initial 3-d system

      strn(1,1) = eps(1)
      strn(2,1) = hn(1)
      strn(3,1) = hn(2)
      strn(4,1) = eps(2)
      strn(5,1) = hn(3)
      strn(6,1) = eps(3)

      alam      = 0.0d0
      ha        = 0.0d0

      nit   = 0
      dsnrm = 1.0d0
      snrm  = 1.0d0

c     Set not converged to start iterations

      notconv  = .true.

      do while ((nit.lt.nitmax) .and. notconv)
        nit = nit + 1
        call modlsd(ii,d,ta,strn,hn(4),h1(4),nh-3,istrt,ds,strs,
     &              alam,ha,isw)

        do i = 1,6
          tau(i) = strs(imap(i))
          do j = 1,6
            dm(j,i) = ds(imap(j),imap(i),1)
          end do ! j
        end do ! i

c       Compute stress

        call invert(dm(4,4),3,6)

        do i = 4,6
          dstrn(i) = - dm(i,4)*strs(2)
     &               - dm(i,5)*strs(3)
     &               - dm(i,6)*strs(5)
        end do ! i

        strn(2,1) = strn(2,1) + dstrn(4)
        strn(3,1) = strn(3,1) + dstrn(5)
        strn(5,1) = strn(5,1) + dstrn(6)
        dsnrm = dstrn(4)**2  + dstrn(5)**2  + dstrn(6)**2
        snrm  = dot(strn,strn,6)

c       Check convergence

        notconv  = dsnrm.gt.tol*snrm

      end do ! while

c     Warning on lack of convergence

      if(notconv) then
        write( *,*) ' *WARNING* NIT-tot:',nit,dsnrm,snrm
      endif

c     Save strains in history

      h1(1) = strn(2,1)
      h1(2) = strn(3,1)
      h1(3) = strn(5,1)

c     Remap moduli and stresses

      do i = 1,3
        do j = 1,3
          tmp(j,i) = dm(j,4)*dm(4,i+3)
     &             + dm(j,5)*dm(5,i+3)
     &             + dm(j,6)*dm(6,i+3)
        end do ! j
      end do ! i
      do i = 1,3
        do j = 1,3
          dm(j,i) = dm(j,i) - tmp(j,1)*dm(4,i)
     &                      - tmp(j,2)*dm(5,i)
     &                      - tmp(j,3)*dm(6,i)
        end do ! j
      end do ! i

      do i = 1,3
        sig(i) = tau(i)
        do j = 1,3
          dd(j,i) = dm(j,i)
        end do ! j
      end do ! i

c     Compute multi-dimensioinal contributions to beam section

      df        = sig(1)*ww

      stress(1) = stress(1) + sig(2)*ww
      stress(2) = stress(2) + sig(3)*ww
      stress(3) = stress(3) + df
      stress(4) = stress(4) + df*zz
      stress(5) = stress(5) - df*yy
      stress(6) = stress(6) + (sig(3)*yy - sig(2)*zz)*ww

      df         = dd(1,1)*ww

      mhook(1,1) = mhook(1,1) +  dd(2,2)*ww
      mhook(1,2) = mhook(1,2) +  dd(2,3)*ww
      mhook(1,3) = mhook(1,3) +  dd(2,1)*ww
      mhook(1,4) = mhook(1,4) +  dd(2,1)*zz*ww
      mhook(1,5) = mhook(1,5) -  dd(2,1)*yy*ww
      mhook(1,6) = mhook(1,6) + (dd(2,3)*yy - dd(2,2)*zz)*ww

      mhook(2,1) = mhook(2,1) +  dd(3,2)*ww
      mhook(2,2) = mhook(2,2) +  dd(3,3)*ww
      mhook(2,3) = mhook(2,3) +  dd(3,1)*ww
      mhook(2,4) = mhook(2,4) +  dd(3,1)*zz*ww
      mhook(2,5) = mhook(2,5) -  dd(3,1)*yy*ww
      mhook(2,6) = mhook(2,6) + (dd(3,3)*yy - dd(3,2)*zz)*ww

      mhook(3,1) = mhook(3,1) +  dd(1,2)*ww
      mhook(3,2) = mhook(3,2) +  dd(1,3)*ww
      mhook(3,3) = mhook(3,3) +  df
      mhook(3,4) = mhook(3,4) +  df*zz
      mhook(3,5) = mhook(3,5) -  df*yy
      mhook(3,6) = mhook(2,6) + (dd(1,3)*yy - dd(1,2)*zz)*ww

      mhook(4,1) = mhook(4,1) +  dd(1,2)*zz*ww
      mhook(4,2) = mhook(4,2) +  dd(1,3)*zz*ww
      mhook(4,3) = mhook(4,3) +  df*zz
      mhook(4,4) = mhook(4,4) +  df*zz*zz
      mhook(4,5) = mhook(4,5) -  df*zz*yy
      mhook(4,6) = mhook(4,6) + (dd(1,3)*yy-dd(1,2)*zz)*zz*ww

      mhook(5,1) = mhook(5,1) -  dd(1,2)*yy*ww
      mhook(5,2) = mhook(5,2) -  dd(1,3)*yy*ww
      mhook(5,3) = mhook(5,3) -  df*yy
      mhook(5,4) = mhook(5,4) -  df*zz*yy
      mhook(5,5) = mhook(5,5) +  df*yy*yy
      mhook(5,6) = mhook(5,6) - (dd(1,3)*yy-dd(1,2)*zz)*yy*ww

      mhook(6,1) = mhook(6,1) + (dd(3,2)*yy - dd(2,2)*zz)*ww
      mhook(6,2) = mhook(6,2) + (dd(3,3)*yy - dd(2,3)*zz)*ww
      mhook(6,3) = mhook(6,2) + (dd(3,1)*yy - dd(2,1)*zz)*ww
      mhook(6,4) = mhook(6,4) + (dd(3,1)*yy-dd(2,1)*zz)*zz*ww
      mhook(6,5) = mhook(6,5) - (dd(3,1)*yy-dd(2,1)*zz)*yy*ww
      mhook(6,6) = mhook(6,6) + (dd(2,2)*zz*zz + dd(3,3)*yy*yy
     &                        - (dd(2,3) + dd(3,2))*yy*zz)*ww

      end
