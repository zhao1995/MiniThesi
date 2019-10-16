c$Id:$
      subroutine bm2trn(s,cs,sn,nst,ndf,itype)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Transform 2-d frame to global frame

c      Inputs:
c         Itype: 1  Transform matrix s(nst,nst)
c                2  Transform vector s(nst,1)

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   nst,ndf,itype, i,j,n
      real*8    cs,sn,t,tol, s(nst,*)

      save

      data      tol/ 1.d-12 /

      if(cs.gt.1.0d0-tol) return

      if(itype.eq.1) then

        do i = 1,nst,ndf
          j = i + 1
          do n = 1,nst
            t      = s(n,i)*cs - s(n,j)*sn
            s(n,j) = s(n,i)*sn + s(n,j)*cs
            s(n,i) = t
          end do ! n
        end do ! i
        do i = 1,nst,ndf
          j = i + 1
          do n = 1,nst
            t      = s(i,n)*cs - s(j,n)*sn
            s(j,n) = s(i,n)*sn + s(j,n)*cs
            s(i,n) = t
          end do ! n
        end do ! i

      else

        do i=1,nst,ndf
          j = i + 1
          t      =  s(i,1)*cs + s(j,1)*sn
          s(j,1) = -s(i,1)*sn + s(j,1)*cs
          s(i,1) =  t
        end do ! i

      endif

      end
