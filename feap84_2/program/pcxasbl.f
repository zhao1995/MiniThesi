c$Id:$
      subroutine pcxasbl(omega,eta,id,f,d,adr,adi,aur,aui,jp,
     &                   k,m,c,jc,ir,nneq)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Assemble complex array in profilel form

c     Inputs:
c        omega     - Frequency
c        eta       - Damping ratio
c        id(ndf,*) - Equation numbers of nodes
c        f(ndf,*)  - Nodal force
c        d(ndf,*)  - Nodal displacement
c        jp(*)     - Column pointer array
c        k(*)      - Stiffness array
c        m(*)      - Mass      array
c        c(*)      - Damping   array
c        jc(*)     - Pointer for sparse storage of columns
c        ir(*)     - Row numbers for each column
c        nneq      - Total number of equations

c     Outputs:
c        adr(*)    - Real      diagonal part of array
c        adi(*)    - Imaginary diagonal part of array
c        aur(*)    - Real      diagonal part of array
c        aui(*)    - Imaginary diagonal part of array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'

      integer   nneq,n,nn, j,jj,jn
      real*8    omega,omega2, eta

      integer   id(*),jp(*),jc(*),ir(*)
      real*8    f(nneq,2),d(nneq), adr(*),adi(*), aur(*),aui(*)
      real*8    k(*),m(*),c(*)

      save

      omega2 = omega*omega

      do nn = 1,nneq
        n = id(nn)
        if(n.gt.0) then
          adr(n) = k(nn) - omega2*m(nn)
          adi(n) = omega*c(nn) + eta*k(nn)
        endif
      end do ! nn

      do nn = 2,nneq
        n  = id(nn)
        if(n.gt.0) then
          jn = jp(n) - n + 1
          do jj = jc(nn-1)+1,jc(nn)
            j = id(ir(jj))
            if(j.gt.0) then
              aur(j+jn) = k(jj+nneq) - omega2*m(jj+nneq)
              aui(j+jn) = omega*c(jj+nneq)
            else
              f(n,1) = f(n,1) - (k(jj+nneq)
     &                        - omega2*m(jj+nneq))*d(ir(jj))
              f(n,2) = f(n,2) - omega *c(jj+nneq) *d(ir(jj))
            endif
          end do ! jj
        else
          do jj = jc(nn-1)+1,jc(nn)
            j = id(ir(jj))
            if(j.gt.0) then
              f(j,1) = f(j,1) - (k(jj+nneq) - omega2*m(jj+nneq))*d(nn)
              f(j,2) = f(j,2) - omega*c(jj+nneq)*d(nn)
            endif
          end do ! jj
        endif
      end do ! nn

      end
