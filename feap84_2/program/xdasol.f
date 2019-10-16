c$Id:$
      subroutine xdasol(fau, fal, iunau, iunal, b, jp, neq,
     &                  energy, ad, au, al, iblk, maxbl, scale)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Solution of symmetric equations stored in profile
c               form.  Coefficient matrix must be decomposed into
c               its triangular factors using xdatri before using
c               xdasol.

c      Inputs:
c         fau       - Name of file containing upper blocks
c         fal       - Name of file containing lower blocks
c         iunau     - Logical unit number for upper blocks
c         iunal     - Logical unit number for lower blocks
c         b(*)      - Right hand side of equations
c         jp(*)     - Pointer to end of row/columns in profile
c         neq       - Number of equations
c         ad(*)     - Diagonals of factored matrix
c         au(*)     - Upper profile of factored matrix
c         al(*)     - Lower profile of factored matrix
c         iblk(3,*) - Block information
c                     iblk(1,I) = first equation in block I
c                     iblk(2,I) = last  equation in block I
c                     iblk(3,I) = Number of block needed to start
c         maxblk    - Number of blocks to store profile
c         scale     - Scale equations

c      Outputs:
c         b(*)      - Solution vector
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'fdata.h'
      include  'iofile.h'
      include  'part0.h'
      include  'pointer.h'
      include  'comblk.h'

      logical   scale
      character fau*(*), fal*(*)
      integer   i,j,ij,is,ie,ista,jr,jh,iunau,iunal,neq,isize,maxbl
      integer   jp(neq), iblk(3, maxbl)
      real*4    etime, tt,tary(2)
      real*8    dot, energy,bd, al(neq),au(neq),ad(neq),b(neq)

      save

c     Find first non-zero entry in right hand side

      do ie = 1,neq
        if(b(ie).ne.0.0d0) go to 100
      end do ! ie

c     Equations have zero right hand side: Solution zero!

      if(ior.gt.0) write(iow,2000)
      if(ior.lt.0) write(*,2000)
      return
100   if(ie.lt.neq) then

c       Scale equations

        if(scale) then
          call pscalb(b,hr(np(234+npart)),ie,neq)
        endif

c       Find first block for forward solution

        do i = 1, maxbl
          if(iblk(2, i).ge. ie) go to 110
        end do ! i

c       Reduce right hand side

110     j = ie + 1
        do ij = i, maxbl
          isize = jp(iblk(2,ij)) - jp(iblk(1,ij)-1)
          call rfile(iunal, fal, ij, isize, al)
          ista  = iblk(1, ij) - 1
          is =  j - 1
          do j = is+1,iblk(2, ij)
            if(jp(j).gt.jp(j-1)) then
              jh = jp(j)   - jp(j-1)
              jr = jp(j-1) - jp(ista) + 1
              b(j) = b(j) - dot(al(jr),b(j-jh),jh)
            endif
          end do ! j
          if(pfr .and. ior.lt.0) then
            tt = etime(tary)
            write(*,2001) ij, iblk(2,ij), tary
          endif
        end do ! ij
      endif

c     Multiply by inverse of diagonal elements and compute energy

      energy = 0.0d0
      do j = 1,neq
        bd     = b(j)
        b(j)   = b(j)*ad(j)
        energy = energy + bd*b(j)
      end do ! j

c     Backsubstitution

      if(neq.gt.1) then
        do ij = maxbl, 1, -1
          isize = jp(iblk(2,ij)) - jp(iblk(1,ij)-1)
          call rfile(iunau, fau, ij, isize, au)
          ista  = iblk(1, ij) - 1
          do j = iblk(2, ij), ista+1, -1
            if(jp(j).gt.jp(j-1)) then
              jh = jp(j)   - jp(j-1)
              jr = jp(j-1) - jp(ista) + 1
              call colred(au(jr),b(j),jh, b(j-jh))
            endif
          end do ! j
          if(pfr .and. ior.lt.0) then
            tt = etime(tary)
            write(*,2001) ij, iblk(2,ij), tary
          endif
        end do ! ij

c       Scale equations

        if(scale) then
          call pscalb(b,hr(np(234+npart)),1,neq)
        endif

      endif

c     Format

2000  format(' *WARNING* Zero right-hand-side vector')

2001  format('    --> Block',i4,' completed: No. Equations =',i8,
     &       '    t=',2f9.2)

      end
