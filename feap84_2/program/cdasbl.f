c$Id:$
      subroutine cdasbl(s,p,ld,jp,ns,neqs,afl,bfl,
     +                  br,bi,alr,ali,aur,aui,adr,adi)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:    Assemble complex arrays for 'CDASOL'
c                  Global arrays stored as real and imaginary parts
c                  Local  arrays stored as double precision complex.
c                  First 'neqs' equations are symmetric, rest are
c                  unsymmetric.

c      Inputs:
c         s(*,*) - Complex element matrix
c         p(*)   - Complex element vector
c         ld(*)  - Global/local equation numbers for assembly
c         ns     - Size of 's', 'p', and 'ld' arrays
c         neqs   - Number of symmetric equations
c         afl    - Assemble 'A' if true
c         bfl    - Assemble 'b' if true

c      Outputs:
c         br(*)  - Real      part of rhs   entries
c         bi(*)  - Imaginary part of rhs   entries
c         alr(*) - Real      part of lower entries
c         ali(*) - Imaginary part of lower entries
c         aur(*) - Real      part of upper entries
c         aui(*) - Imaginary part of upper entries
c         adr(*) - Real      part of diagonal entries
c         adi(*) - Imaginary part of diagonal entries
c-----[--.----+----.----+----.-----------------------------------------]
      implicit    none

      logical     afl, bfl
      integer     i, ii, j, jj, je, ns, neqs
      integer     ld(ns),jp(*)
      real*8      alr(*),aur(*),adr(*),br(*)
      real*8      ali(*),aui(*),adi(*),bi(*)
      complex*16  s(ns,ns),p(ns)

      save

c     Loop through rows

      je = jp(neqs)

      if(afl) then
        do i = 1,ns
          if(ld(i).gt.0) then
            ii = ld(i) + 1

c           Loop through columns

            do j = 1,ns

              if(ld(j).eq.ld(i)) then

                adr(ld(i)) = adr(ld(i)) +  dble(s(i,j))
                adi(ld(i)) = adi(ld(i)) + dimag(s(i,j))

              elseif(ld(j).gt.ld(i)) then
                jj = ii + jp(ld(j)) - ld(j)

                aur(jj) = aur(jj) +  dble(s(i,j))
                aui(jj) = aui(jj) + dimag(s(i,j))

                if(ld(j).gt.neqs) then
                  alr(jj-je) = alr(jj-je) +  dble(s(j,i))
                  ali(jj-je) = ali(jj-je) + dimag(s(j,i))
                endif

              endif
            end do ! j
          endif
        end do ! i
      endif

c     Assemble a vector

      if(bfl) then
        do i = 1,ns
          if(ld(i).gt.0) then
            br(ld(i))  = br(ld(i))  +  dble(p(i))
            bi(ld(i))  = bi(ld(i))  + dimag(p(i))
          end if
        end do ! i
      endif

      end
