c$Id:$
      subroutine m_mult (a,b,c,nra,nca,nrb,ncb,nrc,ncc,fltra,fltrb)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise            July 10, 1996            1.0

c      Acronym: Matrix MULTiplication

c      Purpose: Perform rows by columns matrix multiplication

c      Inputs:
c         a       - First matrix
c         b       - Second matrix
c         nra     - # of rows of A
c         nca     - # of columns of A
c         ncb     - # of columns of B

c      Outputs:
c         c       - Resulting matrix
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      logical   fltra, fltrb
      integer   nra,nca,nrb,ncb,nrc,ncc, kr,kc,kca
      real*8    a(nra,nca),b(nrb,ncb),c(nrc,ncc)

      save

      call pzero (c,nrc*ncc)

      if(.not. fltra .and. .not. fltrb) then
        if(nca.ne.nrb) then
          write(*,*) '   DIMENSIONS DO NOT MATCH'
        endif
        do kr = 1,nrc
          do kc = 1,ncc
            do kca = 1,nca
              c(kr,kc) = c(kr,kc) + a(kr,kca)*b(kca,kc)
            end do
          end do
        end do
      elseif(fltra .and. .not. fltrb) then
        if(nra.ne.nrb) then
          write(*,*) '   DIMENSIONS DO NOT MATCH'
        endif
        do kr = 1,nrc
          do kc = 1,ncc
            do kca = 1,nra
              c(kr,kc) = c(kr,kc) + a(kca,kr)*b(kca,kc)
            end do
          end do
        end do
      elseif(.not. fltra .and. fltrb) then
        if(nca.ne.ncb) then
          write(*,*) '   DIMENSIONS DO NOT MATCH'
        endif
        do kr = 1,nrc
          do kc = 1,ncc
            do kca = 1,nca
              c(kr,kc) = c(kr,kc) + a(kr,kca)*b(kc,kca)
            end do
          end do
        end do
      elseif(fltra .and. fltrb) then
        if(nra.ne.ncb) then
          write(*,*) '   DIMENSIONS DO NOT MATCH'
        endif
        do kr = 1,nrc
          do kc = 1,ncc
            do kca = 1,nra
              c(kr,kc) = c(kr,kc) + a(kca,kr)*b(kc,kca)
            end do
          end do
        end do
      endif

      end
