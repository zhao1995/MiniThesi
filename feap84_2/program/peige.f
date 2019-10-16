c$Id:$
      subroutine peige(s,nst, dr,vflg,prtev)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add print flage 'prtev'                          31/10/2007
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Compute the eigenvalues and vectors for last
c                computed element array (numel array)

c      Inputs:
c         s(nst,*) - Last element array
c         nst      - Dimension of element array
c         vflg     - Flag, compute vectors if true
c         prtev    - Output vectors if true

c      Outputs:
c         dr(*)    - Array of element eigenvalues and vectors
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'pointer.h'
      include  'comblk.h'

      include  'p_int.h'

      logical   vflg,prtev
      integer   nst, i, j, j1, n1,n2,n3
      real*8    s(nst,nst), dr(*)

      save

c     Set pointers

      n1 =  1 +  nst
      n2 = n1 + (nst*(nst+1))/2
      n3 = n2 +  nst*nst

c     Load stiffness terms into triangular matrix

      j1 = -1
      do j = 1,nst
        do i = 1,j
          dr(n1+j1+i) = s(i,j)
        end do ! i
        j1 = j1 + j
      end do ! j

c     Compute eigenpairs for last element computed

      call eisql(dr(n1),dr(1),dr(n3),dr(n2),nst,j1)

c     Move eigenvalues and vectors to 'EIGE' storage for plots.

      do i = 0,nst*nst-1
        hr(np(75)+i) = dr(n2+i)
      end do ! i
      fp(1) = np(75) + nst*nst - 1
      do i = 1,nst
        hr(fp(1)+i) = dr(i)
      end do ! i

c     Output eigenpairs for element

      if(prtev) then
        call mprint (  dr,       1, nst,  1,'Eigenvalue')
        if(vflg) then
          call mprint (dr(n2), nst, nst,nst,'Eigenvect.')
        endif
      endif

      end
