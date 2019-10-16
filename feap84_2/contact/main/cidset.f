c$Id:$
      subroutine cidset(cs0, ics, ip, ndm,ndf)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    25/01/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Robert L. Taylor         January 25, 2013            1.0

c      Acronym: Contact ID SET

c      Purpose: Set flags to be sure contact nodes are active

c      Inputs :
c         cs0(*)  - Surface table
c         ics(*)  - Surface node numbers (packed)
c         ndm     - Mesh dimension
c         ndf     - No. DOF checked

c      Outputs:
c         ip(*)   - Set of nodes (in first entry)
c-----[--.----+----.----+----.-----------------------------------------]

      implicit   none

      include   'c_0.h'
      include   'c_contac.h'
      include   'c_comnd.h'

      integer    ndm,ndf
      integer    ip(ndf+1,*), ics(*)
      real*8     cs0(nr0,n0c1:nc01,*)

      integer    ofs,neps,dnope,nope, ke,kn,ks,k, n0,nod

c     Perform set of nodes on contact surfaces

      do ks = 1,numcs

c       Find dimensioning information

        ofs  = nint(cs0(2,-1,ks))
        neps = nint(cs0(3,-1,ks))
        dnope= nint(cs0(4,-1,ks))
        nope = nint(cs0(2,0,ks))

c       Loop over nodes

        do ke = 1,neps
          n0 = ofs + (ke-1)*dnope - 1
          do kn = 1,nope
            nod = ics(n0+kn)
            if(nod.gt.0) then
              ip(1,nod) = 1
              if(ndf.gt.0) then
                do k = 1,ndm
                  ip(k+1,nod) = 1
                end do ! k
              endif ! ndf > 0
            endif
          end do ! kn
        end do ! ke
      end do ! ks

      end
