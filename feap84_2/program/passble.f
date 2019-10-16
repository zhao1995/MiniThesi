c$Id:$
      subroutine passble(s,p,ld,ix, jp,a,al,b,
     &                   alfl,aufl,bfl,dfl,gfl,rel, nsp,nov,jsw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Assemble global arrays from element arrays

c      Inputs:
c        s(*)   - Element array
c        p(*)   - Element vector
c        ld(*)  - Assembly numbers
c        ix(*)  - Global nodes for elements
c        jp(*)  - Pointer array
c        alfl   - Lower assembly if .true.
c        aufl   - Upper assembly if .true.
c        bfl    - Vector assembly if .true.
c        dfl    - Full vector if .true.
c        gfl    - Real assembly flag
c        rel    - Rigid element if .true.
c        nsp    - Size of 's' and 'p'
c        nov    - Number of elements connected to 'ix'
c        jsw    - Switch parameter

c      Outputs:
c        a(*)   - Diagonal and upper global array
c        al(*)  - Lower global array
c        b(*)   - Global vector
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'compac.h'
      include   'compas.h'
      include   'complx.h'
      include   'eldata.h'
      include   'eldatp.h'
      include   'eqsym.h'
      include   'prstrs.h'
      include   'sdata.h'
      include   'pointer.h'
      include   'comblk.h'

      include   'setups.h'

      logical    alfl,aufl,bfl,dfl,gfl,rel
      integer    jsw,nsp,nov, i,j, ld(*), ix(*),jp(*)
      real*8     p(*),s(nsp,*), a(*),al(*),b(*)

      save

c     Add element projections to total array

      if(jsw.eq.8) then
        if(aufl) then
          call dasble(s,p,ix,mr(np(197)),nel,numnp,aufl,bfl,
     &                b,al,a(numnp+1),a)
        else

c         Prevent accumulation of weights for multiple elements

          nov = nov + 1
          if(nov.gt.1) then
            do i = 1,nen
              p(i) = 0.0d0
            end do ! i
          endif
          call seproj(p,s,p(nen+1),hr(nph),hr(nph+numnp),hr(ner),ix,
     &                nel,nen,numnp)
          if(hpltfl) then
            call heproj(hr(np(304)),hr(np(305)),ix,nel,nen,numnp)
          endif
        endif

c     Stiffness and residual assembly

      elseif(aufl.or.bfl) then

c       Transform and assemble rigid body part

        if(np(167).ne.0) then
          call rasblk(s,p,ld,ix,mr(np(100)),mr(np(167)),hr(np(43)),
     &                aufl,bfl,jp,b,al,a(neq+1),a)
        elseif(rel) then
          call rasbly(s,p,hr(np(44)),hr(np(41)),ld,jp,ix,mr(np(100)),
     &                mr(np(96)),mr(np(99)),hr(np(95)),ndm,ndf,nel,nen,
     &                nsp,alfl,aufl,bfl,dfl, b,al,a(neq+1),a)
        else
          if(cplxfl) then

c           Assemble complex real part

            call dasble(s,p,ld,jp,nsp,neqs,aufl,bfl,b,al,a(2*neq+1),a)

c           Assemble complex imaginary part

            j = neq+1
            if(compfl) then
              i = jcmplx  + 1
              if(kdiag.or.ddiag.or.mdiag.or.udiag) then
                j = i
              endif
            else
              i = jp(neq) + 1
            endif
            call dasble(s(1,nsp+1),p(nsp+1),ld,jp,nsp,neqs,aufl,bfl,
     &                  b(neq+1),al(i),a(i+2*neq),a(j))
          elseif(gfl) then

c           Assemble for real arithmetic

            call dasble(s,p,ld,jp,nsp,neqs,aufl,bfl,b,al,a(neq+1),a)

          endif
        endif
      endif

      end
