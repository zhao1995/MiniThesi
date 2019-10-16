c$Id:$
      subroutine ploadl(id,jp,ld,b,ad,al,au,p,s,x,xl,u,ul,
     1                  prop,ndf,ndm,dfl,aufl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Assemble surface loading from element contributions.
c               Data is input using ploadi routine.

c      Inputs:
c         id(ndf,*)  - Active equation numbers
c         jp(*)      - Pointers for row/columns in tangent array
c         x(ndm,*)   - Nodal coordinates of mesh
c         u(ndf,*)   - Nodal solutions of mesh
c         prop       - Current total proportional load level
c         ndf        - Number dof/node
c         ndm        - Spatial dimension of mesh
c         dfl        - Flag, true if uncompressed vector assembled
c         aufl       - Flag, form tangent if true

c      Scratch:
c         ld(ndf,*)      - Element local/global equation numbers
c         p(*)       - Element vector
c         s(*)       - Element matrix
c         xl(ndm,*)  - Element nodal coordinates
c         ul(ndf,*)  - Element solution/rate parameters

c      Outputs:
c         b(*)       - Includes surface load effect in residual/reaction
c         ad(*)      - Diagonal part of tangent array from surface load
c         al(*)      - Lower part of tangent array
c         au(*)      - Lower part of tangent array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'cornum.h'
      include  'corset.h'
      include  'eldata.h'
      include  'eqsym.h'
      include  'mdata.h'
      include  'part0.h'
      include  'pointer.h'
      include  'sldata.h'
      include  'comblk.h'

      include  'p_int.h'

      logical   dfl,aufl
      integer   ndf,ndm, i,j, ii,jj, nn,nrec,nrot(3),ns1, ns,nv
      integer   id(ndf,*),jp(*),ld(ndf,*)
      real*8    prop, b(*),p(*),xl(ndm,*),x(ndm,*),u(ndf,*),ul(ndf,*)
      real*8    valu(64),s(*),ad(*),al(*),au(*)

      save

c     Compute surface loads for current deformation

      if(numsl.gt.0) then

        jj = 3
        if(aufl) jj = 2
        do nn = 1,numsl
          iel = iels(1,nn)
          nim = iels(2,nn)
          nre = iels(3,nn)
          ma  = iels(4,nn)
          nrec= iels(5,nn)
          ns  = inods(nn)
          nv  = ivals(nn)
          nel = ns            ! variable passed to element in /eldata/
          mct = nv            ! variable passed to element in /eldata/
          ns1 = ndf*nel
          do j = 1,nrec
            do n = 1,nel
              ii = mr(np(37)+nim+n-1)
              do i = 1,ndm
                xl(i,n) = x(i,ii)
              end do ! i
              do i = 1,ndf
                ul(i,n) = u(i,ii)
                ld(i,n) = id(i,ii)
                if(ndfp(i).ne.npart) ld(i,n) = 0
              end do ! i
            end do ! n
            if(anglefl) then
              call pangl(mr(np(37)+nim),nel,hr(np(46)),hr(np(45)),
     &                   nrot(1))
            else
              nrot(1) = 0
            endif
            if(eulerfl) then
              call peule(mr(np(37)+nim),nel,hr(np(243)),hr(np(242)),
     &                   nrot(2))
            else
              nrot(2) = 0
            endif
            if(triadfl) then
              call pltriad(mr(np(37)+nim),nel,hr(np(275)),hr(np(274)),
     &                   nrot(3))
            else
              nrot(3) = 0
            endif
            call pzero(p,ns1)
            call pzero(s,ns1*ns1)
            call pzero(valu,mct)
            call colred(hr(np(241)+nre),-prop,mct, valu)

c           Rotate displacements by 'angl'

            if(nrot(1).gt.0) then
              call ptrans(dal,hr(np(46)),ul,p,s,nel,ndf,ns1,1)
              if(ral(1).ne.0) then
                call ptrans(ral,hr(np(46)),ul,p,s,nel,ndf,ns1,1)
              endif
            endif

c           Rotate displacements by 'euler'

            if(nrot(2).gt.0) then
              call petrans(dal,hr(np(243)),ul,p,s,nel,ndf,ns1,1)
              if(ral(1).ne.0) then
                call petrans(ral,hr(np(243)),ul,p,s,nel,ndf,ns1,1)
              endif
            endif

c           Rotate displacements by 'triad'

            if(nrot(3).gt.0) then
              call pttrans(dal,hr(np(275)),ul,p,s,nel,ndf,ns1,1)
              if(ral(1).ne.0) then
                call pttrans(ral,hr(np(275)),ul,p,s,nel,ndf,ns1,1)
              endif
            endif

c           Form element contributions

            call elmlib(valu,ul,xl,mr(np(37)+nim),valu,s,p,
     &                  ndf,ndm,ns1,iel,7)

c           Rotate element arrays by 'angl'

            if(nrot(1).gt.0) then
              call ptrans(dal,hr(np(46)),ul,p,s,nel,ndf,ns1,jj)
              if(ral(1).ne.0) then
                call ptrans(ral,hr(np(46)),ul,p,s,nel,ndf,ns1,jj)
              endif
            endif

c           Rotate element arrays by 'euler'

            if(nrot(2).gt.0) then
              call petrans(dal,hr(np(243)),ul,p,s,nel,ndf,ns1,jj)
              if(ral(1).ne.0) then
                call petrans(ral,hr(np(243)),ul,p,s,nel,ndf,ns1,jj)
              endif
            endif

c           Rotate displacements by 'triad'

            if(nrot(3).gt.0) then
              call pttrans(dal,hr(np(275)),ul,p,s,nel,ndf,ns1,jj)
              if(ral(1).ne.0) then
                call pttrans(ral,hr(np(275)),ul,p,s,nel,ndf,ns1,jj)
              endif
            endif

            call dasble(s,p,ld,jp,ns1,neqs,aufl,.not.dfl, b,al,au,ad)

c           Assemble an uncompressed vector if needed

            if(dfl) then
              do n = 1,nel
                ii = mr(np(37)+nim+n-1)
                do i = 1,ndf
                  if(id(i,ii).gt.0. and. ndfp(i).eq.npart ) then
                    ld(i,n) = ii*ndf - ndf + i
                  else
                    ld(i,n) = 0
                  endif
                end do
              end do
              call dasble(s,p,ld,jp,ns1,neqs,.false.,dfl, b,al,au,ad)
            endif
            nim = nim + ns
            nre = nre + nv
          end do ! j
        end do ! nn

      end if

c     Compute follower load contributions

      if(nfol.gt.0) then
        nel = 2
        ns1 = ndf*nel
        jj  = 3
        if(aufl) jj = 2
        do nn = 0,nfol-1
          fp(1) = np(129) + nn*2
          fp(2) = np(130) + nn
          do j = 1,2
            ii  = mr(fp(1) + j - 1)
            do i = 1,ndm
              xl(i,j) = x(i,ii)
            end do ! i
            do i = 1,ndf
              ul(i,j) = u(i,ii)
              ld(i,j) = id(i,ii)
              if(ndfp(i).ne.npart) ld(i,j) = 0
            end do ! i
          end do ! j
          if(anglefl) then
            call pangl(mr(fp(1)),nel,hr(np(46)),hr(np(45)),nrot(1))
          else
            nrot(1) = 0
          endif
          if(eulerfl) then
            call peule(mr(np(37)+nim),nel,hr(np(243)),hr(np(242)),
     &                 nrot(2))
          else
            nrot(2) = 0
          endif
          if(triadfl) then
            call pltriad(mr(np(37)+nim),nel,hr(np(275)),hr(np(274)),
     &                 nrot(3))
          else
            nrot(3) = 0
          endif

c         Rotate displacements by 'angl'

          if(nrot(1).gt.0) then
            call ptrans(dal,hr(np(46)),ul,p,s,nel,ndf,ns1,1)
            if(ral(1).ne.0) then
              call ptrans(ral,hr(np(46)),ul,p,s,nel,ndf,ns1,1)
            endif
          endif

c         Rotate displacements by 'euler'

          if(nrot(2).gt.0) then
            call petrans(dal,hr(np(243)),ul,p,s,nel,ndf,ns1,1)
            if(ral(1).ne.0) then
              call petrans(ral,hr(np(243)),ul,p,s,nel,ndf,ns1,1)
            endif
          endif

c         Rotate displacements by 'triad'

          if(nrot(3).gt.0) then
            call pttrans(dal,hr(np(275)),ul,p,s,nel,ndf,ns1,1)
            if(ral(1).ne.0) then
              call pttrans(ral,hr(np(275)),ul,p,s,nel,ndf,ns1,1)
            endif
          endif

c         Form element contributions

          call pfolel(xl,ul,hr(fp(2)),p,s,prop,ndm,ndf,ns1)

c         Rotate element arrays by 'angl'

          if(nrot(1).gt.0) then
            call ptrans(dal,hr(np(46)),ul,p,s,nel,ndf,ns1,jj)
            if(ral(1).ne.0) then
              call ptrans(ral,hr(np(46)),ul,p,s,nel,ndf,ns1,jj)
            endif
          endif

c         Rotate element arrays by 'euler'

          if(nrot(2).gt.0) then
            call petrans(dal,hr(np(243)),ul,p,s,nel,ndf,ns1,jj)
            if(ral(1).ne.0) then
              call petrans(ral,hr(np(243)),ul,p,s,nel,ndf,ns1,jj)
            endif
          endif

c         Rotate displacements by 'triad'

          if(nrot(3).gt.0) then
            call pttrans(dal,hr(np(275)),ul,p,s,nel,ndf,ns1,jj)
            if(ral(1).ne.0) then
              call pttrans(ral,hr(np(275)),ul,p,s,nel,ndf,ns1,jj)
            endif
          endif

          call dasble(s,p,ld,jp,ns1,neqs,aufl,.not.dfl, b,al,au,ad)

c         Assemble an uncompressed vector if needed

          if(dfl) then
            do j = 1,nel
              ii = mr(fp(1) + j - 1)
              do i = 1,ndf
                if(id(i,ii).gt.0. and. ndfp(i).eq.npart ) then
                  ld(i,j) = ii*ndf - ndf + i
                else
                  ld(i,j) = 0
                endif
              end do ! i
            end do ! j
            call dasble(s,p,ld,jp,ns1,neqs,.false.,dfl, b,al,au,ad)
          endif
        end do ! nn
      endif

      end
