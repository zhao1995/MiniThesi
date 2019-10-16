c$Id:$
      subroutine plate2d(d,ul,xl,ix,s,p,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Set 'j' to 0 and then add to nh3 (for int8 use)  17/04/2011
c-----[--.----+----.----+----.-----------------------------------------]
c     Programed by: Ferdinando Auricchio

c     Triangular plate: 3 dofs per node (w, theta-x, theta-y)
c                       2 bubble modes for rotations
c                       2 shear parameters

c     Mixed approach for shear stiffness.
c        Step 1: Condensation of shear terms
c        Step 2: Condensation of bubble terms

c     Three integration points are used.

c     Arguments:
c        d(*)      - specified parameter array
c        ul(ndf,*) - local nodal solution values
c        xl(ndm,*) - local nodal coordinate values
c        ix(*)     - node numbers
c        s(nst,nst) - finite element array (stiffness, mass, geometric
c                                           stiffness)
c        p(nst)     - finite element array (residual, lumped mass)
c        ndf        - number of degree of freedoms at node ( > or = 3 )
c        ndm        - spatial dimension of element         ( > or = 2 )
c        nst        - size of finite element arrays        ( > or = 9 )
c        isw        - solution option
c                   = 1: Input values and store in d(*) array
c                   = 2: Check mesh coordinate and connection inputs
c                        for errors
c                   = 3: Compute element residual (p) and stiffness (s)
c                   = 4: Output element results
c                   = 5: Compute mass (p,s)/geometric stiffness array(s)
c                   = 6: Compute element residual (p)
c                   = 8: Compute nodal projections

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Input parameters set as follows:

c         ndm = 2 (x,y cartesian coordinates at nodes)
c         ndf = 3 (w,theta-x,theta-y, at nodes)
c         nen = 3 nodes (counterclockwise around element)
c                  or
c         nen = 4 nodes (counterclockwise around element)

c**********************************************************************

      implicit  none

      include  'cdata.h'
      include  'eldata.h'
      include  'hdata.h'
      include  'iofile.h'
      include  'mdata.h'
      include  'strnum.h'

      include  'comblk.h'

      integer   ndf,ndm,nst,isw, i,j, tdof

      integer   ix(*)
      real*8    d(*),ul(ndf,nen,*),xl(ndm,*),s(nst,*),p(ndf,*)

      save

c     Go to correct array processor

      if(isw.eq.0 .and. ior.lt.0) then
        write(*,*) '   Plate2d: 2-d Plate Linear Elastic (3 or 4 node)'

c     Input material properties

      elseif(isw.eq.1) then
        write(iow,2000)
        if(ior.lt.0) write(*,2000)
        call inmate(d,tdof,ndf*4,4)

c       Check dimensions and material for errors

        if(ndm.ne.2 .or. ndf.lt.3 .or. nen.lt.3) then
          write(  *,4000) ndm,ndf,nen
          write(ilg,4000) ndm,ndf,nen
          write(iow,4000) ndm,ndf,nen
          call plstop()
        endif

c       Set plot sequence for 3-node triangle

        pstyp = 2

c       Set rotation parameters: theta-x = 2; theta-y = 3

        ea(1,-iel) = 2
        ea(2,-iel) = 3

        istv = 20

c       Deactivate dof in element for dof > 3

        do i = 4,ndf
          ix(i) = 0
        end do ! i

c     History variable manipulation: None currently required

      elseif(isw.eq.12) then

        return
c     Initialize element strains for activation

      elseif(isw.eq.17) then

        j = 0
        do i = 1,nel
          hr(nh3+j  ) = ul(1,i,1)
          hr(nh3+j+1) = ul(2,i,1)
          hr(nh3+j+2) = ul(3,i,1)
          j           = j + 3
        end do ! i

c     Initialize element strains for deactivation

      elseif(isw.eq.18) then

        do i = 0,3*nel-1
          hr(nh3+i) = 0.0d0
        end do ! i

c     Critical time step computation

      elseif(isw.eq.21) then

        call tcritnd(d,xl,ul,ndm,ndf,nel)

c     External node determination

      elseif(isw.eq.26) then

        call pcorner2d()

      else

c       Compute element properties

        if(nel.eq.3) then
          call plate2t(d,ul,xl,ix,s,p,ndf,ndm,nst,isw)
        elseif(nel.eq.4) then
          call plate2q(d,ul,xl,ix,s,p,ndf,ndm,nst,isw)
        else
          write(iow,3000) n
          if(ior.lt.0) write(*,3000) n
          call plstop()
        endif
      endif

c     Formats for input-output

2000  format(/5x,'E l a s t i c   P l a t e   E l e m e n t'/)

3000  format(' *ERROR* Plate Element',i8,' has more than 4-nodes')

4000  format(' *ERROR* Problem control record incorrect:'/
     &    '         Mesh dim. (ndm) = ',i4,' - Should be 2'/
     &    '         DOFs/Node (ndf) = ',i4,' - Should be 3 or more'/
     &    '         Nodes/Elm (nen) = ',i4,' - Should be 3 or more')

      end
