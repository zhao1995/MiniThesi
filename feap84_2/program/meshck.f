c$Id:$
      subroutine meshck(ie,nty,ix,nie,nen,nen1,numnp,numel,nummat,errs)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Revise formats for outputs of errors             12/01/2007
c       2. Move set of final boundary codes to psetid       21/07/2007
c          Eliminate unused arguments.
c       3. Remove 'ndf' from argument list                  13/12/2008
c       4. Add check that material set has been input       15/04/2009
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Check mesh data to ensure nodes/elements input

c      Inputs:
c         ie(nie,*)      - Material set assembly information
c         nty(*)         - Nodal type
c         ix(nen1,*)     - Element nodal connection lists
c         nie            - Dimension of ie array
c         nen            - Maximum number of nodes/element
c         nen1           - Dimension for ix array
c         numnp          - Number of nodes in mesh
c         numel          - Number of elemenst in mesh
c         nummat         - Number of material sets

c      Outputs:
c         errs           - Flag, true if errors detected
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cblend.h'
      include  'corset.h'
      include  'elflag.h'
      include  'iofile.h'
      include  'mdata.h'
      include  'pointer.h'
      include  'comblk.h'

      logical   errs
      integer   i,ii,iel,ma,n,nen,nen1,nie,numnp,numel,nummat
      integer   ie(nie,*),ix(nen1,*),nty(*)

      save

c     Perform mesh checks to ensure nodes/elements input

      errs = .false.
      do n = 1,numel
        ma = ix(nen1,n)
        if (ma.le.0 .or. ma.gt.nummat) then
          write(ilg,2000) ma,n
          write(iow,2000) ma,n
          if(ior.lt.0) write(*,2000) ma,n
          write(  *,2000) ma,n
          errs = .true.
        elseif(ie(nie-1,ma).eq.0) then
          write(ilg,2000) ma,n
          write(iow,2000) ma,n
          if(ior.lt.0) write(*,2000) ma,n
          write(  *,2000) ma,n
          errs = .true.
        else
          do i = 1,nen
            ii = ix(i,n)
            if(ii.gt.numnp .or. ii.lt.0) then
              write(ilg,2001) ii,n
              write(iow,2001) ii,n
              if(ior.lt.0) write(*,2001) ii,n
              write(  *,2001) ii,n
              errs = .true.
            elseif(ii.ne.0 .and. nty(ii).lt.0) then
              write(ilg,2002) ii,n
              write(iow,2002) ii,n
              if(ior.lt.0) write(*,2002) ii,n
              write(  *,2002) ii,n
              errs = .true.
            endif
          end do ! i
        endif
      end do ! n

c     If supernodes used then

      if(numbd.gt.0) then
        if(numsn.gt.0) then
          if(numsd.gt.0) then
            call mshcksn(mr(np(162)),mr(np(164)),numsd,numsn,numbd,errs)
          else
            write(ilg,2003)
            write(iow,2003)
            if(ior.lt.0) write(*,2003)
            errs = .true.
          endif
        else
          write(ilg,2004)
          write(iow,2004)
          if(ior.lt.0) write(*,2004)
          errs = .true.
        endif
      endif

c     Set first and last element numbers for each material type
c     N.B. Limit set by dimension in include file elflag.h

      do ma = 1,min(80,nummat)
        do n = 1,numel
          if(ix(nen1,n).eq.ma) then
            elstart(ma) = n
            go to 100
          endif
        end do ! n
100     do n = numel,1,-1
          if(ix(nen1,n).eq.ma) then
            ellast(ma) = n
            go to 200
          endif
        end do ! n
200     continue
      end do ! ma

c     Check for rotation values

      do ma = 1,3
        dal(ma) = 0
        ral(ma) = 0
      end do ! ma

c     Check for angle condtions

      do ma = 1,nummat
        iel = ie(nie-1,ma)
        if(anglefl) then
          if(iel.gt.0) then
            do i = 1,2
              if(dal(i).eq.ia(i,iel)) then

              elseif(dal(i).eq.0) then
                dal(i) = ia(i,iel)
              else
                write(iow,2005) iel,i,ia(i,iel),dal(i)
                write(ilg,2005) iel,i,ia(i,iel),dal(i)
              endif
            end do ! i
            do i = 1,2
              if(ral(i).eq.ir(i,iel)) then

              elseif(ral(i).eq.0) then
                ral(i) = ir(i,iel)
              else
                write(iow,2005) iel,i,ir(i,iel),ral(i)
                write(ilg,2005) iel,i,ir(i,iel),ral(i)
              endif
            end do ! i
          elseif(iel.lt.0) then
            do i = 1,2
              if(dal(i).eq.ea(i,-iel)) then

              elseif(dal(i).eq.0) then
                dal(i) = ea(i,-iel)
              else
                write(iow,2005) iel,i,ea(i,-iel),dal(i)
                write(ilg,2005) iel,i,ea(i,-iel),dal(i)
              endif
            end do ! i
            do i = 1,2
              if(ral(i).eq.er(i,-iel)) then

              elseif(ral(i).eq.0) then
                ral(i) = er(i,-iel)
              else
                write(iow,2005) iel,i,er(i,-iel),ral(i)
                write(ilg,2005) iel,i,er(i,-iel),ral(i)
              endif
            end do ! i
          endif
        endif
        if(eulerfl .or. triadfl) then
          if(iel.gt.0) then
            do i = 1,3
              if(dal(i).eq.ia3(i,iel)) then

              elseif(dal(i).eq.0) then
                dal(i) = ia3(i,iel)
              else
                write(iow,2005) i,iel,ia3(i,iel),dal(i)
                write(ilg,2005) i,iel,ia3(i,iel),dal(i)
              endif
            end do ! i
            do i = 1,3
              if(ral(i).eq.ir3(i,iel)) then

              elseif(ral(i).eq.0) then
                ral(i) = ir3(i,iel)
              else
                write(iow,2005) iel,i,ir3(i,iel),ral(i)
                write(ilg,2005) iel,i,ir3(i,iel),ral(i)
              endif
            end do ! i
          elseif(iel.lt.0) then
            do i = 1,3
              if(dal(i).eq.ea3(i,-iel)) then

              elseif(dal(i).eq.0) then
                dal(i) = ea3(i,-iel)
              else
                write(iow,2005) iel,i,ea3(i,-iel),dal(i)
                write(ilg,2005) iel,i,ea3(i,-iel),dal(i)
              endif
            end do ! i
            do i = 1,3
              if(ral(i).eq.er3(i,-iel)) then

              elseif(ral(i).eq.0) then
                ral(i) = er3(i,-iel)
              else
                write(iow,2005) iel,i,er3(i,-iel),ral(i)
                write(ilg,2005) iel,i,er3(i,-iel),ral(i)
              endif
            end do ! i
          endif
        endif
      end do ! ma

c     Formats

2000  format(5x,'*ERROR* MESHCK: Data for material set',i5,' on',
     &          ' element ',i10,' not input')

2001  format(5x,'*ERROR* MESHCK: Data for node ',i10,' on element'/
     &       5x,'         ',i10,' greater than maximum or negative')

2002  format(5x,'*ERROR* MESHCK: Data for node ',i10,' on element'/
     &       5x,'         ',i10,' not input')

2003  format(5x,'*ERROR* MESHCK: Blending functions used but no',
     &           ' SIDEs exist')

2004  format(10x,' *ERROR* MESHCK: Blending functions used but no',
     &           ' SNODes exist')

2005  format(10x,' *ERROR* MESHCK: Rotation of dof incompatible'/
     &       10x,'         ELEMENT TYPE =',i3,' DOF =',i3,' Current =',
     &           i3,' Previous =',i3)

      end
