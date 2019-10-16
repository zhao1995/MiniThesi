c$Id:$
      subroutine pouta(lct)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Output of arrays for use with Matlab sparse options.

c      Use:
c                Example:  Output of residual
c                     form
c                     output dr
c                Example:  Output of tangent
c                     tang,,-1
c                     output tang
c                Creates files with name 'dr' and 'tang'.
c                     Format: i  j  a(i,j)

c      Matlab use:
c                load dr
c                b = sparse(dr(:,1),dr(:,2),dr(:,3))
c                load tang
c                a = sparse(tang(:,1),tang(:,2),tang(:,3))

c      Inputs:
c         lct       - Command character parameters

c      Outputs:
c         To files with array name
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'compas.h'
      include   'iodata.h'
      include   'iofile.h'
      include   'part0.h'
      include   'pmatlab.h'
      include   'pointer.h'
      include   'comblk.h'

      logical    pcomp
      character  lct*15,array*15
      integer    i

      save

c     Get array name

      array      = ' '
      array(1:4) = lct(1:4)
      i          = index(array,' ')
      array(i:i) = '_'
      i          = i + 1

c     Tangent terms

      if(pcomp(array(1:4),'tang',4)) then
        call setcount(array,ntang,i)
        open(unit = ios,file = array,status = 'unknown')
        rewind(ios)
        if(ittyp.eq.-1 .or. ittyp.eq.-2) then  ! Blocked or Sparse
          if(max(abs(np(93)),abs(np(94)),abs(np(npart))).eq.0) then
            ntang = ntang - 1
            go to 400
          else
            call pstang(neq,mr(np(93)),mr(np(94)),hr(np(npart)))
          endif
        elseif(ittyp.eq.-3) then               ! Profile
          if(max(abs(np(20+npart)),abs(np(npart))).eq.0) then
            ntang = ntang - 1
            go to 400
          else
            call pptang(neq,mr(np(20+npart)),hr(np(npart)),
     &                  hr(np(npart)+neq))
          endif
        endif
      elseif(pcomp(array,'utan',4)) then
        if(ittyp.eq.-3) then               ! Profile
          if(max(abs(np(20+npart)),abs(np(npart)),
     &                             abs(np(npart+4))).eq.0) then
            go to 400
          else
            call setcount(array,nutan,i)
            open(unit = ios,file = array,status = 'unknown')
            rewind(ios)
            call pptang(neq,mr(np(20+npart)),hr(np(npart)),
     &                  hr(np(npart+4)))
          endif
        endif

c     Mass terms

      elseif(pcomp(array,'lmas',4)) then
        if(abs(np(npart+12)).eq.0) then
          go to 400
        else
          call setcount(array,nlmas,i)
          open(unit = ios,file = array,status = 'unknown')
          rewind(ios)
          call plmass(neq,hr(np(npart+12)))
        endif
      elseif(pcomp(array,'mass',4) .or. pcomp(array,'cmas',4)) then
        if(max(abs(np(90)),abs(np(91)),abs(np(npart+8))).eq.0) then
          go to 400
        else
          call setcount(array,ncmas,i)
          open(unit = ios,file = array,status = 'unknown')
          rewind(ios)
          call psmass(neq,mr(np(90)),mr(np(91)),hr(np(npart+8)),2)
        endif
      elseif(pcomp(array,'umas',4)) then
        if(max(abs(np(90)),abs(np(91)),abs(np(npart+8))).eq.0) then
          go to 400
        else
          call setcount(array,numas,i)
          open(unit = ios,file = array,status = 'unknown')
          rewind(ios)
          call psmass(neq,mr(np(90)),mr(np(91)),hr(np(npart+8)),3)
        endif

c     Damping terms

      elseif(pcomp(array,'damp',4) .or. pcomp(array,'cdam',4)) then
        if(max(abs(np(203)),abs(np(204)),abs(np(npart+16))).eq.0) then
          go to 400
        else
          call setcount(array,ncdam,i)
          open(unit = ios,file = array,status = 'unknown')
          rewind(ios)
          call psmass(neq,mr(np(203)),mr(np(204)),hr(np(npart+16)),2)
        endif
      elseif(pcomp(array,'udam',4)) then
        if(max(abs(np(203)),abs(np(204)),abs(np(npart+16))).eq.0) then
          go to 400
        else
          call setcount(array,nudam,i)
          open(unit = ios,file = array,status = 'unknown')
          rewind(ios)
          call psmass(neq,mr(np(203)),mr(np(204)),hr(np(npart+16)),3)
        endif

c     Residual terms

      elseif(pcomp(array,'dr  ',2) .or. pcomp(array,'form',4)) then
        if(abs(np(26)).eq.0) then
          go to 400
        else
          call setcount(array,nrfrm,i)
          open(unit = ios,file = array,status = 'unknown')
          rewind(ios)
          call prform(neq,hr(np(26)))
        endif
      endif
      close(unit = ios, status = 'keep')

      return

c     Error

400   write(iow,4000) array
      write(ilg,4000) array
      if(ior.lt.0) then
        write(*,4000) array
      endif
      close(unit = ios, status = 'delete')

c     format

4000  format(' *ERROR* Array ',a,' can not be output -- missing data'/)

      end
