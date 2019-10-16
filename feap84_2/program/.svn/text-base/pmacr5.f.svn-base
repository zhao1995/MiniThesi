c$Id:$
      subroutine pmacr5 (lct,ct,j)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Change 'Macro' to 'Command'                      07/01/2008
c       2. Add 'set disp' command                           25/07/2008
c       3. Change 'avg' iter/step to real                   18/10/2008
c       4. Revise format 2008                               02/12/2008
c       5. Set number of total iterations for nstep = 0     04/01/2009
c       6. Output 'show' to file in interactive mode        21/03/2013
c       7. Modify to print "U s e r" for element types      07/08/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Command language instruction subprogram: Part 5

c      Inputs:
c         lct       - Command option
c         ct(3)     - Command parameters
c         j         - Command number in this routine

c      Outputs:
c         Depends on command number j
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include  'augdat.h'
      include  'bdata.h'
      include  'cdata.h'
      include  'comfil.h'
      include  'counts.h'
      include  'ddata.h'
      include  'elcount.h'
      include  'endata.h'
      include  'fdata.h'
      include  'iodata.h'
      include  'iofile.h'
      include  'mdata.h'
      include  'mxsiz.h'
      include  'ndata.h'
      include  'part0.h'
      include  'pdata6.h'
      include  'pglob1.h'
      include  'plflag.h'
      include  'pointer.h'
      include  'prlod.h'
      include  'rdata.h'
      include  'rigid1.h'
      include  'rjoint.h'
      include  'sdata.h'
      include  'tdata.h'
      include  'xtout.h'
      include  'comblk.h'

      include  'p_int.h'

      logical   pcomp, setvar,palloc, aflag
      character lct*15
      integer   i, j, nn,nnp,npr, nlg, ttot
      real*8    ct(3),avg

      save

c     [outm]      Output renumbered input file    - text mode
c     [outm,cont] Output contact slide-line data  - text mode
c     [outm,defo] Output renumbered deformed mesh - text mode
c     [outm,doma] Output mesh for domains         - text mode
c     [outm,elem] Output mesh for single element  - text mode
c     [outm,bina] Output renumbered input file    - binary mode
c     [renu]      Output list of old/renumbered nodes and coords.

      if(j.eq.1 .or. j.eq.2) then

        setvar = palloc(120,'TEMP0',numnp*3+3,1)
        call poutm (mr(np(34)),mr(np(32)),hr(np(25)),
     &              mr(np(31)+ndf*numnp),hr(np(43)),hr(np(40)),
     &              mr(np(33)),hr(np(45)),hr(np(27)),hr(np(38)),
     &              mr(np(29)),hr(np(86)),hr(np(87)),hr(np(88)),
     &              mr(np(120)),mr(np(120)+numnp+1),lct,ct,j)
        setvar = palloc(120,'TEMP0',0,1)

c     [show]      show current solution parameters
c     [show,cont] show user contact types and variables
c     [show,dict] show dictionary of program array allocation
c     [show,elem] show user element types
c     [show,mate] show material type use from last solution
c     [show,part] show partition data.

      elseif(j.eq.3) then

c       Output loaded user element descriptors

        if(pcomp(lct,'elem',4)) then
          write(*,2007)
          do nn = 1,50,20
            do i = nn,min(nn+19,50)
              call elmlib(hr(np(25)),hr(np(41)),hr(np(44)),mr(np(33)),
     &                    hr(np(39)),hr(np(36)),hr(np(35)),ndf,ndm,nst,
     &                    i,0)
            end do ! i
          end do ! nn

c       Contact element descriptors

        elseif(pcomp(lct,'cont',4)) then

          call contact (200)

c       Dictionary prints

        elseif(pcomp(lct,'dict',4)) then

          call pprtd()

c       Material properties

        elseif(pcomp(lct,'mate',4)) then

          write(iow,2005)
          if(ior.lt.0) then
            write(*,2005)
          end if

          do i = 1,10
            if(max(nomats(1,i),nomats(2,i)).gt.0) then
              write(iow,2006) i,nomats(1,i),nomats(2,i)
              if(ior.lt.0) then
                write(*,2006) i,nomats(1,i),nomats(2,i)
              end if
            end if
          end do ! i

          do i = 1,10
            if(max(unmats(1,i),unmats(2,i)).gt.0) then
              write(iow,2006) i+10,unmats(1,i),unmats(2,i)
              if(ior.lt.0) then
                write(*,2006) i+10,unmats(1,i),unmats(2,i)
              end if
            end if
          end do ! i

c       Partition properties

        elseif(pcomp(lct,'part',4)) then

          call shpart()

c       Output array values

        elseif(.not.pcomp(lct,'    ',4)) then

          call outary(lct,ct)

c       Show problem sizes

        else
          if(neq.gt.0) then
            npr = mr(np(20+npart)+neq-1)
            nnp = npr/neq
          else
            npr = 0
            nnp = 0
          end if
          if(nstep.gt.0) then
            ttot = titer
            avg  = dble(titer)/dble(nstep)
          else
            ttot = niter
            avg  = niter
          end if
          if(ior.lt.0) then
            write(*,2008) ndm,ndf,numnp,numel,nummat,neq,
     &                    geqnum,gpart,npr,nnp,
     &                    nrbody,numjts,ttim,rnmax,dt,aengy,tol,resnm,
     &                    prop,augf,npart,noi,nstep,pfr,avg,ttot
          end if
          write(iow,2008) ndm,ndf,numnp,numel,nummat,neq,
     &                    geqnum,gpart,npr,nnp,
     &                    nrbody,numjts,ttim,rnmax,dt,aengy,tol,resnm,
     &                    prop,augf,npart,noi,nstep,pfr,avg,ttot
        end if

c     [scre]en,<on,off> Set plot screen on/off

      elseif(j.eq.4) then

        if(pcomp(lct,'off',3)) then
          screfl = .false.
          if(ior.lt.0) write(*,2009)
        else
          screfl = .true.
          if(ior.lt.0) write(*,2010)
        end if

c     [comm]ent,<message> Echo comment to screen when in batch mode

      elseif(j.eq.5) then
        if(ior.gt.0) write(*,2011) lct

c     [smoo]th,,<number> Smooth coordinates for mesh

      elseif(j.eq.6) then

        nnp = max(1,nint(ct(1)))
        call smooth(hr(np(43)),nnp)

c     [set],temp,dof       - Move 'dof' to 'T' array
c     [set],disp,node,dof  - Move 'node,dof' to 'disp' array
c     [set],comp,#no       - Set component #no for outputs

      elseif(j.eq.7) then

c       Set 'T' array to dof solutions

        if(pcomp(lct,'temp',4)) then
          i = max(1,min(ndf,nint(ct(1))))

          write(iow,2004) i
          if(ior.lt.0) then
            write(*,2004) i
          end if

          nnp = 0
          do nn = i-1,nneq-1,ndf
            hr(np(38)+nnp) = hr(np(40)+nn)
            nnp            = nnp + 1
          end do ! nn

c       Set current solution for node 'nnp', dof = nn into specified
c       displacement array

        elseif(pcomp(lct,'disp',4)) then

          nnp       = max(1,min(numnp,nint(ct(1))))
          nn        = max(1,min(ndf,nint(ct(2))))
          fp(1)     = np(27) + ndf*(numnp + nnp - 1) + nn - 1
          fp(2)     = np(40) + ndf*(nnp - 1) + nn - 1
          hr(fp(1)) = hr(fp(2))
          if(pfr) then
            write(iow,2013) nn,nnp,hr(fp(2))
          endif
          if(ior.lt.0) then
            write(  *,2013) nn,nnp,hr(fp(2))
          endif

c       Set output component to 'nocomp'

        elseif(pcomp(lct,'comp',4)) then
          nocomp = max(1,min(ndf,nint(ct(1))))
          if(pfr) then
            write(iow,2012) nocomp
          endif
          if(ior.lt.0) then
            write(  *,2012) nocomp
          endif

c       Error for a 'set' command

        else
          write(iow,3000) lct
          write(ilg,3000) lct
          if(ior.lt.0) then
            write(*,3000) lct
          end if
        end if

c     [assi],name,num,value - Assign 'value' to 'name(num)'.

      elseif(j.eq.8) then

        call pgetd(lct,fp(1),nlg,npr,aflag)

        if(aflag) then
          nnp = max(1,nint(ct(1))) - 1
          if(npr.eq.1) then
            mr(fp(1)+nnp) = nint(ct(2))
          elseif(npr.eq.2) then
            hr(fp(1)+nnp) = ct(2)
          else
            write(iow,3001) lct,nnp
            write(ilg,3001) lct,nnp
          endif
        else
          write(iow,3001) lct
          write(ilg,3001) lct
        endif

c     [outp]ut,array,- Output arrays for use by 'matlab'
c                      sparse matrix formats
c     array = [tang,utan,lmas,mass,cmas,umas,damp,cdam,udam,dr,form]
c     filename = same as array specified
      elseif(j.eq.9) then

        call pouta(lct)

      end if

c     Formats

2004  format(5x,'-> Solution in degree-of-freedom',i3,
     &          ' moved to T array')
2005  format('   M a t e r i a l   C o u n t s'//
     &       7x,'Type   1-Calls   2-Calls')

2006  format(3i10)

2007  format('   A v a i l a b l e    U s e r    E l e m e n t    ',
     &       'T y p e s',/)

2008  format(/,
     &  '   C u r r e n t    S o l u t i o n    P a r a m e t e r s',/
     &  /,'     Mesh Dimension    =',i12  ,
     &    '  :  Number Dof/Node   =',i12,/
     &    '     Number Nodes      =',i12  ,
     &    '  :  Number Elements   =',i12,/
     &    '     Number Materials  =',i12  ,
     &    '  :  Number Equations  =',i12,/
     &    '     No. Global Eqs.   =',i12  ,
     &    '  :  Global Partition  =',i12,/
     &    '     Profile Terms     =',i12  ,
     &    '  :  Average Column Ht =',i12,/
     &    '     No. Rigid Bodies  =',i12  ,
     &    '  :  Number Joints     =',i12,/
     &    '     Time              =',1p,1e12.4,
     &    '  :  Max. Energy Norm  =',1p,1e12.4,/
     &    '     Dt                =',1p,1e12.4,
     &    '  :  Energy Norm       =',1p,1e12.4,/
     &    '     Solution Tol.     =',1p,1e12.4,
     &    '  :  Rel.Residual Norm =',1p,1e12.4,/
     &    '     Proportional Load =',1p,1e12.4,
     &    '  :  Augment Factor    =',1p,1e12.4,/
     &    '     Partition Number  =',i12  ,
     &    '  :  Time Integration  =',i12,/
     &    '     No. Time Steps    =',i12  ,
     &    '  :  Command Printing  =',l12,/
     &    '     Avg. Iter/Step    =',1p,1e12.4,
     &    '  :  No. Iterations    =',i12,/)

2009  format('   Plot Screen = OFF'/)
2010  format('   Plot Screen = ON'/)
2011  format('   Solution at: ',a)

2012  format('   Output component =',i3/)

2013  format('   Specified displacement d(',i2,',',i8,')  =',1p,1e15.6)

3000  format(' *WARNING* PMACR5: No match on SET command: Option = ',a)

3001  format(' *ERROR* PMACR5: No match on ASSIgn to array :',a:,
     &       ' for element number =',i10)

      end
