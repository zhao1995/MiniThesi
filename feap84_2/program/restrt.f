c$Id:$
      subroutine restrt(fres,ndm,ndf,nneq,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Dimension 'fres*(*)'                             12/10/2007
c       2. Add base proportional load 'F0'                  29/12/2008
c       3. Add 'prldv' array and 'propo'                    14/07/2009
c       4. Add 'tx' to argument on call to umshlib          26/09/2011
c       5. Check 20 user macros                             25/01/2012
c       6. Increase allocation for F0 to 4*nneq             10/01/2013
c       7. Add 'frotas' /crotas.h for rotational parameters 20/08/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Read/save restart files for resolutions

c      Inputs:
c         fres    - Name of restart file to read/save
c         ndm     - Spatial dimension of mesh
c         ndf     - Number dof/node
c         nneq    - Total dumber of parameters in solutions
c         isw     - Switch: = 1 for read; =2 for save.

c      Outputs:
c         u(*)    - Solution state read
c         none    - from/to pointers
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'arclel.h'
      include  'arcler.h'
      include  'cdata.h'
      include  'counts.h'
      include  'crotas.h'
      include  'ddata.h'
      include  'dyndat.h'
      include  'evdata.h'
      include  'fdata.h'
      include  'gltran.h'
      include  'iodata.h'
      include  'iofile.h'
      include  'ndata.h'
      include  'part7.h'
      include  'pointer.h'
      include  'print.h'
      include  'prlod.h'
      include  'prld1.h'
      include  'rdata.h'
      include  'tdata.h'
      include  'tdato.h'
      include  'umac1.h'
      include  'comblk.h'

      include  'p_point.h'

      logical   exst,sfl,fl9,setvar,ralloc,walloc,flrota
      character fres*(*),yorn*1,lct*15, tx(8)*15
      integer   i,ndm,ndf,nneq,isw, nnpo,nnlo,nnmo,ndmo,ndfo,nh2
      integer   mmo,mmt,mmr
      real*8    ctl(3)

      save

      data      ctl /3*0.0d0/

c     Check file status

1     inquire(file=fres,exist=exst)
      if(.not.exst.and.isw.eq.1) then
        write(iow,3002) fres
        if(ior.lt.0) then
          write(*,3002) fres
          call pprint(
     &       '           Specify new name for restart file? (y or n) >')
          read (*,1000) yorn
          if(yorn.eq.'y' .or. yorn.eq.'Y') then
            call pprint('           New Restart File Name >')
            read (*,1000) fres
            goto  1
          endif
        endif
        return
      endif

c     Open file

      if(exst) then
        open (ios,file=fres,form='unformatted',status='old')
      else
        open (ios,file=fres,form='unformatted',status='new')
      endif
      rewind ios

c     Read restart files

      if(isw.eq.1) then

c       Control information

        read(ios) nnpo,nnlo,nnmo,ndmo,ndfo,fl(9)
        if((nnpo.eq.numnp).and.(nnlo.eq.numel).and.(nnmo.eq.nummat)
     &          .and.(ndmo.eq.ndm).and.(ndfo.eq.ndf)) then

c         Solution parameters

          read(ios) theta,nrk,nrc,nrm,nrt,noi,numint,alpha,gtan,
     &              nstep,niter,naugm,titer,taugm,iaugm,iform,
     &              ttim,dt,dtold,rnmax,prop,rlnew,c0,cs01,cs02,
     &              ds0,r,det0,xn,fl9,mf,mq,flrota

c         Input displacement solution

          setvar = ralloc( 40,'U    ',ios)
          write(iow,2000) 'I n p u t',nstep,ttim,dt,'input'
          if(ior.lt.0) then
            write(*,2000) 'I n p u t',nstep,ttim,dt,'input'
          endif

c         Eigenpairs

          if(mq.gt.0) then
            setvar = ralloc( 76,'EVAL ',ios)
            setvar = ralloc( 77,'EVEC ',ios)
            write(iow,2001) 'input',mf,mq
            if(ior.lt.0) then
              write(*,2001) 'input',mf,mq
            endif
          endif

c         Transient data

          write(iow,2002) prop,rlnew
          if(ior.lt.0) then
            write(*,2002) prop,rlnew
          endif
          if(fl9) then
            setvar = ralloc( 42,'VEL  ',ios)
            write(iow,2003) 'input',noi
            if(ior.lt.0) then
              write(*,2003) 'input',noi
            endif
          endif

c         Current load state

          setvar = ralloc( 30,'FTN  ',ios)
          write(iow,2004) 'input'
          if(ior.lt.0) then
            write(*,2004) 'input'
          endif

c         Base proportional load (newf) and prldv/propo

          setvar = ralloc( 28,'F0   ',ios)
          read(ios) prldv,propo
          write(iow,2006) 'input'
          if(ior.lt.0) then
            write(*,2006) 'input'
          endif

c         History data

          setvar = ralloc( 49,'H    ',ios)
          write(iow,2005) 'input'
          if(ior.lt.0) then
            write(*,2005) 'input'
          endif
          refl   = .true.
          fl(11) = .false.

c         Contact data

          call contact (306)

c         Rotational parameters

          if(flrota) then
            setvar = ralloc( 81,'MO   ',ios)
            setvar = ralloc( 82,'MR   ',ios)
            setvar = ralloc( 83,'MT   ',ios)
            write(iow,2007) 'input'
            if(ior.lt.0) then
              write(*,2007) 'input'
            endif
          endif

c         User data

          urest = 1
          lct   = 'restart'
          do i = 1,10
            call umshlib(i,tx,prt)
          end do ! i
          do i = 1,20
            call umaclib(i,lct,ctl,prt)
          end do ! i
          urest = 0

        else
          write(iow,3001)
          if(ior.lt.0) then
            write(*,3001)
          endif
        endif

c     Save information for restart

      elseif(isw.eq.2) then

c       Control information

        write(ios) numnp,numel,nummat,ndm,ndf,fl(9)

c       Solution parameters

        fl9 = flp(9,1) .or. flp(9,2) .or. flp(9,3) .or. flp(9,4)
     &                 .or. fl(9)
        write(ios) theta,nrk,nrc,nrm,nrt,noi,numint,alpha,gtan,
     &             nstep,niter,naugm,titer,taugm,iaugm,iform,
     &             ttim,dt,dtold,rnmax,prop,rlnew,c0,cs01,cs02,
     &             ds0,r,det0,xn,fl9,mf,mq,frotas

c       Solution state

        setvar = walloc( 40,'U    ',nneq*3,2, ios)
        write(iow,2000) 'O u t p u t',nstep,ttim,dt,'output'
        if(ior.lt.0) then
          write(*,2000) 'O u t p u t',nstep,ttim,dt,'output'
        endif

c       Eigenpairs

        if(mq.gt.0) then
          setvar = walloc( 76,'EVAL ', mq    , 2, ios)
          setvar = walloc( 77,'EVEC ', mq*neq, 2, ios)
          write(iow,2001) 'output',mf,mq
          if(ior.lt.0) then
            write(*,2001) 'output',mf,mq
          endif
        endif

c       Transient data

        if(fl9) then
          setvar = walloc( 42,'VEL  ',nrt*nneq,2, ios)
          write(iow,2003) 'output',noi
          if(ior.lt.0) then
            write(*,2003) 'output',noi
          endif
        endif

c       Current load state

        setvar = walloc( 30,'FTN  ', 4*nneq, 2, ios)
        write(iow,2004) 'output'
        if(ior.lt.0) then
          write(*,2004) 'output'
        endif

c       Base proportional load (newf) and prldv/propo

        setvar = walloc( 28,'F0   ',4*nneq, 2, ios)
        write(ios) prldv,propo
        write(iow,2006) 'output'
        if(ior.lt.0) then
          write(*,2006) 'output'
        endif

c       History data

        call pgetd('H   ',point,nh2, i,sfl)
        if(.not.sfl) then
          nh2 = 1
        endif
        setvar = walloc( 49,'H    ', nh2, 2, ios)
        write(iow,2005) 'output'
        if(ior.lt.0) then
          write(*,2005) 'output'
        endif

c       Contact data

        call contact (307)

c       Rotational parameters

        if(frotas) then
          call pgetd('MO  ',point,mmo , i,sfl)
          if(sfl) then
            setvar = walloc( 81,'MO   ', mmo, 1,ios)
          else
            write(*,*) ' --> ERROR: MO =',mmo,i
          endif
          call pgetd('MR  ',point,mmr, i,sfl)
          if(sfl) then
            setvar = walloc( 82,'MR   ', mmr, 2,ios)
          else
            write(*,*) ' --> ERROR: MR =',mmr,i
          endif
          call pgetd('MT  ',point,mmt, i,sfl)
          if(sfl) then
            setvar = walloc( 83,'MT   ', mmt, 2,ios)
          else
            write(*,*) ' --> ERROR: MT =',mmt,i
          endif
          write(iow,2007) 'output'
          if(ior.lt.0) then
            write(*,2007) 'output'
          endif
        endif

c       User data

        urest = 2
        lct   = 'restart'
        do i = 1,10
          call umshlib(i,tx,prt)
        end do ! i
        do i = 1,20
          call umaclib(i,lct,ctl,prt)
        end do ! i
        urest = 0

      endif

c     Close file

      close(ios)

c     Formats

1000  format(a)

2000  format('   R e s t a r t   ',a,'   D a t a'/
     &       10x,'Time step number  =',i8/
     &       10x,'Time at restart   =',1p,1e12.5/
     &       10x,'Time increment    =',1p,1e12.5/
     &       10x,'Displacements ',a)

2001  format(10x,'Eigenpairs ',a,' for',i4,' modes',i4,' total values')

2002  format(10x,'Proportional load =',1p,1e12.5/
     &       10x,'Arc-length   load =',1p,1e12.5)

2003  format(10x,'Transient states ',a,' (noi =',i2,')')

2004  format(10x,'Force vector ',a)

2005  format(10x,'History data ',a)

2006  format(10x,'Prop load data ',a)

2007  format(10x,'Rotational parmeters ',a)

3001  format(' *ERROR* Incorrect information in a restart')

3002  format(' *ERROR* Restart file ',a17,' does not exist')

      end
