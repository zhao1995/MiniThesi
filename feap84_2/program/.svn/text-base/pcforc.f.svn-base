c$Id:$
      subroutine pcforc(x,f,ntyp,ndm,ndf,numnp,prt,prth,name,fext)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       Increase allowable ndm + ndf to 30                  05/07/2007
c       1. Add read of 'ldnum,ldprp,ldflg'; 'pload1.h'      09/03/2009
c          Add read of 'spnum,spflg'
c       2. Set length of file extension to 8                25/01/2011
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set nodal forces or displacments based on coordinate
c               inputs

c      Inputs:
c         x(ndm,*)  - Nodal coordinates
c         ntyp(*)   - Node type (negative for inactive node)
c         ndm       - Spatial dimension of mesh
c         ndf       - Number dof/node
c         numnp     - Number of nodes in mesh
c         prt       - Output generated results if true
c         prth      - Output title/header data if true
c         name      - Name of basic file to use for data
c         fext      - Extender, used to define data type

c      Outputs:
c         f(ndf,*)  - Generated force or displacement values
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'comfil.h'
      include  'cdat2.h'
      include  'conval.h'
      include  'iodata.h'
      include  'iofile.h'
      include  'ioincl.h'
      include  'pload1.h'
      include  'trdata.h'

      character type*4,name*5,fnam*132,fext*8,ftype*15
      logical   prt ,prth, errck,tinput,pcomp, lsave, trflg, clflg
      logical   oprt,oprth
      integer   ndm,ndf,numnp, iosave,i, n,nbc, ntyp(*)
      real*8    dotx,xmn, tmn,xl(3),f(ndf,numnp),x(ndm,numnp),td(31)
      real*8    sindeg,cosdeg

      save

      data      trflg /.false./

c     Open file for reads

      fnam = fsav
      call addext(fnam,fext,128,8)
      call opnfil(fext,fnam,-2,ios,lsave)
      if(lsave) then
        iosave = ior
        ior    = ios

        do i = 0,36
          do n = 1,26
            vvsave(n,i) = vvv(n,i)
          end do ! n
        end do ! i
        do i = 1,3
          x0sav(i) = x0(i)
        end do ! i
        oprt  = prt
        oprth = prth
        read(ior,1000) type,fincld(isf),irecrd(isf),prt,prth
        read(ior,1001) vvv
        read(ior,1001) tr,xr,trdet,x0
        read(ior,1002) ldnum,ldprp,spnum,ldflg,spflg
      else
        write(iow,3000)
        write(ilg,3000)
        call plstop()
      endif

c     Output the header

      call prtitl(prth)
      if(prt) write(iow,2000) (n,name,n=1,ndf)

c     Read input of boundary coordinates for nodal forced values
c         td(i) - coordinates :   1 =< i =< ndm
c         td(i) - forced value: ndm  < i =< ndm + ndf

200   if(ior.lt.0) then
        write(*,3001) ndf
        call pprint('   >')
      endif
      call pzero(td,ndf+ndm)
      ftype = ' '
      i     = min(15,ndf + ndm)
      errck = tinput(ftype,1,td,i)
      if(errck) go to 200

c     Find closest node to input coordinates

      if(pcomp(ftype,'cart',4)) then
        trflg = .false.
      elseif(pcomp(ftype,'pola',4)) then
        trflg = .true.
      elseif(pcomp(ftype,'node',4)) then
        if(ndf.gt.15) then
          i     = ndf - 15
          errck = tinput(ftype,0,td(16),i)
        endif
        clflg = .false.
        do n = 1,numnp
          if(ntyp(n).ge.0) then
            if(trflg) then
              call pdegree(x(2,n), sindeg,cosdeg)
              xl(1) = x0(1) + x(1,n)*cosdeg
              xl(2) = x0(2) + x(1,n)*sindeg
              do i = 3,ndm
                xl(i) = x0(i) + x(i,n)
              end do ! i
            else
              do i = 1,ndm
                xl(i) = x(i,n)
              end do ! i
            endif
            tmn = dotx(td(1),xl,ndm)
            if(clflg) then
              if(tmn.lt.xmn) then
                xmn = tmn
                nbc = n
              endif
            else
              xmn   =  tmn
              nbc   =  n
              clflg = .true.
            endif
          endif ! ntyp(n) > 0
        end do ! n

c       Set forced values

        if(clflg) then
          if(pcomp(type,'add',3)) then
            do n = 1,ndf
              f(n,nbc)  = f(n,nbc) + td(ndm+n)
              td(ndm+n) = f(n,nbc)
            end do ! n
          else
            do n = 1,ndf
              f(n,nbc) = td(ndm+n)
            end do ! n
          endif

c         Output nodal forced values

          if(prt) write(iow,2001) nbc,(td(ndm+n),n=1,ndf)
        endif

        go to 200
      elseif(pcomp(ftype,'    ',4)) then
        close(ior,status = 'delete')
        ior = iosave
      else
        go to 200
      endif

      do i = 0,36
        do n = 1,26
          vvv(n,i) = vvsave(n,i)
        end do ! n
      end do ! i
      do i = 1,3
        x0(i) = x0sav(i)
      end do ! i
      prt  = oprt
      prth = oprth

c     Format

1000  format(a4,2x,a12,i8,2l5)
1001  format(1p,4e20.12)
1002  format(3i8,2l3)

2000  format('  C o o r d i n a t e    N o d a l    ',
     &      'V a l u e s'//4x,'node',6(i6,'-',a5):/(8x,6(i6,'-',a5)))

2001  format(i8,1p,6e12.4:/(8x,1p,6e12.4))

3000  format(' *ERROR* PCFORC: Surface loading file',a,' does not',
     &       ' exist')

3001  format(' Input: x(i),i=1,ndm),(f(i),i=1,',i2,')')

      end
