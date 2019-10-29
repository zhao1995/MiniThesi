c$Id: umacr9.f,v 1.1 2008/12/12 23:58:47 rlt Exp $
      subroutine umacr9(lct,ctl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2009: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Remove 'prt' from argument list                  09/07/2009
c       2. Inserted two damage variables (lines 242-264)    15/04/2018
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  User interface for adding solution command language
c                instructions.

c      Inputs:
c         lct       - Command character parameters
c         ctl(3)    - Command numerical parameters

c      Outputs:
c         N.B.  Users are responsible for command actions.  See
c               programmers manual for example.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'comfil.h'
      include  'counts.h'
      include  'iofile.h'
      include  'pointer.h'
      include  'sdata.h'
      include  'cdat1.h' !for 'nie'
      include  'strnum.h'
      include  'umac1.h'
      include  'hdata.h'
      include  'tdata.h'
      include  'swvars.h'

      include  'comblk.h'

      logical   pcomp
      character lct*15
      character fname*128
      integer*8 i,ii,node,h2Pointer,iHVar,h2PtrOffs,iVtk,iNp, idf,j,k
      real*8, intent(inout)::    ctl(3)
      integer,parameter::numHVars=9!number of history vars
      integer,parameter::nIntPts=4!nmb int pts per el
      character namesHVars(numHVars)*25
      integer   dimHVars(numHVars), nScalHVars
      logical setvar,ualloc
      character(len=20):: filename
      integer xpoint,xlen,xpre
      logical flg
      real*8 eAv,omegaIP(numel,nIntPts),omegaEL(numel),omegaTot
      real*8 dblDmy, dblDmy2
!       real*8 epsBTBar(3,3),xNds(numnp,ndm),uBar(ndm*numnp)
!       integer mat2plot
      real*8, parameter::w22=dsqrt(2.d0)/2.d0
      
      double precision, dimension(:), allocatable::hVrsAvg
      
      save


!       namesHVars(1)='lambdahat';     dimHVars(1)=4
!       namesHVars(2)='udelta';    dimHVars(2)=12
!       namesHVars(3)='DD';        dimHVars(3)=4
!       namesHVars(4)='tt';       dimHvars(4)=4
!       namesHVars(2)='udelta';    dimHVars(2)=12
!       namesHVars(3)='DD';        dimHVars(3)=4
!       namesHVars(4)='tt';       dimHvars(4)=4
      
      namesHVars(1)='Omegap';   dimHVars(1)=1
      namesHVars(2)='S2PKMatr';     dimHVars(2)=6
      namesHVars(3)='PiH';     dimHVars(3)=1
      namesHVars(4)='PP';     dimHVars(4)=1
      namesHVars(5)='DMatr';     dimHVars(5)=1
      namesHVars(6)='E_GL';     dimHVars(6)=6
      namesHVars(7)='inelStep';     dimHVars(7)=1
      namesHVars(8)='DFbr';     dimHVars(8)=1
      namesHVars(9)='inelStpFb';     dimHVars(9)=1
      
!       write(*,*) 'ttim=', ttim
      nScalHVars=0
      do i=1,numHVars
        nScalHVars = nScalHVars+dimHVars(i)
      end do
      
!       call cpu_time(dblDmy)
!       write(*,*) 'cpu time=', dblDmy-uuSW
      allocate(hVrsAvg(nScalHVars))

!       write(*,*) ctl(1)
!       ctl(1)=ctl(1)+1

!       write(*,*) 'NOTE: Using PVIE requires to possibly modify'
!       write(*,*) '  macr9.f according to your history variables.'

c     Set command word
!       write(*,*) 'ndf=', ndf
!       do iNp=0,numnp-1
!         if( dabs(hr( np(43) + iNp*ndm )-12000.d0).lt.1.d-5 ) then
!           do idf=0,ndf-1
!             write(*,*) hr( np(26) + iNp*ndf + idf )
!           end do
!         endif
!       end do
  
      if(pcomp(uct,'mac9',4)) then      ! Usual    form
        uct = 'pvie'  
        write(*,*) 'iiSW used by umacr9'
        iiSW=0        
!         write(*,*) 'Input file requires line <<rpod init>>'
!         read(*,*)
        ! Specify 'name'
!         call palloc(239, 'IVTKR', 8, 1)
      elseif(urest.eq.1) then           ! Read  restart data

      elseif(urest.eq.2) then           ! Write restart data

      else                              ! Perform user operation
        iiSW=iiSW+1
!       call pgetd('IVTK ', xpoint, xlen, xpre, flg)
!       write(*,*) 'IVTK ', xpoint, xlen, xpre, flg
!       if(.not.flg) setvar=ualloc(2, 'IVTK ', 1, 1)
!       mr(np(336))=0
      
      
      if(.not.pcomp(lct,'time',4)) then
!         mr(up(2)) = mr(up(2))+1
!         iVtk=mr(up(2))
        iVtk=iiSW
      if(iVtk<10) then
        write(filename, "(A4,I1,A4)") lct, iVtk,".vtu"
      elseif(iVtk<100) then
        write(filename, "(A4,I2,A4)") lct, iVtk,".vtu"
      elseif(iVtk<1000) then
        write(filename, "(A4,I3,A4)") lct, iVtk,".vtu"
      elseif(iVtk<10000) then
        write(filename, "(A4,I4,A4)") lct, iVtk,".vtu"
      else
        write(*,*) 'ERROR. Can not do more than 9999 output files.'
        stop
      endif !if(iVtk<10)
! 	open(unit=14, file=filename, status="unknown")
        open(unit=99,file=filename,access='sequential')
        write(*,*) 'Saving PARAVIEW data to ',filename,'at t=',ttim
!         write(*,*) 'WARNING: hist var. avg for square or cubic els.'
      else
        i = index(fplt,' ')
        fname(1:128) = ' '
        fname(1:i-1)=fplt(1:i-1)
        fname(i:i+4) = '00000'
        if (nstep.le.9) then
          write(fname(i+4:i+4),'(i1)') nstep
        elseif (nstep.le.99) then
          write(fname(i+3:i+4),'(i2)') nstep
        elseif (nstep.le.999) then
          write(fname(i+2:i+4),'(i3)') nstep
        elseif (nstep.le.9999) then
          write(fname(i+1:i+4),'(i4)') nstep
        elseif (nstep.le.99999) then
          write(fname(i:i+4),'(i5)') nstep
        endif
        call addext(fname,'vtu',128,3)
        open(unit=99,file=fname,access='sequential')
        write(*,*) 'Saving PARAVIEW data to ',fname,'at t=',ttim
      end if

c     Write out top header 
      write(99,1000) 

      write(99,1010) numnp,numel               ! Start Mesh/Piece Section
      write(99,1020)                           ! Start Point/Node data
      write(99,1030) 'Float64','nodes',3

      do i = 1,numnp
       do ii = 1,ndm
         write(99,5000) hr(np(43)+(i-1)*ndm+(ii-1))
       end do
       if (ndm==2) then
        write(99,5000) 0.d0
       endif
      end do

      write(99,1031)                           ! Close Node data
      write(99,1021)                           ! Close Points section

      write(99,1050)                           ! Start Cell Section
      write(99,1030) 'Int32','connectivity',1  ! Start Elements

c     Offsets memory allocation
      call palloc(111,'TEMP1',numel+1,1)
      mr(np(111)) = 0;
!       call palloc(113, 'TEMP3', 1, 1)


      do i = 1,numel
       mr(np(111)+i) = mr(np(111)+i-1)
       do ii = 1,nen
        node = mr(np(33) + ii-1 + nen1*(i-1))
        if (node .ne. 0) then
           write(99,6000) node-1
           mr(np(111)+i) = mr(np(111)+i) + 1
        endif
       end do
      end do

      write(99,1031)                          ! Close Elements


      write(99,1030) 'Int32','offsets',1      ! Start Offsets
      do i = 1,numel
       write(99,6000) mr(np(111)+i)
      end do

      write(99,1031)                           ! Close Offsets

      write(99,1030) 'UInt8','types',1         ! Start Element types
      do i = 1,numel
       if (mr(np(111)+i)-mr(np(111)+i-1) .eq. 2) then          ! 2 node line
           write(99,6000) 3
       elseif (mr(np(111)+i)-mr(np(111)+i-1) .eq. 3) then      ! 3 node triangle
           write(99,6000) 5
       elseif (mr(np(111)+i)-mr(np(111)+i-1) .eq. 4) then      ! 4 node quad
           write(99,6000) 9
       elseif (mr(np(111)+i)-mr(np(111)+i-1) .eq. 8) then      ! 8 node brick
           write(99,6000) 12
       elseif (mr(np(111)+i)-mr(np(111)+i-1) .eq. 27) then      ! 27 node brick (plot as 20 node)
           write(99,6000) 25
       endif
      end do
      call palloc(111,'TEMP1',0,1)

      write(99,1031)                             ! Close Element types
      write(99,1051)                             ! Close Cell Section

      write(99,1060)                             ! Start Point Data

      write(99,1030) 'Float64','disp', 3 !ndm    ! Start Displacements
      do i = 1,numnp
       do ii =0,ndf-3
         write(99,5000) hr(np(40)+(i-1)*ndf+ii)
       end do
       if (ndm==2) then
       write(99,5000) 0 
       endif
       !write(99,5000) hr(np(40)+(i-1)*ndf+2)
! ! ! ! ! !        write(99,5000) hr(np(40)+(i-1)*ndf+ndf-1)
       if (ndf==2) then
         write(99,5000) 0.d0
       endif
      end do
      write(99,1031)                             ! Close Displacements


      write(99,1030) 'Float64','dam',2        ! Start Damage
      do i = 1,numnp
       do ii =ndf-2,ndf-1
         write(99,5000) hr(np(40)+(i-1)*ndf+ii)
       end do
      end do
      write(99,1031)                             ! Close Damage
      
      
      
      
      
      
      
      
      
C C       TEST STUFF
C       write(99,1030) 'Int32','intP',1 
C       do i = 1,numnp
C       write(99,6000) i
C       end do
C       write(99,1031)


!       write(*,*) 'nh1=', nh1, 'nh2=', nh2, 'ht1=', ht1, 'ht2=', ht2
      write(99,1061)                             ! Close Point Data Section

      write(99,1070)				!Cell Data
      h2PtrOffs=0
      do i=1,nummat
!        write(*,*) "Offs h.-vars(ma=",i,"):",
!      1  mr(np(32)-1+(i-1)*nie+nie-3)      
      end do
      omegaTot=0.d0
      do iHVar=1,nScalHVars; hVrsAvg(iHVar)=0.d0; end do
      do iHVar=1,numHVars
	write(99,1030) 'Float64',namesHVars(iHVar),dimHVars(iHVar)        ! Start 
	do i = 0,(numel-1)
          h2Pointer=np(49)+mr(np(33)-1+i*nen1+nen+2)
     1     +h2PtrOffs!49:address of history-vars, 33:address of ix
!           do j=1,nIntPts
!             omegaIP(j)=hr(h2Pointer+6+nScalHVars*(j-1))
!           end do
          if(iHVar.eq.1) then
            omegaEL(i+1)=0.d0
            do j=1,nIntPts
              omegaIP(i+1,j)=hr(h2Pointer+nScalHVars*(j-1))
              omegaEL(i+1)=omegaEL(i+1)+omegaIP(i+1,j)
            end do
            omegaTot=omegaTot+omegaEL(i+1)
          endif
! 	write(*,*) "offs. h-Vars.:", mr(np(33)-1+i*nen1+nen+2)
          do ii =0,dimHVars(iHVar)-1
            eAv=0.d0
            do j=0,nIntPts-1
              eAv=eAv+hr(h2Pointer+ii+nScalHVars*j)*omegaIP(i+1,j+1)
            end do
            hVrsAvg(1+h2PtrOffs+ii) = hVrsAvg(1+h2PtrOffs+ii)+eAv
            eAv=eAv/omegaEL(i+1)
!             if(mr(np(33)+nen1*(i+1)-1).ne.mat2plot) then
!               eav=0.d0
!             end if
            if(dabs(eAv).lt.1.d-89) eAv=0.d0!1.d-16
            write(99,5000) eAv
            
          end do         
	end do
	
	write(99,1031)                             ! Close 
	h2PtrOffs = h2PtrOffs+dimHVars(iHVar)
      end do
      do iHVar=1,nScalHVars
        hVrsAvg(iHVar)=hVrsAvg(iHVar)/omegaTot
      end do
!       write(*,*) 'omegaTot=', omegaTot
      write(99,1071)				!Cell Data

      write(99,1011)                             ! Close Mesh/Piece


c     Close the XML file
      write(99,1001)

      
      
      
!       do i=1,3; epsBTBar(i,i)=hVrsAvg(i); end do
!       epsBTBar(1,2)=hVrsAvg(6)*w22; epsBTBar(1,3)=hVrsAvg(5)*w22
!       epsBTBar(2,3)=hVrsAvg(4)*w22; epsBTBar(2,1)=epsBTBar(1,2)
!       epsBTBar(3,1)=epsBTBar(1,3); epsBTBar(3,2)=epsBTBar(2,3)
!       
!       do i=0,numnp-1; do j=1,ndm;
!         uBar(ndm*i+1:ndm*i)=matmul( epsBTBar(1:ndm,1:ndm),xNds(i+1,:) )
!       end do; end do
!       open(unit=294,file='AFT.dat',access='append')
!       write(294,*) hr(np(40):np(40)+numnp*ndf-1) - uBar
!       close(294)
!       
!       
!       open(unit=384, file="nSnapshots.dat", status='old', action='read')
!       read(384,*) i
!       close(384)
!       open(unit=462,file='nSnapshots.dat',access='sequential')
!       write(462,*) i+1
!       close(462)
      
      
      
      
      
      
      
      
      close(99)
!       if(iVtk.eq.1) then
!         open(unit=98,file='histVars.dat',access='sequential')
!         write(98,*) ' '
!         close(98)      
!       endif
!       open(unit=98,file='histVars.dat',access='append')
!       write(98,5000) ttim
!       do iHVar=1,nScalHVars; eAv=hVrsAvg(iHVar); write(98,5000) eAv
!       end do
!       write(98,*) '/'
!       close(98)
      
      
      
      
      
      
! ! ! ! !       
! ! ! ! !       
! ! ! ! !       do i = 1,numnp
! ! ! ! !        do ii = 1,ndm
! ! ! ! !          write(99,5000) hr(np(43)+(i-1)*ndm+(ii-1))
! ! ! ! !        end do
! ! ! ! !       end do
! ! ! ! !       
! ! ! ! !       
! ! ! ! !        do i = 1,numel
! ! ! ! ! !        mr(np(111)+i) = mr(np(111)+i-1)
! ! ! ! !        do ii = 1,nen
! ! ! ! !         node = mr(np(33) + ii-1 + nen1*(i-1))
! ! ! ! !         if (node .ne. 0) then
! ! ! ! !            write(99,6000) node-1
! ! ! ! ! !            mr(np(111)+i) = mr(np(111)+i) + 1
! ! ! ! !         endif
! ! ! ! !        end do
! ! ! ! !       end do
      
      
      
!       
!       
!       open(unit=982,file='rgnsStrctrd.dat',access='sequential')
!       do i=1,numel
!         dblDmy=0.d0; dblDmy2=0.d0
!         do j=1,4
!           node=mr(np(33) + j-1 + nen1*(i-1))
!           dblDmy=dblDmy+hr(np(43)+(node-1)*ndm)
!           dblDmy2=dblDmy2+hr(np(43)+(node-1)*ndm+1)
!         end do
!         dblDmy=dblDmy/4.d0; dblDmy2=dblDmy2/4.d0
!         k=int(dblDmy*4.d0)+4*int(dblDmy2*4.d0)
!         write(982,*) k
!         write(982,*) k
!         write(982,*) k
!         write(982,*) k
!       end do      
!       close(982)
!       
      
      endif

1000  format('<?xml version="1.0"?>',/
     * '<VTKFile type="UnstructuredGrid" version="0.1">',/
     * '<UnstructuredGrid>')
1001  format('</UnstructuredGrid> </VTKFile>')


1010  format('<Piece NumberOfPoints="',i10,'" NumberOfCells="',i10,'">')
1011  format('</Piece>')


1020  format('<Points>')
1021  format('</Points>')

1030  format('<DataArray type="',a,'" 
     * Name="',a,'"
     * NumberOfComponents="',i2,'" format="ascii">')
1031  format('</DataArray>')

1060  format('<PointData>')
1061  format('</PointData>')

1070  format('<CellData>')
1071  format('</CellData>')

1050  format('<Cells>')
1051  format('</Cells>')

5000  format(e15.5,' ',$)
6000  format(i6,' ',$)

      deallocate(hVrsAvg)

      end
