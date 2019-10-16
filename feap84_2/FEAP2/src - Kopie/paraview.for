      subroutine pvpost(x,ix,s,u,v,a,fldyn1,fldyn2,nen1,ndm,ndf,n1,n2,
     +n3)
c-----------------------------------------------------------------------
c
c.... Purpose: store data fields for postprocessing in ParaView Format
c
c     Inputs:
c     x(ndm,*)   - nodal coordinates
c     ix(nen1,*) - node-element-relation
c     s(numnp,*) - nodal stresses (not for beam elements)
c     u(ndf,*)   - nodal displacements
c     nen1       - nen+4 nen=maximum number of nodes on element
c     ndm        - coordinate dimension
c     ndf        - number d.o.f at each node (1-6)
c     n1         - type of element or time step
c     n3         - 1=init / 2=next / 3=eigv
c
c     npm(numnp)         = assign material to each node
c     npm1(numnp,nummat) = defines node numbers per material set
c                          0,1,2,... per material set
c
c     [parv,init,n1,n2] initialize relevant data  n3=1
c                       initialize Rfile.pvd (head file)
c                       write first (initial) time step to Rfile.mxx.t0001.vtu
c                       n2 = switch for file separation per material set
c                          separation = 0 (default)
c                          no separation <> 0 or nummat > 99
c
c     [parv,next,n1]    write all further time steps to Rfile.mxx.tyyyy.vtu
c                       n3=2
c
c     [parv,eigv,n1,n2] write eigenvector n2 to Rfile.EV'n2'.vtu
c                       n3=3
c
c     Outputs:
c     Rfile.pvd           collection file includes all material sets per timestep (head file)
c     Rfile.mxx.tyyyy.vtu data of each material set and timestep
c     Rfile.EVxx.vtu      data of eigenvector xx
c
c     ParaView 4.01 supports
c     eltyp numbered as given in 'VTK Formats for Version 4.2' by www.kitware.com)
c     1D: vtk_line                    (vtktyp= 3)
c     2D: vtk_triangle                (vtktyp= 5)
c         vtk_quad                    (vtktyp= 9)
c         vtk_quadratic_quad          (vtktyp=23)
c         vtk_biquadratic_quad        (vtktyp=28)
c     3D: vtk_tetra                   (vtktyp=10)
c         vtk_hexahedron              (vtktyp=12)
c         vtk_quadratic_hexahedron    (vtktyp=25)
c
c--------------------------------------------------------------------
c     Comments WW and open in paraview.for or paraview.exe
c            - for beams no stress resultants!
c            - Names of stresses/displacements + indices
c            - store data binary
c
c            - 16 node VTK_HIGHER_ORDER_QUAD                (vtktyp=62)
c            - 27 node VTK_TRIQUADRATIC_HEXAHEDRON          (vtktyp=29)
c            - 64 node brick- VTK_HIGHER_ORDER_HEXAHEDRON   (vtktyp=67)
c
c          27 node hexahedron(29) eingebaut, fast richtig, Fehler in paraview??
c                                 auskommentiert
c          64 node hexahedron(67) eingebaut, Bild f. Knotennummern fehlt,
c                                 geht nicht, auskommentiert.
c
c
c     S. Lauterbach IBS KIT 11/12
c
c     SLau/WW KIT 10/13 Modifikation def.mesh = u(1-3)
c     SLau    KIT 11/13 Modifikation Eigenvektor
c
c--------------------------------------------------------------------
      USE cdata
      USE comfil
      USE fileno
      USE iofile
      USE pnodn
      USE strnam
      USE subdt
      USE tdata
      implicit double precision(a-h,o-z)
      logical     ex,flgow,flgel,fldyn1,fldyn2
      character   antw*3,eltyp*69
      character   fname1*229,fname2*229,fname3*229,fname4*229
      integer     vtktyp(nummat)
      integer, allocatable, dimension (:) :: npm
      integer, allocatable, dimension (:,:) :: npm1
      dimension   x(ndm,numnp),ix(nen1,numel),s(numnp,*)
      dimension   u(ndf,numnp),v(ndf,numnp),a(ndf,numnp)
      dimension   eltyp(16)
      dimension   numnp1(nummat),numel1(nummat),node(nummat)
c
c

      save ityp,itime
      save flgow
      save iswm
c
      data ip1/23/
      data ip2/24/
      data ntyp/15/
c.... at present same as for tecplot
      data eltyp/
     1'plane stress ndf = 2  ndm = 2,3 nen = 3,4,8,9                 ',
     2'beam2d/axish ndf = 3  ndm = 2   nen = 2                       ',
     3'plate        ndf = 3  ndm = 2,3 nen = 3,4,8,9                 ',
     4'plain strain ndf = 3  ndm = 2   nen = 4,9   2D cosserat       ',
     5'brick        ndf = 3  ndm =   3 nen = 8,20,27                 ',
     6'brick        ndf = 4  ndm =   3 nen = 8,20,27 th.-mech. coupl.',
     7'brick        ndf = 6  ndm =   3 nen = 8,20,27 3d cosserat     ',
     8'shells5      ndf = 5  ndm =   3 nen = 3,4,9                   ',
     9'beam3d       ndf = 6  ndm =   3 nen = 2                       ',
     +'shells6      ndf = 6  ndm =   3 nen = 3,4,9                   ',
     1'beam3d       ndf = 7  ndm =   3 nen = 2     incl. warping     ',
     2'shell        ndf = 7  ndm =   3 nen = 4,9   incl. warping beam',
     3'plain strain ndf = 1  ndm =   2 nen = 4,9,16 2D phase field   ',
     4'plain strain ndf = 5  ndm =   2 nen = 4,9,16 2D phase field   ',
     5'brick        ndf = 7  ndm =   3 nen = 8,27 3D phase field     ',
     6'tetraeder    ndf = 3  ndm =   3 nen = 4                       '/
c
      mstv = iabs(istv) ! plot only mstv stress values
      if(n3.eq.1) then
        ityp = n1
        iswm = n2
        if(nummat.ge.100) then
          iswm = 1
          call drawmess('file separation for materials disabled
     + (limited to 99 sets)',1,0)
        end if
      end if
c
c.... zero arrays
      node   = 1
      numnp1 = 0
      numel1 = 0
c
      if(iswm.eq.0)then
        allocate (npm(numnp))
        allocate (npm1(numnp,nummat))
        npm    = 0
        npm1   = 0
        do i = 1,nummat,1
c....   number of elements per material set: numel1(i)
          do j = 1,numel,1
            if(ix(nen1,j).eq.i) then
              numel1(i)=numel1(i)+1
              do k = 1,nen,1
                if(ix(k,j).ne.0) npm(ix(k,j))=i
              end do
            end if
          end do
c....   number of nodes per material set: numnp1(i)
          do j = 1,numnp,1
            if(npm(j).eq.i) then
              numnp1(i) = numnp1(i) + npm(j)/i
              npm1(j,npm(j))=node(npm(j))
              node(npm(j))=node(npm(j))+1
            end if
          end do
        end do
      end if
c
      goto (1,2,3) n3
c
c.... initial state
1     itime = 1
      if(ityp.lt.1.or.ityp.gt.ntyp)then
        write(*,'(a23)') 'Element not implemented'
        write(*,'(a)')   'available elements'
        do i=1,ntyp
          write(*,'(i3,a,a)') i,' ',eltyp(i)
        end do
        return
      end if
c
c.... write .vtu file(s)
      if(iswm.ne.0)then
        fname1 = fres
        call vtufile(ip1,fname1,1,itime,flgow,iret,iswm)
        call pvdata_1(ip1,x,ix,s,u,v,a,fldyn1,fldyn2,mstv,1,numnp,
     +       numel,ityp)
        close(ip1)
      else
        do i = 1,nummat,1
          fname1 = fres
          call vtufile(ip1,fname1,i,itime,flgow,iret,iswm)
          if(iret.eq.1) return
          call pvdata(ip1,x,ix,s,u,v,a,fldyn1,fldyn2,mstv,i,npm,npm1,
     +    numnp1(i),numel1(i),ityp,iswm)
          close(ip1)
        end do
      end if
c
c.... write .pvd file  (define which files are handled together)
      fname2 = fres
      call pvdfile(ip2,fname1,fname2,itime,n3,iswm)
      return
c
c.... write all data for each material at time step
2     if(n1.ne.0) itime = n1
      itime = itime + 1
      if(itime.ge.10000) then
        call drawmess('Time steps only up to 9999 possible',1,0)
        return
      end if
c
      fname2 = fres
      call addpv(fname2,'pvd',3)
      inquire(FILE=fname2,EXIST=ex)
      if(.not.ex) then
        call drawmess(
     +  'PVD file needs to be initialized ! [parv,init,n1]',1,0)
        return
      end if
c
c.... write .vtu file(s)
      if(iswm.ne.0)then
        fname1 = fres
        call vtufile(ip1,fname1,1,itime,flgow,iret,iswm)
        call pvdata_1(ip1,x,ix,s,u,v,a,fldyn1,fldyn2,mstv,1,numnp,
     +                numel,ityp)
        close(ip1)
      else
        do i = 1,nummat,1
          fname1 = fres
          call vtufile(ip1,fname1,i,itime,flgow,iret,iswm)
          if(iret.eq.1) return
          call pvdata(ip1,x,ix,s,u,v,a,fldyn1,fldyn2,mstv,i,npm,npm1,
     +                numnp1(i),numel1(i),ityp,iswm)
          close(ip1)
        end do
      end if
c
c.... write .pvd file
      fname2 = fres
      call pvdfile(ip2,fname1,fname2,itime,n3,iswm)
      return
c
c.... write eigenvalue n1
3     fname3 = fres
      call vtueigv(ip1,fname3,n2,iret)
      if(iret.eq.1) return
      call pvdata_2(ip1,x,ix,u,numnp,numel,n1,n2)
      close (ip1)

c
      return
      end
c======================================================================
c
c
c
      subroutine vtufile(ip1,fname1,nummat1,itime,flgow,iret,iswm)
c-----------------------------------------------------------------------
c     defines name for .vtu file of each timestep
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer ip1,nummat1,itime
      character fname1*229,cnum1*3,cnum2*5,antw*3
      logical flgow,ex
      iret=0
c
c.... create .vtu file
      write(cnum1,'(i3)') nummat1
      if(cnum1(1:1).eq.' ') cnum1(1:1) = 'm'
      if(cnum1(2:2).eq.' ') cnum1(2:2) = '0'
      if(iswm.eq.0) call addpv(fname1,cnum1,3)
      write(cnum2,'(i5)') itime
      if(cnum2(1:1).eq.' ') cnum2(1:1) = 't'
      if(cnum2(2:2).eq.' ') cnum2(2:2) = '0'
      if(cnum2(3:3).eq.' ') cnum2(3:3) = '0'
      if(cnum2(4:4).eq.' ') cnum2(4:4) = '0'
      call addpv(fname1,cnum2,5)
      call addpv(fname1,'vtu',3)
      if(.not.flgow) then
        inquire(FILE=fname1,EXIST=ex)  ! Check if file already exists
        if(ex) then
          write(*,1001) 'File ',fname1,' exists! '
          write(*,1003) ; read(*,'(a)') antw ! overwrite y/n/all
          if(antw.eq.'y')then
            open(ip1,file=fname1,status='unknown',form='formatted')
            rewind(ip1)
          else if(antw.eq.'all')then
            open(ip1,file=fname1,status='unknown',form='formatted')
            rewind(ip1)
            flgow=.true.
          else
            iret=1
            return
          end if
        else     ! File does not exist
          open(ip1,file=fname1,status='unknown',form='formatted')
        end if
      else
        open(ip1,file=fname1,status='unknown',form='formatted')
      end if
c
      return
c
1001  format(a10,/,a229,/,a24)
1003  format('Overwrite [y,n,all]',' ? ',$)
c
      end
c
c
c
      subroutine vtueigv(ip1,fname3,ieigv,iret)
c-----------------------------------------------------------------------
c     defines name for .vtu file for (each) eigenvector
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      character cnum*4,fname3*229,antw*1
      integer ip1,ieigv,iret
      logical ex
      iret=0
c
      write(cnum,'(i4)') ieigv
      if(cnum(3:3).eq.' ') cnum(3:3) = '0'
      cnum(1:2) = 'EV'
      call addpv(fname3,cnum,4)
      call addpv(fname3,'vtu',3)
c
      inquire(FILE=fname3,EXIST=ex)                 ! Check if file already exists
      if(ex) then
        write(*,1001) 'File ',fname3,' exists! '
        write(*,1002) ; read(*,'(a)') antw          ! overwrite y/n
        if(antw.eq.'y') then
          open(ip1,file=fname3,status='unknown',form='formatted')
          rewind(ip1)
        else
          iret=1
          return
        endif
      else
        open(ip1,file=fname3,status='unknown',form='formatted')
      endif
c
      return
c
1001  format(a10,/,a229,/,a24)
1002  format('Overwrite [y,n]',' ? ',$)
c
      end
c
c
c
      subroutine pvdfile(ip2,fname1,fname2,itime,n3,iswm)
c-----------------------------------------------------------------------
c     defines name for .pvd file
c-----------------------------------------------------------------------
      USE cdata
      USE comfil
      implicit double precision (a-h,o-z)
      logical flgow,ex
      integer i,j,jj,ip2,itime,len1,len2,n3
      character fname1*229,fname2*229,fname3*229,fname4*229
      character cnum1*3,cnum2*5,antw*3
c
c
c.... create .pvd file
      call addpv(fname2,'pvd',3)
      inquire(FILE=fname2,EXIST=ex)
      if(ex) then
        if(n3.eq.1) then
          write(*,1001) 'File ',fname2,' exists! '
          write(*,1002) ; read(*,'(a)') antw ! overwrite y/n/all
          if(antw.eq.'y')then
            open(ip2,file=fname2,status='unknown',form='formatted')
             rewind(ip2)
           else
            return
          end if
        elseif(n3.eq.2) then
          open(ip2,file=fname2,status='unknown',form='formatted')
        endif
      else     ! File does not exist
        if(n3.eq.1) open(ip2,file=fname2,status='unknown',
     +                   form='formatted')
        if(n3.eq.2) call drawmess('PVD file needs to be initialized
     +                             [parv,init,n1]',1,0)
      end if
      write(ip2,1100)
c
      len1 = ipos1(fname1,229)
      len2 = ipos (fname1,229)
      jj = 0
      fname3 = ''
      do j = len1+1,len2,1
        jj = jj + 1
        write(fname3(jj:jj),'(a)') fres(j:j)
      end do
      do i = 1,itime,1
        nummat1 = nummat
        if(iswm.ne.0) nummat1 = 1
        do j = 1,nummat1,1
          fname4 = fname3
          write(cnum1,'(i3)') j
          if(cnum1(1:1).eq.' ') cnum1(1:1) = 'm'
          if(cnum1(2:2).eq.' ') cnum1(2:2) = '0'
          if(iswm.eq.0) call addpv(fname4,cnum1,3)
          write(cnum2,'(i5)') i
          if(cnum2(1:1).eq.' ') cnum2(1:1) = 't'
          if(cnum2(2:2).eq.' ') cnum2(2:2) = '0'
          if(cnum2(3:3).eq.' ') cnum2(3:3) = '0'
          if(cnum2(4:4).eq.' ') cnum2(4:4) = '0'
          call addpv(fname4,cnum2,5)
          call addpv(fname4,'vtu',3)
          fname4(jj+1:jj+3) = '"/>'
          write(ip2,1101) cnum2(2:5),cnum1(2:3),fname4
        end do
      end do
      write(ip2,1102)
      close(ip2)
c
      return
c
1001  format(a10,/,a229,/,a24)
1002  format('Overwrite [y,n]',' ? ',$)
c
1100  format('<?xml version="1.0"?>',/,/,
     +       '<VTKFile type="Collection" version="0.1"',
     +            ' byte_order="LittleEndian">',/,
     +       2x,'<Collection>')
1101  format(4x,'<DataSet timestep="',a,'" part="',a,'" file="',a)
1102  format(2x,'</Collection>',/,
     +       '</VTKFile>')
c
      end
c
c
c

      subroutine pvdata(ip1,x,ix,s,u,v,a,fldyn1,fldyn2,mstv,nummat1,
     +                  npm,npm1,numnp1,numel1,ityp,iswm)
c--------------------------------------------------------------------
c
c.... Purpose: store data fields for postprocessing in ParaView Format
c
c--------------------------------------------------------------------
      USE cdata
      USE pnodn
      USE sdata
      USE strnam
      implicit double precision(a-h,o-z)
      logical   fldyn1,fldyn2
      integer   ityp,vtktyp,numnp1,numel1,nummat1
      dimension x(ndm,numnp),ix(nen1,numel),s(numnp,*)
      dimension u(ndf,numnp),v(ndf,numnp),a(ndf,numnp)
      dimension npm(numnp),npm1(numnp,nummat),ixx(nen,numel)
      dimension n20(20)
c     dimension n27(27)
c     dimension n56(56)
c
c
c
      data xit/-999.d0/
c.... Node relation FEAP-Paraview for 20 node hexahedron
      data n20/1,2,3,4,5,6,7,8,13,14,15,16,18,19,20,21,9,10,11,12/

c.... Node relation FEAP-Paraview for 27 node hexahedron
c      data n27/1,2,3,4,5,6,7,8,13,14,15,16,18,19,20,21,9,10,11,12,
c     +23,24,25,26,17,22,27/ ! n20 ist enthalten!!

c.... Node relation FEAP-Paraview for 64 node hexahedron
c      data n56/ 1, 4,16,13, 2, 3, 8,12,15,14, 9, 5,
c     +         49,52,64,61,50,51,56,60,63,62,57,53,
c     +         17,20,32,29,18,19,24,28,31,30,25,21,
c     +         33,36,48,45,34,35,40,44,47,46,41,37,
c     +          6, 7,11,10,54,55,59,58/
c
c.... header
      write(ip1,2001) numnp1,numel1

c.... topology
        write(ip1,2020) 'Float32','position',3
        do i=1,numnp,1
          x1 = x(1,i)
          if(x1.eq.xit) call prxtie(gtie,i,x,x1,ndm) ! tied nodes
          x2 = x(2,i)
          x3 = 0.d0
          if(ndm.eq.3) x3=x(3,i)
          if(npm(i).eq.nummat1) write(ip1,2022) x1,x2,x3
        end do
        write(ip1,2030)
        write(ip1,2002)

c....   Node-element connection
        write(ip1,2020) 'Int32','connectivity',1
        do i=1,numel,1
          nel1 = 0.d0
          do k=1,nen,1
            if(ix(k,i).ne.0) nel1 = nel1 + 1
          end do
          nel2=nel1
          if(nel1.eq.27) nel2=20
c          if(nel1.eq.64) nel2=56
          if(ix(nen1,i).eq.nummat1) then
            do j=1,nel2,1
              ixx(j,i)=npm1(ix(j,i),nummat1)
              if(nel1.eq.20) ixx(j,i)=npm1(ix(n20(j),i),nummat1) ! 20 node brick
              if(nel1.eq.27) ixx(j,i)=npm1(ix(n20(j),i),nummat1) ! 27 node brick
c              if(nel1.eq.27) ixx(j,i)=npm1(ix(n27(j),i),nummat1) ! 27 node brick
c              if(nel1.eq.64) ixx(j,i)=npm1(ix(n56(j),i),nummat1) ! 64 node brick
            end do
            write(ip1,2021) (ixx(j,i)-1,j=1,nel2)
          end if
        end do
        write(ip1,2030)

c....   Node offsets
        write(ip1,2020) 'Int32','offsets',1
        noff = 0
        do i=1,numel,1
          nel1 = 0.d0
          do k=1,nen,1
            if(ix(k,i).ne.0) nel1 = nel1 + 1
          end do
          nel2=nel1
          if(nel1.eq.27) nel2=20
c          if(nel1.eq.64) nel2=56
          if(ix(nen1,i).eq.nummat1) then
            noff = noff + nel2
            write(ip1,2021) noff
          end if
        end do
        write(ip1,2030)

c....   VTKtypes
        write(ip1,2020) 'UInt8','types',1
        do i=1,numel,1
c....     vtktyp of each element: vtktyp(i)
          nel1 = 0.d0
          do k=1,nen,1
            if(ix(k,i).ne.0) nel1 = nel1 + 1
          end do
          if(nel1.eq.2) vtktyp=3                   ! 2 node element
          if(nel1.eq.3) vtktyp=5                   ! 3 node triangle
          if(nel1.eq.4.and.ityp.ne.16) vtktyp=9    ! 4 node quad
          if(nel1.eq.4.and.ityp.eq.16) vtktyp=10   ! 4 node tetrahedron
          if(nel1.eq.8)then
            vtktyp=23                              ! 8 node quad
            if(ityp.eq.5.or.ityp.eq.6.or.ityp.eq.7.or.ityp.eq.15)
     +         vtktyp=12                           ! 8 node hexahedron
          end if
          if(nel1.eq.9)  vtktyp=28                 !  9 node quad
          if(nel1.eq.20) vtktyp=25                 ! 20 node hexahedron
          if(nel1.eq.27) vtktyp=25                 ! 27 node hexahedron
c          if(nel1.eq.27) vtktyp=29                ! 27 node hexahedron
c          if(nel1.eq.64) vtktyp=67                ! 64 node hexahedron
c
          if(ix(nen1,i).eq.nummat1) write(ip1,2021) vtktyp
        end do
        write(ip1,2030)
      write(ip1,2003)
c
c.... write Pointdata
      write(ip1,2004)

c....   Displacements
        write(ip1,2020) 'Float32','Displacements',ndf
        do i=1,numnp,1
          if(npm(i).eq.nummat1) write(ip1,2022) (u(j,i),j=1,ndf)
        end do
        write(ip1,2030)
c
c....   Deformed mesh = u(1-3)
        write(ip1,2020) 'Float32','Deformed mesh',3
        do i=1,numnp,1
          x1 = x(1,i)
          if(x1.eq.xit) call prxtie(gtie,i,x,x1,ndm) ! tied nodes
          x2 = x(2,i)
          x3 = 0.d0
          if(ndm.eq.3) x3=x(3,i)
          u1=u(1,i)
          u2=u(2,i)
          u3=0.d0
          if(ndm.eq.3) u3=u(3,i)
          if(ityp.eq.3) then  ! plate
            u1=0.d0
            u2=0.d0
            u3=u(1,i)
          end if
          if(npm(i).eq.nummat1) write(ip1,2022) u1,u2,u3
        end do
        write(ip1,2030)
c
        if(fldyn1) then
c....     Velocity
          write(ip1,2020) 'Float32','Velocity',ndf
          do i=1,numnp,1
            if(npm(i).eq.nummat1) write(ip1,2022) (v(j,i),j=1,ndf)
          end do
          write(ip1,2030)

          if(fldyn2) then
c....       Acceleration
            write(ip1,2020) 'Float32','Acceleration',ndf
            do i=1,numnp,1
              if(npm(i).eq.nummat1) write(ip1,2022) (a(j,i),j=1,ndf)
            end do
            write(ip1,2030)
          end if
        end if
c
c....   Stresses
        write(ip1,2020) 'Float32','Stress',mstv
        do i=1,numnp,1
          if(npm(i).eq.nummat1) write(ip1,2022) (s(i,j),j=1,mstv)
        end do
        write(ip1,2030)
      write(ip1,2005)
c
c.... write CellData
      write(ip1,2006)
c....   Material sets
        write(ip1,2020) 'Int32','Material',1
        do i=1,numel,1
          if(ix(nen1,i).eq.nummat1) write(ip1,2021) ix(nen1,i)
        end do
        write(ip1,2030)
      write(ip1,2007)
c
c.... write end
      write(ip1,2008)
c
2001  format('<?xml version="1.0"?>',/,/,
     +       '<VTKFile type="UnstructuredGrid" version="0.1"',
     +            ' byte_order="LittleEndian">',/,
     +       2x,'<UnstructuredGrid>',/,
     +       4x,'<Piece NumberOfPoints="',i10,'" NumberOfCells="',i10,
     +            '">',/,
     +       6x,'<Points>')
2002  format(6x,'</Points>',/,
     +       6x,'<Cells>')
2003  format(6x,'</Cells>')
2004  format(6x,'<PointData>')
2005  format(6x,'</PointData>')
2006  format(6x,'<CellData>')
2007  format(6x,'</CellData>')
2008  format(4x,'</Piece>',/,
     +       2x,'</UnstructuredGrid>',/,
     +       '</VTKFile>')

2020  format(8x,'<DataArray type="',a,'" Name="',a,'"',
     +          ' NumberOfComponents="',i2,'" format="ascii">')
2021  format(10x,20(i10,2x))
2022  format(10x,25(e13.7,2x))
2030  format(8x,'</DataArray>')
c
      return
      end
c
c
c
      subroutine pvdata_1(ip1,x,ix,s,u,v,a,fldyn1,fldyn2,mstv,nummat1,
     +                  numnp1,numel1,ityp)
c--------------------------------------------------------------------
c
c.... Purpose: store data fields for postprocessing in ParaView Format
c
c--------------------------------------------------------------------
      USE cdata
      USE pnodn
      USE sdata
      USE strnam
      implicit double precision(a-h,o-z)
      logical   fldyn1,fldyn2
      integer   ityp,vtktyp,numnp1,numel1,nummat1
      dimension x(ndm,numnp),ix(nen1,numel),s(numnp,*)
      dimension u(ndf,numnp),v(ndf,numnp),a(ndf,numnp)
      dimension ixx(nen,numel)
      dimension n20(20)
c     dimension n27(27)
c     dimension n56(56)
c
c
c
      data xit/-999.d0/
c.... Node relation FEAP-Paraview for 20 node hexahedron
      data n20/1,2,3,4,5,6,7,8,13,14,15,16,18,19,20,21,9,10,11,12/

c.... Node relation FEAP-Paraview for 27 node hexahedron
c      data n27/1,2,3,4,5,6,7,8,13,14,15,16,18,19,20,21,9,10,11,12,
c     +23,24,25,26,17,22,27/ ! n20 ist enthalten!!

c.... Node relation FEAP-Paraview for 64 node hexahedron
c      data n56/ 1, 4,16,13, 2, 3, 8,12,15,14, 9, 5,
c     +         49,52,64,61,50,51,56,60,63,62,57,53,
c     +         17,20,32,29,18,19,24,28,31,30,25,21,
c     +         33,36,48,45,34,35,40,44,47,46,41,37,
c     +          6, 7,11,10,54,55,59,58/
c
c.... header
      write(ip1,2001) numnp1,numel1

c.... topology
        write(ip1,2020) 'Float32','position',3
        do i=1,numnp,1
          x1 = x(1,i)
          if(x1.eq.xit) call prxtie(gtie,i,x,x1,ndm) ! tied nodes
          x2 = x(2,i)
          x3 = 0.d0
          if(ndm.eq.3) x3=x(3,i)
          write(ip1,2022) x1,x2,x3
        end do
        write(ip1,2030)
        write(ip1,2002)

c....   Node-element connection
        write(ip1,2020) 'Int32','connectivity',1
        do i=1,numel,1
          nel1 = 0.d0
          do k=1,nen,1
            if(ix(k,i).ne.0) nel1 = nel1 + 1
          end do
          nel2=nel1
          if(nel1.eq.27) nel2=20
c          if(nel1.eq.64) nel2=56
          write(ip1,2021) (ix(j,i)-1,j=1,nel2)
        end do
        write(ip1,2030)

c....   Node offsets
        write(ip1,2020) 'Int32','offsets',1
        noff = 0
        do i=1,numel,1
          nel1 = 0.d0
          do k=1,nen,1
            if(ix(k,i).ne.0) nel1 = nel1 + 1
          end do
          nel2=nel1
          if(nel1.eq.27) nel2=20
c          if(nel1.eq.64) nel2=56
          noff = noff + nel2
          write(ip1,2021) noff
        end do
        write(ip1,2030)

c....   VTKtypes
        write(ip1,2020) 'UInt8','types',1
        do i=1,numel,1
c....     vtktyp of each element: vtktyp(i)
          nel1 = 0.d0
          do k=1,nen,1
            if(ix(k,i).ne.0) nel1 = nel1 + 1
          end do
          if(nel1.eq.2) vtktyp=3                   ! 2 node element
          if(nel1.eq.3) vtktyp=5                   ! 3 node triangle
          if(nel1.eq.4.and.ityp.ne.15) vtktyp=9    ! 4 node quad
          if(nel1.eq.4.and.ityp.eq.15) vtktyp=10   ! 4 node tetrahedron
          if(nel1.eq.8)then
            vtktyp=23                              ! 8 node quad
            if(ityp.eq.5.or.ityp.eq.6.or.ityp.eq.7.or.ityp.eq.14)
     +         vtktyp=12                           ! 8 node hexahedron
          end if
          if(nel1.eq.9)  vtktyp=28                 !  9 node quad
          if(nel1.eq.20) vtktyp=25                 ! 20 node hexahedron
          if(nel1.eq.27) vtktyp=25                 ! 27 node hexahedron
c          if(nel1.eq.27) vtktyp=29                ! 27 node hexahedron
c          if(nel1.eq.64) vtktyp=67                ! 64 node hexahedron
c
          write(ip1,2021) vtktyp
        end do
        write(ip1,2030)
      write(ip1,2003)
c
c.... write Pointdata
      write(ip1,2004)

c....   Displacements
        write(ip1,2020) 'Float32','Displacements',ndf
        do i=1,numnp,1
          write(ip1,2022) (u(j,i),j=1,ndf)
        end do
        write(ip1,2030)
c
c....   Deformed mesh = u(1-3)
        write(ip1,2020) 'Float32','Deformed mesh',3
        do i=1,numnp,1
          x1 = x(1,i)
          if(x1.eq.xit) call prxtie(gtie,i,x,x1,ndm) ! tied nodes
          x2 = x(2,i)
          x3 = 0.d0
          if(ndm.eq.3) x3=x(3,i)
          u1=u(1,i)
          u2=u(2,i)
          u3=0.d0
          if(ndm.eq.3) u3=u(3,i)
          if(ityp.eq.3) then  ! plate
            u1=0.d0
            u2=0.d0
            u3=u(1,i)
          end if
          write(ip1,2022) u1,u2,u3
        end do
        write(ip1,2030)
c
        if(fldyn1) then
c....     Velocity
          write(ip1,2020) 'Float32','Velocity',ndf
          do i=1,numnp,1
            write(ip1,2022) (v(j,i),j=1,ndf)
          end do
          write(ip1,2030)

          if(fldyn2) then
c....       Acceleration
            write(ip1,2020) 'Float32','Acceleration',ndf
            do i=1,numnp,1
              write(ip1,2022) (a(j,i),j=1,ndf)
            end do
            write(ip1,2030)
          end if
        end if
c
c....   Stresses
        write(ip1,2020) 'Float32','Stress',mstv
        do i=1,numnp,1
          write(ip1,2022) (s(i,j),j=1,mstv)
        end do
        write(ip1,2030)
      write(ip1,2005)
c
c.... write CellData
      write(ip1,2006)
c....   Material sets
        write(ip1,2020) 'Int32','Material',1
        do i=1,numel,1
          write(ip1,2021) ix(nen1,i)
        end do
        write(ip1,2030)
      write(ip1,2007)
c
c.... write end
      write(ip1,2008)
c
2001  format('<?xml version="1.0"?>',/,/,
     +       '<VTKFile type="UnstructuredGrid" version="0.1"',
     +            ' byte_order="LittleEndian">',/,
     +       2x,'<UnstructuredGrid>',/,
     +       4x,'<Piece NumberOfPoints="',i10,'" NumberOfCells="',i10,
     +            '">',/,
     +       6x,'<Points>')
2002  format(6x,'</Points>',/,
     +       6x,'<Cells>')
2003  format(6x,'</Cells>')
2004  format(6x,'<PointData>')
2005  format(6x,'</PointData>')
2006  format(6x,'<CellData>')
2007  format(6x,'</CellData>')
2008  format(4x,'</Piece>',/,
     +       2x,'</UnstructuredGrid>',/,
     +       '</VTKFile>')

2020  format(8x,'<DataArray type="',a,'" Name="',a,'"',
     +          ' NumberOfComponents="',i2,'" format="ascii">')
2021  format(10x,20(i10,2x))
2022  format(10x,25(e13.7,2x))
2030  format(8x,'</DataArray>')
c
      return
      end
c
c
c
      subroutine pvdata_2(ip1,x,ix,u,numnp1,numel1,ityp,n2)
c--------------------------------------------------------------------
c
c.... Purpose: store data fields for postprocessing in ParaView Format
c
c--------------------------------------------------------------------
      USE cdata
      USE pnodn
      USE sdata
      USE strnam
      implicit double precision(a-h,o-z)
      logical   fldyn1,fldyn2
      integer   ityp,vtktyp,numnp1,numel1,nummat1
      dimension x(ndm,numnp),ix(nen1,numel),u(ndf,numnp)
      dimension ixx(nen,numel)
      dimension n20(20)
c     dimension n27(27)
c     dimension n56(56)
c
c
c
      data xit/-999.d0/
c.... Node relation FEAP-Paraview for 20 node hexahedron
      data n20/1,2,3,4,5,6,7,8,13,14,15,16,18,19,20,21,9,10,11,12/

c.... Node relation FEAP-Paraview for 27 node hexahedron
c      data n27/1,2,3,4,5,6,7,8,13,14,15,16,18,19,20,21,9,10,11,12,
c     +23,24,25,26,17,22,27/ ! n20 ist enthalten!!

c.... Node relation FEAP-Paraview for 64 node hexahedron
c      data n56/ 1, 4,16,13, 2, 3, 8,12,15,14, 9, 5,
c     +         49,52,64,61,50,51,56,60,63,62,57,53,
c     +         17,20,32,29,18,19,24,28,31,30,25,21,
c     +         33,36,48,45,34,35,40,44,47,46,41,37,
c     +          6, 7,11,10,54,55,59,58/
c
c.... header
      write(ip1,2001) numnp1,numel1

c.... topology
        write(ip1,2020) 'Float32','position',3
        do i=1,numnp,1
          x1 = x(1,i)
          if(x1.eq.xit) call prxtie(gtie,i,x,x1,ndm) ! tied nodes
          x2 = x(2,i)
          x3 = 0.d0
          if(ndm.eq.3) x3=x(3,i)
          write(ip1,2022) x1,x2,x3
        end do
        write(ip1,2030)
        write(ip1,2002)

c....   Node-element connection
        write(ip1,2020) 'Int32','connectivity',1
        do i=1,numel,1
          nel1 = 0.d0
          do k=1,nen,1
            if(ix(k,i).ne.0) nel1 = nel1 + 1
          end do
          nel2=nel1
          if(nel1.eq.27) nel2=20
c          if(nel1.eq.64) nel2=56
          write(ip1,2021) (ix(j,i)-1,j=1,nel2)
        end do
        write(ip1,2030)

c....   Node offsets
        write(ip1,2020) 'Int32','offsets',1
        noff = 0
        do i=1,numel,1
          nel1 = 0.d0
          do k=1,nen,1
            if(ix(k,i).ne.0) nel1 = nel1 + 1
          end do
          nel2=nel1
          if(nel1.eq.27) nel2=20
c          if(nel1.eq.64) nel2=56
          noff = noff + nel2
          write(ip1,2021) noff
        end do
        write(ip1,2030)

c....   VTKtypes
        write(ip1,2020) 'UInt8','types',1
        do i=1,numel,1
c....     vtktyp of each element: vtktyp(i)
          nel1 = 0.d0
          do k=1,nen,1
            if(ix(k,i).ne.0) nel1 = nel1 + 1
          end do
          if(nel1.eq.2) vtktyp=3                   ! 2 node element
          if(nel1.eq.3) vtktyp=5                   ! 3 node triangle
          if(nel1.eq.4.and.ityp.ne.15) vtktyp=9    ! 4 node quad
          if(nel1.eq.4.and.ityp.eq.15) vtktyp=10   ! 4 node tetrahedron
          if(nel1.eq.8)then
            vtktyp=23                              ! 8 node quad
            if(ityp.eq.5.or.ityp.eq.6.or.ityp.eq.7.or.ityp.eq.14)
     +         vtktyp=12                           ! 8 node hexahedron
          end if
          if(nel1.eq.9)  vtktyp=28                 !  9 node quad
          if(nel1.eq.20) vtktyp=25                 ! 20 node hexahedron
          if(nel1.eq.27) vtktyp=25                 ! 27 node hexahedron
c          if(nel1.eq.27) vtktyp=29                ! 27 node hexahedron
c          if(nel1.eq.64) vtktyp=67                ! 64 node hexahedron
c
          write(ip1,2021) vtktyp
        end do
        write(ip1,2030)
      write(ip1,2003)
c
c.... write Pointdata
      write(ip1,2004)

c....   Eigenvector
        write(ip1,2023) 'Float32',n2,'. Eigenvector',ndf
        do i=1,numnp,1
          write(ip1,2022) (u(j,i),j=1,ndf)
        end do
        write(ip1,2030)
c
c....   Deformed mesh = u(1-3)
        write(ip1,2020) 'Float32','Deformed mesh',3
        do i=1,numnp,1
          x1 = x(1,i)
          if(x1.eq.xit) call prxtie(gtie,i,x,x1,ndm) ! tied nodes
          x2 = x(2,i)
          x3 = 0.d0
          if(ndm.eq.3) x3=x(3,i)
          u1=u(1,i)
          u2=u(2,i)
          u3=0.d0
          if(ndm.eq.3) u3=u(3,i)
          if(ityp.eq.3) then  ! plate
            u1=0.d0
            u2=0.d0
            u3=u(1,i)
          end if
          write(ip1,2022) u1,u2,u3
        end do
        write(ip1,2030)
      write(ip1,2005)
c
c.... write end
      write(ip1,2008)
c
2001  format('<?xml version="1.0"?>',/,/,
     +       '<VTKFile type="UnstructuredGrid" version="0.1"',
     +            ' byte_order="LittleEndian">',/,
     +       2x,'<UnstructuredGrid>',/,
     +       4x,'<Piece NumberOfPoints="',i10,'" NumberOfCells="',i10,
     +            '">',/,
     +       6x,'<Points>')
2002  format(6x,'</Points>',/,
     +       6x,'<Cells>')
2003  format(6x,'</Cells>')
2004  format(6x,'<PointData>')
2005  format(6x,'</PointData>')
2008  format(4x,'</Piece>',/,
     +       2x,'</UnstructuredGrid>',/,
     +       '</VTKFile>')

2020  format(8x,'<DataArray type="',a,'" Name="',a,'"',
     +          ' NumberOfComponents="',i2,'" format="ascii">')
2021  format(10x,20(i10,2x))
2022  format(10x,25(e13.7,2x))
2023  format(8x,'<DataArray type="',a,'" Name="'i2,a,'"',
     +          ' NumberOfComponents="',i2,'" format="ascii">')
2030  format(8x,'</DataArray>')
c
      return
      end
c
c
c
      subroutine addpv(fnam,fext,k)
c-----------------------------------------------------------------------
c      Purpose: adds character string to file fnam->fnam.fext
c      Input:
c         fnam(229)  -  character string without extension fnam
c         fext(k)    -  extension to add
c      Output:
c         fnam(229)  -  character string with extension    fnam.fext
c--------------------------------------------------------------------------
      character fnam*229
      character*1 fnam1(229),fext1(k),fext(k)
      do i = 1,229
        fnam1(i) = fnam(i:i)
      end do
      do i = 1,k
        fext1(i) = fext(i)
      end do
      iposl = ipos(fnam1,229)
      iposx = ipos(fext1,k)
      ii = iposl + 1
      do i = ii,229
          fnam1(i) = ' '
      end do
      fnam1(ii) = '.'
      if((ii+iposx).gt.229) then
        call drawmess('Filename + Extension is to long(<229!)',1,0)
        return
      end if
      do i = 1,iposx
          fnam1(ii+i) = fext1(i)
      end do
      do i = 1,229
        fnam(i:i) = fnam1(i)
      end do
      return
      end