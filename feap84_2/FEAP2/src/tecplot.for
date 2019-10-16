      subroutine tecpost(x,ix,s,u,nen1,ndm,ndf,ktec,ityp1)
c----------------------------------------------------------------------
c
c.... Purpose: store data fields for postprocessing in TecPlot Format
c
c     Inputs:                                                          |
c     x(ndm,*)   = nodal cordinates
c     ix(nen1,*) = node-element-relation
c     s(numnp,*) = nodal stresses (not for beam elements
c     u(ndf,*)   = nodal displacements
c     nen1       = nen+4 nen=number of nodes/element
c     ndm        = dimension of problem
c     ndf        = dof in problem
c     ktec       = n1
c     ityp1      = n2 type of element, see below array eltyp
c
c     [tec,init,n2]     initialize outputfile rfile.tec,  n1=1 n2=ityp define element
c                       write nodal coordinates+element-node-relations
c
c     [tec,write]       write all values (def)            n1=3
c     [tec,eigv,n1,n2]  write eigenvector n1 instead of displacements, scaling factor n2
c
c     [tec,close]       close file                        n1=2
c
c     Outputs
c     rfile.tec
c
c     Tecplot support triangle, quadrilateral, brick
c
c     modifications:
c     write only in first step
c     1)   in col 4-6 data for x,y,z
c     1)   in col 3-4 data for x,y
c     2)   connectivity data
c
c----------------------------------------------------------------------
      USE cdata
      USE comfil
      USE iofile
      USE pnodn
      USE strnam
      USE tdata
      implicit double precision (a-h,o-z)
      logical ex, first
      character*1 antw,eltyp*69,eltyp1*69,text*111
      character*229 fname1
      character     strsust*15
      dimension eltyp(16)  ! Element-Typ
      dimension x(ndm,*)   ! Knoten Koordinaten
      dimension ix(nen1,*) ! Element-Knoten-Beziehung
      dimension s(numnp,*) ! Knoten Spannungen
      dimension u(ndf,*)   ! Knoten Verschiebungen
      dimension strsust(25)! Name der Knotenspannung

      save numel1,first,ityp

      data ip1/23/
      data ntyp/16/
      data eltyp/
     1'plane stress ndf = 2  ndm = 2,3 nen = 3,4,9                   ',
     2'beam2d/axish ndf = 3  ndm = 2   nen = 2                       ',
     3'plate        ndf = 3  ndm = 2,3 nen = 3,4,9                   ',
     4'plain strain ndf = 3  ndm = 2   nen = 4,9   2D cosserat       ',
     5'brick        ndf = 3  ndm =   3 nen = 8,27                    ',
     6'brick        ndf = 4  ndm =   3 nen = 8,27  th.-mech. coupling',
     7'brick        ndf = 6  ndm =   3 nen = 8,27  3d cosserat       ',
     8'shells5      ndf = 5  ndm =   3 nen = 3,4,9                   ',
     9'beam3d       ndf = 6  ndm =   3 nen = 2                       ',
     +'shells6      ndf = 6  ndm =   3 nen = 3,4,9                   ',
     1'beam3d       ndf = 7  ndm =   3 nen = 2     incl. warping     ',
     2'shell        ndf = 7  ndm =   3 nen = 4,9   incl. warping beam',
     3'plain strain ndf = 1  ndm =   2 nen = 4,9,16 2D phase field   ',
     4'plain strain ndf = 5  ndm =   2 nen = 4,9,16 2D phase field   ',
     5'brick        ndf = 7  ndm =   3 nen = 8,27 3D phase field     ',
     6'tetraeder    ndf = 3  ndm =   3 nen = 4                       '/

      data xit/-999.0d0/
      mstv = iabs(istv) ! plot only mstv stress values max=24!

      goto (1,2,3) ktec
c
1     first=.true.
      ityp = ityp1

      if(ityp.lt.1.or.ityp.gt.ntyp) then
        write(*,'(a23)') 'Element not implemented'
        write(*,'(a)')   'available elements'
        do i = 1,ntyp
          write(*,'(i3,a,a)')   i,' ',eltyp(i)
        end do
        return
      end if

c...  open file for TecPlot Data
      fname1 = fres
      call addext (fname1,'tec ')
      inquire(FILE=fname1,EXIST=ex)   ! Pruefen ob Datei schon existiert
      if(ex) then                    ! => Datei existiert
        write (*,1001) 'File ', fname1, ' exists! '
1001    format(a10,/,a229,/,a24)
        write (*,1002)
1002    format('Overwrite [y,n]',' ? ',$)
        read (*,'(a)') antw ! overwrite y/n
        if(antw.eq.'y') then
          open(unit=ip1,file=fname1,status='unknown',form='formatted')
          rewind(ip1)
        else
          return
        end if
      else     ! => Datei existiert nicht
          open(unit=ip1,file=fname1,status='unknown',form='formatted')
      end if
      text='Files for Postprocessing opened, Element: '
      eltyp1=eltyp(ityp)
      text(42+1:42+69)=eltyp1(1:69)
      call drawmess(text,1,-3)
      return
c...  close file for TecPlot Data
2     close(ip1)
      call drawmess('Files for Postprocessing closed',1,-3)
      return
c
c.... write all data for step
3     continue

c...  count number of active plot elements
      numel1 = 0
      do  80 i=1, numel
        ma = ix(nen1,i)
        if(iplma(ma).eq.0) goto 80
        numel1 = numel1 + 1
80    continue
      if(first) then
c.....  Tec_Plot Header rausschreiben
c....   name of stress
        do i = 1,25
          if(strsus(i).eq.'               ') then
            write(strsust(i),'(a9,i3)' ) '  STRESS ',i
          else
            strsust(i) = strsus(i)
          end if
        end do
        write (ip1,'(a14)',err=99) 'TITLE = "FEAP"'

        if(ityp.eq.1) then  ! 3/4 node plane stress element 2d/3d
          write(ip1,'(1(a50),5(a9),24(a19))')
     +            'VARIABLES = "X + U(x)", "Y + U(y)", "Z+0",     ',
     +            '"X",     ',
     +            '"Y",     ',
     +            '"Z",     ',
     +            '"disp x",',
     +            '"disp y" ',
     +            (', '//'"'//strsust(i)//'"',i=1,mstv)

        else if(ityp.eq.2) then  !  2 node beam element 2d/axishell
          write(ip1,'(1(a34),5(a9))')
     +            'VARIABLES = "X + U(x)","Y + U(y)",',
     +            '"X",     ',
     +            '"Y",     ',
     +            '"Disp x",',
     +            '"Disp y",',
     +            '"Rot z"  '

        else if(ityp.eq.3) then  !  3/4 node plate element 2d/3d
          write(ip1,'(1(a33),6(a9),24(a19))')
     +            'VARIABLES = "X+0", "Y+0", "Z + U(z)",',
     +            '"X",     ',
     +            '"Y",     ',
     +            '"Z",     ',
     +            '"Disp z",',
     +            '"Rot x", ',
     +            '"Rot y"  ',
     +            (', '//'"'//strsust(i)//'"',i=1,mstv)

        elseif(ityp.eq.4) then  !  cosserat element 2d
          write(ip1,'(1(a47),6(a9),24(a19))')
     +            'VARIABLES = "X + U(x)", "Y + U(y)", "Z + 0   ",',
     +            '"X",     ',
     +            '"Y",     ',
     +            '"Z",     ',
     +            '"Disp x",',
     +            '"Disp y", ',
     +            '"Rot  z"  ',
     +            (', '//'"'//strsust(i)//'"',i=1,mstv)

        else if(ityp.eq.5) then  ! 8/27 node hexahedral element
          write(ip1,'(1(a47),6(a9),24(a19))')
     +            'VARIABLES = "X + U(x)", "Y + U(y)", "Z + U(z)",',
     +            '"X",     ',
     +            '"Y",     ',
     +            '"Z",     ',
     +            '"disp x",',
     +            '"disp y",',
     +            '"disp z" ',
     +            (', '//'"'//strsust(i)//'"',i=1,mstv)

        else if(ityp.eq.6) then  ! 3D-brick element for electromechanical/thermomechanical coupling
          write(ip1,'(1(a47),7(a9),24(a19))')
     +           'VARIABLES="X+U(x)", "Y+U(y)", "Z+U(z)   ",',
     +           '"X",      ',
     +           '"Y",      ',
     +           '"Z",      ',
     +           '"disp x", ',
     +           '"disp y", ',
     +           '"disp z", ',
     +           '"d Phi"   ',
     +           (', '//'"'//strsust(i)//'"',i=1,mstv)

        else if(ityp.eq.7) then  ! 3D-brick element for cosserat theory
          write(ip1,'(1(a47),9(a9),24(a19))')
     +           'VARIABLES="X+U(x)", "Y+U(y)", "Z+U(z)   ",',
     +           '"X",      ',
     +           '"Y",      ',
     +           '"Z",      ',
     +           '"disp x", ',
     +           '"disp y", ',
     +           '"disp z", ',
     +           '"rot x",  ',
     +           '"rot y",  ',
     +           '"rot z",  ',
     +           (', '//'"'//strsust(i)//'"',i=1,mstv)

        else if(ityp.eq.8) then  ! 3/4 node 5-parameter shell element
          write(ip1,'(1(a47),8(a9),24(a19))')
     +            'VARIABLES = "X + U(x)", "Y + U(y)", "Z + U(z)",',
     +            '"X",     ',
     +            '"Y",     ',
     +            '"Z",     ',
     +            '"disp x",',
     +            '"disp y",',
     +            '"disp z",',
     +            '"rot x", ',
     +            '"rot y"  ',
     +            (', '//'"'//strsust(i)//'"',i=1,mstv)

        else if(ityp.eq.9) then  ! 2 node 3D-beam element
          write(ip1,'(1(a47),9(a9))')
     +            'VARIABLES = "X + U(x)", "Y + U(y)", "Z + U(z)",',
     +            '"X",     ',
     +            '"Y",     ',
     +            '"Z",     ',
     +            '"disp x",',
     +            '"disp y",',
     +            '"disp z",',
     +            '"rot x", ',
     +            '"rot y"  ',
     +            '"rot z"  '

        else if(ityp.eq.10) then  ! 3/4 node 6-parameter shell element
          write(ip1,'(1(a47),9(a9),24(a19))')
     +            'VARIABLES = "X + U(x)", "Y + U(y)", "Z + U(z)",',
     +            '"X",     ',
     +            '"Y",     ',
     +            '"Z",     ',
     +            '"disp x",',
     +            '"disp y",',
     +            '"disp z",',
     +            '"rot x", ',
     +            '"rot y"  ',
     +            '"rot z"  ',
     +            (', '//'"'//strsust(i)//'"',i=1,mstv)

        else if(ityp.eq.11) then  ! 2 node 3D-beam element incl. warping
          write(ip1,'(1(a47),10(a9))')
     +            'VARIABLES = "X + U(x)", "Y + U(y)", "Z + U(z)",',
     +            '"X",     ',
     +            '"Y",     ',
     +            '"Z",     ',
     +            '"disp x",',
     +            '"disp y",',
     +            '"disp z",',
     +            '"rot x", ',
     +            '"rot y"  ',
     +            '"rot z"  ',
     +            '"warp "  '

        else if(ityp.eq.12) then  ! 4 node 3D-shell element incl. warping dof of beam3D7dof
          write(ip1,'(1(a47),10(a9),24(a19))')
     +            'VARIABLES = "X + U(x)", "Y + U(y)", "Z + U(z)",',
     +            '"X",     ',
     +            '"Y",     ',
     +            '"Z",     ',
     +            '"disp x",',
     +            '"disp y",',
     +            '"disp z",',
     +            '"rot x", ',
     +            '"rot y"  ',
     +            '"rot z"  ',
     +            '"warp "  ',
     +            (', '//'"'//strsust(i)//'"',i=1,mstv)

        elseif(ityp.eq.13) then  !plain element, phase field theory, ndf=1
          write(ip1,'(1(a47),8(a9),24(a19))') 
     +            'VARIABLES = "X + U(x)", "Y + U(y)", "Z + 0   ",', 
     +            '"X",     ',
     +            '"Y",     ',
     +            '"Z",     ',
     +            '"Phi   "  ',
     +            (', '//'"'//strsust(i)//'"',i=1,mstv)
        elseif(ityp.eq.14) then  !plain element, phase field theory, ndf=5
          write(ip1,'(1(a47),8(a9),24(a19))') 
     +            'VARIABLES = "X + U(x)", "Y + U(y)", "Z + 0   ",', 
     +            '"X",     ',
     +            '"Y",     ',
     +            '"Z",     ',
     +            '"disp x",',
     +            '"disp y", ',
     +            '"Pol  x"  ',
     +            '"Pol  y"  ',
     +            '"Phi   "  ',
     +            (', '//'"'//strsust(i)//'"',i=1,mstv)
        else if(ityp.eq.15) then  !brick element, phase field theory, ndf=7
          write(ip1,'(1(a47),10(a9),24(a19))')
     +           'VARIABLES="X+U(x)", "Y+U(y)", "Z+U(z)   ",', 
     +           '"X",      ',
     +           '"Y",      ',
     +           '"Z",      ',
     +           '"disp x", ',
     +           '"disp y", ',
     +           '"disp z", ',
     +           '"Pol x",  ',
     +           '"Pol y",  ',
     +           '"Pol z",  ',
     +           '"Phi  ",  ',
     +           (', '//'"'//strsust(i)//'"',i=1,mstv)
        else if(ityp.eq.16) then  ! 4 node tetrahedron element
          write(ip1,'(1(a47),6(a9),24(a19))') 
     +            'VARIABLES = "X + U(x)", "Y + U(y)", "Z + U(z)",', 
     +            '"X",     ',
     +            '"Y",     ',
     +            '"Z",     ',
     +            '"disp x",',
     +            '"disp y",',
     +            '"disp z" ',
     +            (', '//'"'//strsust(i)//'"',i=1,mstv)
        else
          write(*,'(a23)') 'Element not implemented'
          return
        end if

c....   Zone Header
        if    (nen.eq.2) then                   ! 2 node element (Beam2D/3D/Axishell)
          write(ip1,'(a18,i8,a4,i8,a25,a15,e12.5)',err=99)
     +         'ZONE T="FEAP", N=',numnp,', E=',numel1,
     +         ', F=FEPOINT, ET=TRIANGLE,',', SOLUTIONTIME=',ttim

        else if(nen.eq.3) then                   ! 3 node element (Beam, Plate etc)
          write(ip1,'(a18,i8,a4,i8,a25,a15,e12.5)',err=99)
     +         'ZONE T="FEAP", N=',numnp,', E=', numel1,
     +         ', F=FEPOINT, ET=TRIANGLE',', SOLUTIONTIME=',ttim

        else if(nen.eq.4.and.ityp.ne.15) then               ! Quadrilateral-Element
          write(ip1,'(a18,i8,a4,i8,a30,a15,e12.5)',err=99)
     +    'ZONE T="FEAP", N=',numnp,', E=',numel1,
     +    ', F=FEPOINT, ET=QUADRILATERAL',', SOLUTIONTIME=',ttim

        else if(nen.eq.4.and.ityp.eq.15) then ! Tetrahedron-Element
          write(ip1,'(a18,i8,a4,i8,a30,a15,e12.5)',err=99)
     +    'ZONE T="FEAP", N=',numnp,', E=',numel1,
     +    ', F=FEPOINT, ET=TETRAHEDRON',', SOLUTIONTIME=',ttim

        else if(nen.eq.8) then               ! hexahedral element
          write(ip1,'(a18,i8,a4,i8,a22,a15,e12.5)',err=99)
     +    'ZONE T="FEAP", N=',numnp,', E=',numel1,
     +    ', F=FEPOINT, ET=BRICK',', SOLUTIONTIME=',ttim

        else if(nen.eq.9) then               ! 4 quadrilateral elements (from 9 node element)
          write(ip1,'(a18,i8,a4,i8,a30,a15,e12.5)',err=99)
     +    'ZONE T="FEAP", N=',numnp,', E=',4*numel1,
     +    ', F=FEPOINT, ET=QUADRILATERAL',', SOLUTIONTIME=',ttim

        else if(nen.eq.16) then               ! 9 quadrilateral elements (from 16 node element)
          write(ip1,'(a18,i8,a4,i8,a30,a15,e12.5)',err=99)
     +    'ZONE T="FEAP", N=',numnp,', E=',9*numel1,
     +    ', F=FEPOINT, ET=QUADRILATERAL',', SOLUTIONTIME=',ttim
       else if(nen.eq.27) then               ! 8 hexahedral elements (from 27 node element)
          write(ip1,'(a18,i8,a4,i8,a22,a15,e12.5)',err=99)
     +    'ZONE T="FEAP", N=',numnp,', E=',8*numel1,
     +    ', F=FEPOINT, ET=BRICK',', SOLUTIONTIME=',ttim

        end if

      else  ! other steps
c....   Zone Header
        if    (nen.eq.2) then                     ! 2 node element
         if(ndm.eq.2) then
          write(ip1,'(a18,i8,a4,i8,a25,a49,a15,e12.5)',err=99)
     +     'ZONE T="FEAP", N=',numnp,', E=',numel1,
     +     ', F=FEPOINT, ET=TRIANGLE',
     +     ', VARSHARELIST=([3,4]=1), CONNECTIVITYSHAREZONE=1',
     +     ', SOLUTIONTIME=',ttim
         else if(ndm.eq.3) then
          write(ip1,'(a18,i8,a4,i8,a25,a51,a15,e12.5)',err=99)
     +     'ZONE T="FEAP", N=',numnp,', E=',numel1,
     +     ', F=FEPOINT, ET=TRIANGLE',
     +     ', VARSHARELIST=([4,5,6]=1), CONNECTIVITYSHAREZONE=1',
     +     ', SOLUTIONTIME=',ttim
         end if
        else if(nen.eq.3) then                    ! 3 node element
          write(ip1,'(a18,i8,a4,i8,a25,a51,a15,e12.5)',err=99)
     +     'ZONE T="FEAP", N=',numnp,', E=',numel1,
     +     ', F=FEPOINT, ET=TRIANGLE',
     +     ', VARSHARELIST=([4,5,6]=1), CONNECTIVITYSHAREZONE=1',
     +     ', SOLUTIONTIME=',ttim

        else if(nen.eq.4.and.ityp.ne.15) then     ! Quadrilateral-Element
          write(ip1,'(a18,i8,a4,i8,a30,a51,a15,e12.5)',err=99)
     +     'ZONE T="FEAP", N=',numnp,', E=', numel1,
     +     ', F=FEPOINT, ET=QUADRILATERAL',
     +     ', VARSHARELIST=([4,5,6]=1), CONNECTIVITYSHAREZONE=1',
     +     ', SOLUTIONTIME=',ttim

        else if(nen.eq.4.and.ityp.eq.15) then     ! Tetrahedron-Element
          write(ip1,'(a18,i8,a4,i8,a30,a51,a15,e12.5)',err=99)
     +     'ZONE T="FEAP", N=',numnp,', E=',numel1,
     +     ', F=FEPOINT, ET=TETAHEDRON',
     +     ', VARSHARELIST=([4,5,6]=1), CONNECTIVITYSHAREZONE=1',
     +     ', SOLUTIONTIME=',ttim

        else if(nen.eq.8) then                    ! hexahedral element
          write(ip1,'(a18,i8,a4,i8,a22,a51,a15,e12.5)',err=99)
     +     'ZONE T="FEAP", N=',numnp,', E=',numel1,
     +     ', F=FEPOINT, ET=BRICK',
     +     ', VARSHARELIST=([4,5,6]=1), CONNECTIVITYSHAREZONE=1',
     +     ', SOLUTIONTIME=',ttim

        else if(nen.eq.9) then                    ! 4 quadrilateral elements (from 9 node element)
          write(ip1,'(a18,i8,a4,i8,a30,a51,a15,e12.5)',err=99)
     +     'ZONE T="FEAP", N=',numnp,', E=',4*numel1,
     +     ', F=FEPOINT, ET=QUADRILATERAL',
     +     ', VARSHARELIST=([4,5,6]=1), CONNECTIVITYSHAREZONE=1',
     +     ', SOLUTIONTIME=',ttim

        else if(nen.eq.16) then                    ! 9 quadrilateral elements (from 16 node element)
          write(ip1,'(a18,i8,a4,i8,a30,a51,a15,e12.5)',err=99)
     +     'ZONE T="FEAP", N=',numnp,', E=',9*numel1,
     +     ', F=FEPOINT, ET=QUADRILATERAL',
     +     ', VARSHARELIST=([4,5,6]=1), CONNECTIVITYSHAREZONE=1',
     +     ', SOLUTIONTIME=',ttim

        else if(nen.eq.27) then                   ! 8 hexahedral elements (from 27 node element)
          write(ip1,'(a18,i8,a4,i8,a22,a51,a15,e12.5)',err=99)
     +     'ZONE T="FEAP", N=',numnp,', E=',8*numel1,
     +     ', F=FEPOINT, ET=BRICK',
     +     ', VARSHARELIST=([4,5,6]=1), CONNECTIVITYSHAREZONE=1',
     +     ', SOLUTIONTIME=',ttim
        end if
      end if

c.....  write nodal coordinates
        if(ityp.eq.1) then ! 3/4/9 node plane stress element 2d/3d
c         X+ux,Y+uy,Z,[X,Y,Z],ux,uy
          do i=1,numnp
            xi = x(1,i)
            if(xi.eq.xit) call prxtie(gtie,i,x,xi,ndm) ! tied nodes
            yi = x(2,i)
            zi = 0.d0
            if(ndm.eq.3) zi = x(3,i)
            if(first) then
              write(ip1,1000)
     +         xi+u(1,i),yi+u(2,i),zi,
     +         xi,       yi       ,zi,
     +         u(1,i),u(2,i),
     +        (s(i,k),k=1,mstv)
            else
              write(ip1,1000)
     +         xi+u(1,i),yi+u(2,i),zi,
     +         u(1,i),u(2,i),
     +        (s(i,k),k=1,mstv)
            end if
          end do

        else if(ityp.eq.2) then  !  2 node beam element 2d/axishell
c         X+ux,Y+uy,[X,Y],ux,uy,rz
          do i=1,numnp
            xi = x(1,i)
            if(xi.eq.xit) call prxtie(gtie,i,x,xi,ndm) ! tied nodes
            yi = x(2,i)
            if(first) then
              write(ip1,1000)
     +         xi+u(1,i),yi+u(2,i),
     +         xi,       yi,
     +         u(1,i),u(2,i),u(3,i)
            else
              write(ip1,1000)
     +         xi+u(1,i),yi+u(2,i),
     +         u(1,i),u(2,i),u(3,i)
            end if
          end do

        else if(ityp.eq.3) then  !  3/4/9 node plate element 2d/3d
c         X,Y,Z+uz,[X,Y,Z],uz,rx,ry
          do i=1, numnp
            xi = x(1,i)
            if(xi.eq.xit) call prxtie(gtie,i,x,xi,ndm)
            yi = x(2,i)
            zi = 0.d0
            if(ndm.eq.3) zi=x(3,i)
            if(first) then
              write(ip1,1000)
     +         xi,yi,zi+u(1,i),
     +         xi,yi,zi,
     +         u(1,i),u(2,i),u(3,i),
     +         (s(i,k),k=1,mstv)
            else
              write(ip1,1000)
     +         xi,yi,zi+u(1,i),
     +         u(1,i),u(2,i),u(3,i),
     +         (s(i,k),k=1,mstv)
            end if
          end do
        else if(ityp.eq.4) then  ! 2D element for cosserat theory
          do i=1,numnp
            xi = x(1,i)
            if(xi.eq.xit) call prxtie(gtie,i,x,xi,ndm)
            yi = x(2,i)
            zi = 0.0d0
            if(first) then
              write(ip1,1000)
     +         xi+u(1,i),yi+u(2,i),zi,
     +         xi,yi,zi,
     +         u(1,i),u(2,i),u(3,i),
     +         (s(i,k),k=1,mstv)
            else
              write(ip1,1000)
     +         xi+u(1,i),yi+u(2,i),zi,
     +         u(1,i),u(2,i),u(3,i),
     +         (s(i,k),k=1,mstv)
            end if
          end do
        else if(ityp.eq.5) then  ! hexahedral element 8/27 nodes
c         X+ux,Y+uy,Z+uz,[X,Y,Z],ux,uy,uz
          do i=1,numnp
            xi = x(1,i)
            if(xi.eq.xit) call prxtie(gtie,i,x,xi,ndm)
            yi = x(2,i)
            zi = x(3,i)
            if(first) then
              write(ip1,1000)
     +         xi+u(1,i),yi+u(2,i),zi+u(3,i),
     +         xi,yi,zi,
     +         u(1,i),u(2,i),u(3,i),
     +         (s(i,k),k=1,mstv)
            else
              write(ip1,1000)
     +         xi+u(1,i),yi+u(2,i),zi+u(3,i),
     +         u(1,i),u(2,i),u(3,i),
     +         (s(i,k),k=1,mstv)
            end if
          end do

        else if(ityp.eq.6) then  ! 3D-brick element for electromechanical/thermomechanical coupling
          do i=1,numnp
            xi = x(1,i)
            if(xi.eq.xit) call prxtie(gtie,i,x,xi,ndm)
            yi = x(2,i)
            zi = x(3,i)
            if(first) then
              write(ip1,1000)
     +         xi+u(1,i),yi+u(2,i),zi+u(3,i),
     +         xi,yi,zi,
     +         u(1,i),u(2,i),u(3,i),u(4,i),
     +         (s(i,k),k=1,mstv)
            else
              write(ip1,1000)
     +         xi+u(1,i),yi+u(2,i),zi+u(3,i),
     +         u(1,i),u(2,i),u(3,i),u(4,i),
     +         (s(i,k),k=1,mstv)
            end if
          end do

        else if(ityp.eq.7) then  ! 3D-brick element for cosserat theory
          do i=1,numnp
            xi = x(1,i)
            if(xi.eq.xit) call prxtie(gtie,i,x,xi,ndm)
            yi = x(2,i)
            zi = x(3,i)
            if(first) then
              write(ip1,1000)
     +         xi+u(1,i),yi+u(2,i),zi+u(3,i),
     +         xi,yi,zi,
     +         u(1,i),u(2,i),u(3,i),u(4,i),u(5,i),u(6,i),
     +         (s(i,k),k=1,mstv)
            else
              write(ip1,1000)
     +         xi+u(1,i),yi+u(2,i),zi+u(3,i),
     +         u(1,i),u(2,i),u(3,i),u(4,i),u(5,i),u(6,i),
     +         (s(i,k),k=1,mstv)
            end if
          end do

        else if(ityp.eq.8) then  ! 3/4 node 5-parameter shell element
c         X+ux,Y+uy,Z+uz,[X,Y,Z],ux,uy,uz,rx,ry
          do i=1,numnp
            xi = x(1,i)
            if(xi.eq.xit) call prxtie(gtie,i,x,xi,ndm)
            yi = x(2,i)
            zi = x(3,i)
            if(first) then
              write(ip1,1000)
     +         xi+u(1,i),yi+u(2,i),zi+u(3,i),
     +         xi,yi,zi,
     +         u(1,i),u(2,i),u(3,i),u(4,i),u(5,i),
     +         (s(i,k),k=1,mstv)
            else
              write(ip1,1000)
     +         xi+u(1,i),yi+u(2,i),zi+u(3,i),
     +         u(1,i),u(2,i),u(3,i),u(4,i),u(5,i),
     +         (s(i,k),k=1,mstv)
            end if
          end do

        else  if(ityp.eq.9) then  ! 2 node 3D-beam element
c         X+ux,Y+uy,Z+uz,[X,Y,Z],ux,uy,uz,rx,ry,rz
          do i=1,numnp
            xi = x(1,i)
            if(xi.eq.xit) call prxtie(gtie,i,x,xi,ndm)
            yi = x(2,i)
            zi = x(3,i)
            if(first) then
              write(ip1,1000)
     +         xi+u(1,i),yi+u(2,i),zi+u(3,i),
     +         xi,yi,zi,
     +         u(1,i),u(2,i),u(3,i),u(4,i),u(5,i),u(6,i)
            else
              write(ip1,1000)
     +         xi+u(1,i),yi+u(2,i),zi+u(3,i),
     +         u(1,i),u(2,i),u(3,i),u(4,i),u(5,i),u(6,i)
            end if
          end do
      else  if(ityp.eq.10) then  ! 3/4 node 6-parameter shell element
c         X+ux,Y+uy,Z+uz,[X,Y,Z],ux,uy,uz,rx,ry,rz
          do i=1,numnp
            xi = x(1,i)
            if(xi.eq.xit) call prxtie(gtie,i,x,xi,ndm)
            yi = x(2,i)
            zi = x(3,i)
            if(first) then
              write(ip1,1000)
     +         xi+u(1,i),yi+u(2,i),zi+u(3,i),
     +         xi,yi,zi,
     +         u(1,i),u(2,i),u(3,i),u(4,i),u(5,i),u(6,i),
     +         (s(i,k),k=1,mstv)
            else
              write(ip1,1000)
     +         xi+u(1,i),yi+u(2,i),zi+u(3,i),
     +         u(1,i),u(2,i),u(3,i),u(4,i),u(5,i),u(6,i),
     +         (s(i,k),k=1,mstv)
            end if
          end do
      else if(ityp.eq.11) then  ! 2 node 3D-beam element incl. warping
c         X+ux,Y+uy,Z+uz,[X,Y,Z],ux,uy,uz,rx,ry,rz,warp
          do i=1,numnp
            xi = x(1,i)
            if(xi.eq.xit) call prxtie(gtie,i,x,xi,ndm)
            yi = x(2,i)
            zi = x(3,i)
            if(first) then
              write(ip1,1000)
     +         xi+u(1,i),yi+u(2,i),zi+u(3,i),
     +         xi,yi,zi,
     +         u(1,i),u(2,i),u(3,i),u(4,i),u(5,i),u(6,i),u(7,i)
            else
              write(ip1,1000)
     +         xi+u(1,i),yi+u(2,i),zi+u(3,i),
     +         u(1,i),u(2,i),u(3,i),u(4,i),u(5,i),u(6,i),u(7,i)
            end if
          end do
      else  if(ityp.eq.12) then  ! 4 node 3D-shell element incl. warping dof of beam3D7dof
c         X+ux,Y+uy,Z+uz,[X,Y,Z],ux,uy,uz,rx,ry,rz,warp
          do i=1,numnp
            xi = x(1,i)
            if(xi.eq.xit) call prxtie(gtie,i,x,xi,ndm)
            yi = x(2,i)
            zi = x(3,i)
            if(first) then
              write(ip1,1000)
     +         xi+u(1,i),yi+u(2,i),zi+u(3,i),
     +         xi,yi,zi,
     +         u(1,i),u(2,i),u(3,i),u(4,i),u(5,i),u(6,i),u(7,i),
     +         (s(i,k),k=1,mstv)
            else
              write(ip1,1000)
     +         xi+u(1,i),yi+u(2,i),zi+u(3,i),
     +         u(1,i),u(2,i),u(3,i),u(4,i),u(5,i),u(6,i),u(7,i),
     +         (s(i,k),k=1,mstv)
            end if
          end do
          else if(ityp.eq.13) then  !plain element, 2-D phase field theory, ndf=1
c         X+ux,Y+uy,Z+uz,[X,Y,Z],ux,uy,Px,Py,Phi        
          do i=1,numnp
            xi = x(1,i) 
c            if(xi.eq.xit) call prxtie(gtie,i,x,xi,ndm)!IM schädlich, funktioniert bei nel=16 nicht, ww??
            yi = x(2,i) 
            zi = 0.0d0
            if(first) then
              write(ip1,1000) 
     +         xi,yi,zi, 
     +         xi,yi,zi,
     +         u(1,i),
     +         (s(i,k),k=1,mstv)
            else  
              write(ip1,1000) 
     +         xi,yi,zi, 
     +         u(1,i),
     +         (s(i,k),k=1,mstv)
            end if
          end do
      else  if(ityp.eq.14) then  !plain element, 2-D phase field theory, ndf=5
c         X+ux,Y+uy,Z+uz,[X,Y,Z],ux,uy,Px,Py,Phi        
          do i=1,numnp
            xi = x(1,i) 
            if(xi.eq.xit) call prxtie(gtie,i,x,xi,ndm)
            yi = x(2,i) 
            zi = 0.0d0
            if(first) then  
              write(ip1,1000) 
     +         xi+u(1,i),yi+u(2,i),zi, 
     +         xi,yi,zi,
     +         u(1,i),u(2,i),u(3,i),u(4,i),u(5,i),
     +         (s(i,k),k=1,mstv)
            else  
              write(ip1,1000) 
     +         xi+u(1,i),yi+u(2,i),zi, 
     +         u(1,i),u(2,i),u(3,i),u(4,i),u(5,i),
     +         (s(i,k),k=1,mstv)
            end if
          end do
      else  if(ityp.eq.15) then  !brick element, 3-D phase field theory, ndf=7
c         X+ux,Y+uy,Z+uz,[X,Y,Z],ux,uy,uz,Px,Py,Pz,Phi        
          do i=1,numnp
            xi = x(1,i) 
            if(xi.eq.xit) call prxtie(gtie,i,x,xi,ndm)
            yi = x(2,i) 
            zi = x(3,i) 
            if(first) then  
              write(ip1,1000) 
     +         xi+u(1,i),yi+u(2,i),zi+u(3,i), 
     +         xi,yi,zi,
     +         u(1,i),u(2,i),u(3,i),u(4,i),u(5,i),u(6,i),u(7,i),
     +         (s(i,k),k=1,mstv)
            else  
              write(ip1,1000) 
     +         xi+u(1,i),yi+u(2,i),zi+u(3,i), 
     +         u(1,i),u(2,i),u(3,i),u(4,i),u(5,i),u(6,i),u(7,i),
     +         (s(i,k),k=1,mstv)
            end if
          end do
        else if(ityp.eq.16) then  ! tetahedron element 4 nodes
c         X+ux,Y+uy,Z+uz,[X,Y,Z],ux,uy,uz
          do i=1,numnp
            xi = x(1,i)
            if(xi.eq.xit) call prxtie(gtie,i,x,xi,ndm)
            yi = x(2,i)
            zi = x(3,i)
            if(first) then
              write(ip1,1000)
     +         xi+u(1,i),yi+u(2,i),zi+u(3,i),
     +         xi,yi,zi,
     +         u(1,i),u(2,i),u(3,i),
     +         (s(i,k),k=1,mstv)
            else
              write(ip1,1000)
     +         xi+u(1,i),yi+u(2,i),zi+u(3,i),
     +         u(1,i),u(2,i),u(3,i),
     +         (s(i,k),k=1,mstv)
            end if
          end do

      end if
c
c.....node-element relation
      if(first) then
        do 300 i=1, numel
          ma = ix(nen1,i)
          if(iplma(ma).eq.0) goto 300
          if      (nen.eq.2) then  ! beam to 3 node element like 1,2,1
            write (ip1,'(3i8)',err=99) ix(1,i),ix(2,i),ix(1,i)

          else if(nen.eq.3) then  ! 3 node element
            write (ip1,'(3i8)',err=99) (ix(jel,i),jel=1,nen)

          else if(nen.eq.4) then  ! 4 node element
            if((ix(3,i).eq.0).and.(ix(4,i).eq.0)) then ! beam
              write (ip1,'(4i8)',err=99) ix(1,i),ix(2,i),ix(2,i),ix(1,i)
            else
              write (ip1,'(4i8)',err=99) (ix(jel,i),jel=1,nen)
            end if

          else if(nen.eq.8) then !  8 node element
            write (ip1,'(8i8)',err=99) (ix(jel,i),jel=1,nen)

          else if (nen.eq.9) then ! 9 node element divided into 4x4 node elements
            write (ip1,'(8i8)',err=99) ix(1,i),ix(5,i),ix(9,i),ix(8,i)
            write (ip1,'(8i8)',err=99) ix(5,i),ix(2,i),ix(6,i),ix(9,i)
            write (ip1,'(8i8)',err=99) ix(8,i),ix(9,i),ix(7,i),ix(4,i)
            write (ip1,'(8i8)',err=99) ix(9,i),ix(6,i),ix(3,i),ix(7,i)

          else if (nen.eq.16) then ! 16 node element divided into nine 4-node elements 
          write (ip1,'(8i8)',err=99) ix(1,i),ix(5,i),ix(13,i),ix(12,i)
          write (ip1,'(8i8)',err=99) ix(5,i),ix(6,i),ix(14,i),ix(13,i)
          write (ip1,'(8i8)',err=99) ix(6,i),ix(2,i),ix(7,i),ix(14,i)
          write (ip1,'(8i8)',err=99) ix(12,i),ix(13,i),ix(16,i),ix(11,i)
          write (ip1,'(8i8)',err=99) ix(13,i),ix(14,i),ix(15,i),ix(16,i)
          write (ip1,'(8i8)',err=99) ix(14,i),ix(7,i),ix(8,i),ix(15,i)
          write (ip1,'(8i8)',err=99) ix(11,i),ix(16,i),ix(10,i),ix(4,i)
          write (ip1,'(8i8)',err=99) ix(16,i),ix(15,i),ix(9,i),ix(10,i)
          write (ip1,'(8i8)',err=99) ix(15,i),ix(8,i),ix(3,i),ix(9,i)

          else if (nen.eq.27) then ! 27 node element divided into 8x8 node elements (Christian Fell TUD)
         write (ip1,'(8i8)',err=99) ix( 1,i),ix(13,i),ix(17,i),ix(16,i),
     +                              ix( 9,i),ix(23,i),ix(27,i),ix(26,i)
         write (ip1,'(8i8)',err=99) ix(13,i),ix( 2,i),ix(14,i),ix(17,i),
     +                              ix(23,i),ix(10,i),ix(24,i),ix(27,i)
         write (ip1,'(8i8)',err=99) ix(16,i),ix(17,i),ix(15,i),ix( 4,i),
     +                              ix(26,i),ix(27,i),ix(25,i),ix(12,i)
         write (ip1,'(8i8)',err=99) ix(17,i),ix(14,i),ix( 3,i),ix(15,i),
     +                              ix(27,i),ix(24,i),ix(11,i),ix(25,i)
         write (ip1,'(8i8)',err=99) ix( 9,i),ix(23,i),ix(27,i),ix(26,i),
     +                              ix( 5,i),ix(18,i),ix(22,i),ix(21,i)
         write (ip1,'(8i8)',err=99) ix(23,i),ix(10,i),ix(24,i),ix(27,i),
     +                              ix(18,i),ix( 6,i),ix(19,i),ix(22,i)
         write (ip1,'(8i8)',err=99) ix(26,i),ix(27,i),ix(25,i),ix(12,i),
     +                              ix(21,i),ix(22,i),ix(20,i),ix( 8,i)
         write (ip1,'(8i8)',err=99) ix(27,i),ix(24,i),ix(11,i),ix(25,i),
     +                              ix(22,i),ix(19,i),ix( 7,i),ix(20,i)

          endif
300     continue
      end if
      first=.false.
      return
c.... errors
99                 write(iow,1010)
      if(ior.lt.0) write(*  ,1010)
      return
c.... format statements
1000  format(8(E15.7,2x),/,
     +       8(E15.7,2x),/,
     +       8(E15.7,2x),/,
     +       8(E15.7,2x)) 
1010  format(' ** ERROR ** on a tape write command for =')
      end
c

