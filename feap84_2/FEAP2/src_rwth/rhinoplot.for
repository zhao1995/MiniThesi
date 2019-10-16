      subroutine rhinopost(x,ix,s,u,nen1,ndm,ndf,ktec,ityp1,loadc)
c----------------------------------------------------------------------
c
c.... Purpose: store data fields for postprocessing in RHINO Format
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
c     loadc      = n3 load case actual values are assigned to
c
c     [rhin,init,n2]  initialize outputfile 'filename'.res,  n1=1 n2=ityp define element 
c                       write geometry to file 'filename'.rhin and close file 
c
c     [rhin,write,n3]     write all values (def+stre) to load case n3      n1=3
c     ([rhin,eigv,n1,n2]  write eigenvector n1 instead of displacements, scaling factor n2 )
c
c     [rhin,close]       close file  out.post.res          n1=2
c
c
c----------------------------------------------------------------------
      USE cdata
      USE comfil
      USE iofile
      USE pnodn
      USE strnam
      USE isogeo      
      implicit double precision (a-h,o-z)
      logical ex, first
      character*1 antw,eltyp*69,eltyp1*69,text*111,resultnames*25
      character*229 fname1
      character     strsust*15
      dimension eltyp(2)  ! Element-Typ
      dimension x(ndm,*)   ! Knoten Koordinaten
      dimension ix(nen1,*) ! Element-Knoten-Beziehung
      dimension s(numnp,*) ! Knoten Spannungen
      dimension u(ndf,*)   ! Knoten Verschiebungen
      dimension strsust(25)! Name der Knotenspannung
      dimension resultnames(10)
c      include 'cdata.h'
c      include 'comfil.h'
c      include 'iofile.h'
c      include 'pnodn.h'
c      include 'strnam.h' 
c      include 'isogeo.h'
c      common m(100000)
      save numel1,first,ityp,fname1 
   
      data ip1/23/
      data ntyp/2/
      data eltyp/ 
     1'isogeometric Reissner-Mindlin Shell ndf = 5,6  ndm = 4        ',
     2'dummy                                                         '/ 
      
      data resultnames/
     1'"Normal force n11"', 
     2'"Normal force n22"',  
     3'"Normal force n12"',
     4'"Normal force n12"',
     5'"Bending Moment m11"',
     6'"Bending Moment m22"', 
     7'"Bending Moment m12"',
     8'"Shear force q13"',
     9'"Shear force q23"', 
     1'"Displacement"'/

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

c...  open file for RHINO geometry data
c      fname1 = fres
c      call addext (fname1,'rhin ')
      fname1 = 'out.georhino.txt'
      inquire(FILE=fname1,EXIST=ex)   ! Pruefen ob Datei schon existiert
      if(ex) then                    ! => Datei existiert
        write (*,1001) 'File ', fname1, ' exists! '
1001    format(a10,/,a229,/,a34)
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
c      text='Files for Postprocessing opened, Element: '  
      eltyp1=eltyp(ityp)
c      text(42+1:42+69)=eltyp1(1:69) 
c      call drawmess(text,1,-3)
      write(*,1001) 'File ', fname1, ' contains the geometry for RHINO'
c...  write geometry
c
c      
c.....Geometry header
      write(ip1,'(a7)') 'ND-COOR'
c.... Write down list of control points
1011  format(a5,i8,a5,ES22.14E3,a5,ES22.14E3,a5,ES22.14E3) 
      do i=1,numnp
        xi = x(1,i) 
csk        if(xi.eq.xit) call rhinotie(m(nip),i,x,xi,ndm) ! tied nodes
        if(xi.eq.xit) call rhinotie(gtie,i,x,xi,ndm) ! tied nodes
        yi = x(2,i) 
        zi = x(3,i)  
        write(ip1,1011) 'NODE ',i,' X ',xi,' Y ',yi,' Z ',zi
      end do
c...  Write all neccessary patch information
      write(ip1,'(a1)') ''
      write(ip1,'(a44)') '!############################################'
      write(ip1,'(a44)') '!#########      NURBS-BLOCK       ###########'
      write(ip1,'(a44)') '!############################################'
c...  loop over all patches
      nj1 = 0
      nj2 = 0
      nj = 0
      do i = 1,NURnpatch
        NURn = nNURnmpq(i,1)
        NURm = nNURnmpq(i,2)
        NURp = nNURnmpq(i,3)
        NURq = nNURnmpq(i,4)
        write(ip1,'(a12,i8,a11)') 'NURBS_PATCH ',i,' : NURBS_2D'
        write(ip1,'(a27,i8)')     ' CTRL_PTS = CTRL_PTS_NODES ',i
        write(ip1,'(a9,i8)') ' NCTRL = ',NURn-1
        write(ip1,'(a9,i8)') ' MCTRL = ',NURm-1
        write(ip1,'(a9,i8)') ' PDEG =  ',NURp
        write(ip1,'(a9,i8)') ' QDEG =  ',NURq
c...    write knot vectors
        write(ip1,'(a9$)') ' UKNOT = '
        do j1 = 1,NURn + NURp
          write(ip1,'(ES22.14E3,a2$)') rNURknv1(j1+nj1),', '
        end do
        nj1 = nj1 + NURn + NURp + 1
        write(ip1,'(ES22.14E3)') rNURknv1(nj1)
        
        write(ip1,'(a9$)') ' VKNOT = '
        do j2 = 1,NURm + NURq
          write(ip1,'(ES22.14E3,a2$)') rNURknv2(j2+nj2),', '
        end do
        nj2 = nj2 + NURm + NURq + 1
        write(ip1,'(ES22.14E3)') rNURknv2(nj2)
        write(ip1,'(a42)')'!=========================================='
c....   write Control points of patch under consideration and the corresponding weights        
        write(ip1,'(a15,i8)')     'CTRL_PTS_NODES ',i
        do j = 1,NURn*NURm
c          xi = x(1,j) 
c          if(xi.eq.xit) call rhinotie(m(nip),j,x,xi,ndm) ! tied nodes
           write(ip1,'(a9,i8,a3,ES22.14E3)') ' NODE_ID ',j+nj,' W ',
     +     x(4,j+nj)
        end do
        nj = nj + NURn*NURm
        write(ip1,'(a1)') ''
      end do
      
c...  close geometry file  
      close(ip1)

      
c...  open results file      
c      fname1 = fres
c      call addext (fname1,'\out.post.res')
      fname1 = 'out.post.res'
      inquire(FILE=fname1,EXIST=ex)   ! Pruefen ob Datei schon existiert
      if(ex) then                    ! => Datei existiert
        write (*,1001) 'File ', fname1, ' exists! '
        write (*,1002)

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
      write(*,1001) 'File ', fname1, ' opened to write results'
      
      return  
c...  close file for TecPlot Data
2     close(ip1)
      write(*,1001) 'File ', fname1, ' for results closed'
      return
c
c.... write all data for step
3     continue

      if(first) then   
c....   Rhino Post Results file Header rausschreiben
c....   name of stress is in resultnames and strsus
        !do i = 1,25
        !  if(strsus(i).eq.'               ') then
        !    strsust(i) = resultnames(i)
        !  else
        !    strsust(i) = strsus(i) 
        !  end if
        !end do
        write (ip1,'(a27)') 'Rhino Post Results File 1.0'
      end if

c...  write deformations in control points
      write(ip1,'(a7,a25,a13,i3,a8,a7)') 
     +'Result ',resultnames(10),' "Load Case" ',loadc,
     +' Vector ','OnNodes' 
      if (ityp.eq.1) then  !isogeometric Reissner-Mindlin shell
        write(ip1,'(a6)') 'Values'
        do i=1,numnp
c     !     write(ip1,'(a2,i8,a1,ES22.14E3,a1,ES22.14E3,a1,ES22.14E3,a1,
c     !+                ES22.14E3,a1,ES22.14E3,a1,ES22.14E3)')
c     !+         '  ',i,' ',u(1,i),' ',u(2,i),' ',u(3,i),'  ',u(4,i),' '
c     !+                   ,u(5,i),' ',u(6,i)
        
          write(ip1,'(a2,i8,a1,ES22.14E3,a1,ES22.14E3,a1,ES22.14E3)')
     +         '  ',i,' ',u(1,i),' ',u(2,i),' ',u(3,i)
        
        
        end do
        write(ip1,'(a10)') 'End Values'
      end if


c...  write stresses in control points
      
      do j = 1,mstv
        write(ip1,'(a7,a25,a13,i3,a8,a7)') 
     +    'Result ',resultnames(j),' "Load Case" ',loadc,
     +    ' Scalar ','OnNodes' 
        if (ityp.eq.1) then  !isogeometric Reissner-Mindlin shell
        write(ip1,'(a6)') 'Values'
        do i=1,numnp
          write(ip1,'(a2,i8,a1,ES22.14E3)')
     +         '  ',i,' ',s(i,j)
        
        
          end do
          write(ip1,'(a10)') 'End Values'
        end if  
      end do
      

      first=.false.      
      return
c.... errors
99                 write(iow,1010) 
      if(ior.lt.0) write(*  ,1010) 
      return
c.... format statements
1000  format(32(E15.7,2x)) 
1010  format(' ** ERROR ** on a tape write command for =')
      end
c
      subroutine rhinotie(ip,ni,x,x1,ndm)
c-----------------------------------------------------------------------
c.... look for tied nodes, compare prttie                              |
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer*4 ni,n
      dimension ip(*),x(*)
      n  = ip(ni)
      x1 = x((n-1)*ndm+1)  
      return
      end
