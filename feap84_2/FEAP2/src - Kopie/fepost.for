      subroutine fepost(isw,b,ix,ndm,ndf,nen1,k1,indexn,numnpn)
c----------------------------------------------------------------------
c.... Purpose: store data fields for postprocessing
c
c     Inputs:
c         isw = 1      - initialize outputfile, write system data
c         isw = 2      - write displacements
c         isw = 3      - write nodal reactions
c         isw = 4      - write eigenvector n1
c         isw = 5      - write nodal stresses
c         isw = 6      - write nodal warping values
c         isw = 7      - close files              
c         b(*)         - array to write
c         ix(nen1,*)   - Element nodal connections of mesh
c         ndm          - Spatial dimension of mesh
c         ndf          - Number dof/node
c         nen1         - Dimension for ix array
c         k1           - No. of eigenvector
c         indexn(numnp)- array to correct node numbers
c                        new node number i is: i-indexn(i)        
c
c     Outputs:  -     
c
c     Comments: 
c       1)numnpn = numnp: original mesh, indexn(i) = 0 
c         when tie is used:
c         list of nodes: tied nodes have coordinate x1=-999
c         list of elements with active nodes, thus tied nodes are not used
c         for disp etc:  tied nodes have correct values 
c
c       2)numnpn < numnp: new mesh, indexn(i) .ne. 0 
c        when tie is used:
c        renumbering of nodes and nodes on elements
c        tied nodes are completely densed out
c        list of nodes: renumbered, withhout tied nodes
c        list of elements: renumbered, without tied nodes
c        for disp etc:  without tied nodes  
c
c       1) and 2) works only for continous(!) node numbering
c
c
c       Store data in File Rname.pos
c       
c     W. Wagner BS KIT 08/10
c----------------------------------------------------------------------
      USE arcl
      USE cdata
      USE comfil
      USE fdata
      USE iofile
      USE mdata
      USE pdata3
      USE prlod
      implicit double precision (a-h,o-z)
      character*229 fname1
      character*4 ct(7)
      character*80 text
      dimension b(*),ix(nen1,*),indexn(numnp)
      data ip1/23/,
     +      ct/'init','disp','reac','eigv','stre','warp','clos'/

      propl = prop

      if(arcf) propl = propl*rlnew

      goto(100,200,300,400,500,600,700) isw
c.... open file and write system
100   fname1 = fres
      call addext(fname1,'pos ')
      open(ip1,file=fname1,status='unknown',form='formatted')
      rewind ip1
      call drawmess('File for Postprocessing opened',1,-2)
c     basic data 
      write(ip1,1001,err=990) numnpn,numel,nummat,ndm,ndf,nen
c
c     nodal coordinates:  field x (ndm,numnp)
      knold = 0
      do k = 1,numnp
         call newnode(k,kn,knold,indexn,numnp)
         if(kn.gt.0) then
           ia = (k-1)*ndm+1
           ie = ia+ndm-1
           write(ip1,9000,err=990) (b(i),i=ia,ie)
         end if  
      end do  
c
c     element connections:  field ix (nen1,numel) 
c     1-nen = nodes,nen+1=nt1,nen+2=nt2,nen+3=nt3,nen1=nen+4=mat.number)

      do k = 1,numel
        write(ip1,9001,err=990) (ix(i,k)-indexn(ix(i,k)),i=1,nen),
     +       (ix(i,k),i=nen+1,nen1)
      end do

      return
c
c.... save current displacement field u(ndf,numnp)
200   write(text,1002) propl
      write( ip1,2000) text

      knold = 0
      do k = 1,numnp
         call newnode(k,kn,knold,indexn,numnp)
         if(kn.gt.0) then
           ia = (k-1)*ndf+1
           ie = ia+ndf-1
           write(ip1,9002,err=990)  (b(i),i=ia,ie)
         end if  
      end do  

      return
c
c.... save current reaction field r(ndf,numnp)
300   write(text,1003) propl
      write( ip1,2000) text

      knold = 0
      do k = 1,numnp
         call newnode(k,kn,knold,indexn,numnp)
         if(kn.gt.0) then
           ia = (k-1)*ndf+1
           ie = ia+ndf-1
           write(ip1,9003,err=990)  (b(i),i=ia,ie)
         end if  
      end do   

      return
c
c.... save current eigenvector field k1 ev(ndf,numnp)
400   write(text,1004) propl,k1
      write( ip1,2000) text

      knold = 0
      do k = 1,numnp
         call newnode(k,kn,knold,indexn,numnp)
         if(kn.gt.0) then
           ia = (k-1)*ndf+1
           ie = ia+ndf-1
           write(ip1,9004,err=990)  (b(i),i=ia,ie)
         end if  
      end do  

      return
c
c.... save current nodal stress state s(numnp*npstr-1) npstr=27
500   write(text,1005) propl
      write( ip1,2000) text

      npstr1 = npstr-1
      ii2    = numnp*npstr1
      knold = 0
      do k = 1,numnp
         call newnode(k,kn,knold,indexn,numnp)
         if(kn.gt.0) then
           ii1 = (k-1)+1
           write(ip1,9005,err=990) (b(i),i=ii1,ii2,numnp)
         end if
      end do  

      return
c
c.... save nodal warping values =  stre,5
600   write(text,1006) propl
      write( ip1,2000) text

      knold = 0
      do k = 1,numnp
         call newnode(k,kn,knold,indexn,numnp)
         if(kn.gt.0) then
           write(ip1,9006,err=990)  b(k)
         end if
      end do   

      return
c
c.... close file
700   close(ip1)
      call drawmess('File for Postprocessing closed',1,-2)
      return
c
c.... errors
990                write(iow,1000) ct(isw)
      if(ior.lt.0) write(*  ,1000) ct(isw)
      return
c
c.... format statements
1000  format(' ** ERROR ** on a tape write command for =',a4)
1001  format('numnp=',i5,' numel=',i5,' nummat=',i5,' ndm =',i5,
     +       ' ndf =',i5,' nen =',i5)
1002  format('disp   at prop. load = ',e12.5)
1003  format('reac   at prop. load = ',e12.5)
1004  format('eigv   at prop. load = ',e12.5,' eigv no. ',i3)
1005  format('stre   at prop. load = ',e12.5)
1006  format('warp   at prop. load = ',e12.5)
2000  format(a80)
c.... formats for ascii file   
9000  format('nodes   ',3f14.8)  ! coordinates, old node ,new node  
9001  format('elmts   ',12i6)    ! nel<=9 
9002  format('displ.  ',7e12.5)  ! ndf<=7 
9003  format('react.  ',7e12.5)  ! ndf<=7  
9004  format('eigv.   ',7e12.5)  ! ndf<=7  
9005  format('stre    ',6e12.5)  ! npstr = 27      
9006  format('warp    ',1e15.8)  ! 
      end
c
c-----------------------------------------------------------------------
      subroutine newmesh(x,ndm,numnp,indexn,numnpn)
c-----------------------------------------------------------------------
c
c     Purpose: Calculate new mesh
c 
c     Inputs:
c        x(ndm,numnp)  - nodal coordinates 
c
c     Outputs:
c        numnpn         - new total number of nodes
c        indexn(numnp)  - new node number i is: i-indexn(i)        
c  
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension indexn(numnp)  
      dimension x(ndm,numnp)
      data bl/-999.d0/
      
c     nodal coordinates:
      numnpn=0
      do i=1,numnp
        xi=x(1,i) 
        if(xi.eq.bl) then
          indexn(i)=indexn(i-1)+1
        else
          if (i.eq.1) then
            indexn(i)=0
          else
            indexn(i)=indexn(i-1)
          endif
          numnpn=numnpn+1
        end if
      end do

c      do i = 1,numnp
c         newnode= i-indexn(i)
c         write(*,*) 'inode ',i,' index ',indexn(i),' new node ', newnode
c      end do     


      return
      end

c                
c----------------------------------------------------------------------
c
      subroutine newnode(k,kn,knold,indexn,numnp)
c----------------------------------------------------------------------
c
c     Purpose: Find new node number with respect to TIE
c
c     Inputs: 
c       k      - old node number
c       knold  - last new node number
c
c     Outputs
c       kn     - new node number  
c
c     W. Wagner BS UKA 04/09
c----------------------------------------------------------------------
c
      integer indexn(numnp)
      integer k,kn,knold 
      
      kk = indexn(k) 
      kn = k -kk
      if(kn.eq.knold) then 
        kn = -1      ! no new number, skip node
      else 
        knold = kn   ! new node number
      end if
      return
      end
   
