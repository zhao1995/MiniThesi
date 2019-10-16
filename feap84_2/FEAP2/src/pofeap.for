      subroutine pofeap
c--------------------------------------------------------------------
c
c      Purpose Plot first FEAP logo on Graphic screen
c
c     W. Wagner IBNM UH 3/92
c--------------------------------------------------------------------
      USE pdata2
      USE pdata8
      USE vdata
      implicit double precision(a-h,o-z)
      integer*2 imx(21),imy(21),igx(8),igy(8),ixl(2)
      dimension xt(4),yt(5,4)
      character yy*63
      data imx/26,18,36,19, 0, 0,50,46,36,18, 0, 0,22,46,36,19, 9,
     1          13,38,27, 0/
      data imy/ 0,18,36,46,50, 0, 0,19,36,18,26,13, 9,19,36,46,22,
     1           0, 0,27,38/
      data igx/ 0, 0,25,50,50, 0,50,50/
      data igy/ 0,40,50,40,20,20,20, 0/
      data ixl/50,175/
      data xt /0.23,0.23,0.42,0.32/
      data yt /
     1       0.299,0.232,0.166,0.099,0.032,
     2       0.299,0.232,0.166,0.099,0.032,
     3       0.332,0.280,0.228,0.175,0.123,
     4       0.332,0.280,0.228,0.175,0.123/
c.... plot meshes 0=no/1=yes
      imesh = 0
c
      if(iplot.eq.0) return ! after startgr only true for DOS

      istruc = 1
      fopn = .true.
      xl  = 0.02
      yl  = 0.4
      siz = 0.9
      size = 200.0/siz
      call pfeap(0.05d0,0.4d0,0.9d0,2)
c.... write text on screen
      call drawtxt(1,xt(idev),yt(1,idev),3,1,63,
     1'A  F i n i t e  E l e m e n t   A n a l y s i s   P r o g r a m')
c
      call drawtxt(1,xt(idev),yt(2,idev),5,1,63,
     1'          (C)   R. L. T A Y L O R   1 9 8 4 - 2 0 1 5          ')
c
      call drawtxt(1,xt(idev),yt(3,idev),6,1,63,
     1'U n i v e r s i t y  o f  C a l i f o r n i a,  B e r k e l e y')
c
      call drawtxt(1,xt(idev),yt(4,idev),4,1,63,
cww  1'            Parallel Version: (C) W. Wagner                    ')
     1'(C) Institut f. Baustatik  Karlsruhe Institute of Technology   ')
c
      write(yy,2000) versn(1),versn(2)
      call drawtxt(1,xt(idev),yt(5,idev),7,1,63,yy)
c2000 format(7x,a16,'      Last Update  ',a16,5x)
2000  format(14x,a16,6x,a16,12x)
c
c.... draw a border around the mesh
      call pppcol(2)
      call ppbox(0.0003d0, 0.0003d0, 1.279d0, 0.9697d0, 3)
c.... plot meshes
      if(imesh.eq.1) then
c      mesh m
          call pppcol(4)
          dixl = xl*size + ixl(1)
          iyl  = 65
          call dplot((imx(1)+dixl)/size,(imy(1)+iyl)/size+yl,3,0)
          do 110 i = 2,21
           call dplot((imx(i)+dixl)/size,(imy(i)+iyl)/size+yl,2,0)
110       continue
          call dplot((imx(1)+dixl)/size,(imy(1)+iyl)/size+yl,4,0)
c      mesh g
          dixl = xl*size + ixl(2)
          call dplot((igx(1)+dixl)/size,(igy(1)+iyl)/size+yl,3,0)
          do 111 i = 2,8
           call dplot((igx(i)+dixl)/size,(igy(i)+iyl)/size+yl,2,0)
111       continue
          call dplot((igx(1)+dixl)/size,(igy(1)+iyl)/size+yl,4,0)
      endif
c
c.... logos institut
c     call ibnm  ! Hannover
c     call unih
      call unika ! Karlsruhe
c
      call plclos
      return
      end

      subroutine plotini
c--------------------------------------------------------------------
c
c      Purpose: define some initial values for plot
c
c--------------------------------------------------------------------
      USE hpgl1
      USE pback
      USE pdata2
      USE pdata11
      USE pdatap
      USE plinet
      USE plotter
      implicit double precision(a-h,o-z)
      deltax = 0.0
      deltay = 0.0
      ipgl   = 1
      ihpgl  = 1
      imono  = 0
      iback  = 0
      ilinp  = 1
      ibor   = 1
      nexte  = 0
      ibps   = 0
      icps   = 1
      ipan   = 0
      iprin  = 1
      return
      end
