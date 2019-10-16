      SUBROUTINE GOPEN
     &       (XMIN,XMAX,XTICK,YMIN,YMAX,YTICK,XNAME,YNAME,IGRID,ICROSS)
c-----------------------------------------------------------------------
c     Purpose: switches the screen to graphics mode,
c              defines a coordinate system , draws a frame,axes and grid 
c
c     Inputs:
c      XMIN:   lowest value on the x axis
c      XMAX:   highest value on the x axis
c      XTICK:  separation of tick marks on the x axis,
c              XTICK= 0.D0: no tick marks
c              XTICK=-1.D0: automatic determination of separation
c      YMIN:   lowest value on the y axis
c      YMAX:   highest value on the y axis
c      YTICK:  separation of tick marks on the y axis,
c              YTICK= 0.D0: no tick marks
c              YTICK=-1.D0: automatic determination of separation
c      XNAME:  labelling of the x axis
c        ' ' : no labelling of the axis and the tick marks
c      YNAME:  labelling of the y axis
c        ' ' : no labelling of the axis and the tick marks
c      IGRID:  IGRID=0 : grid lines not shown
c              IGRID=1 : grid lines shown
c      ICROSS: ICROSS=0 : axes not shown
c              ICROSS=1 : axes shown, intersecting at the origin
c
c-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O,Q-Z)
      CHARACTER*(*) XNAME,YNAME
      integer open

      COMMON/gxINFO/OPEN,MAXCOL,XAUF,YAUF,XCHHI,XCHWID

      open = 0
      CALL GINIT(1)

c.... move   right 0.05 up 0.1 ww 3/95
c     XWMIN=0.15D0/1.3d0
      XWMIN=0.20D0/1.3d0

c     IF (YNAME.EQ.' '.OR.YTICK.EQ.0.D0) XWMIN=0.D0
      IF (YNAME.EQ.' '.OR.YTICK.EQ.0.D0) XWMIN=0.05D0

c     XWMAX=0.92D0/1.3d0
      XWMAX=0.97D0/1.3d0

c     YWMIN=0.15D0/1.3d0
      YWMIN=0.25D0/1.3d0

c     IF (XTICK.EQ.0.D0) YWMIN=0.075D0/1.3d0
      IF (XTICK.EQ.0.D0) YWMIN=0.125D0/1.3d0

c     IF (XNAME.EQ.' ')  YWMIN=0.D0
      IF (XNAME.EQ.' ')  YWMIN=0.1D0

c     YWMAX=0.87D0/1.3d0
      YWMAX=0.97D0/1.3d0

      xchhi  = 0.02
      xchwid = 0.02

      CALL GWINDO(XMIN,XMAX,YMIN,YMAX,XWMIN,XWMAX,YWMIN,YWMAX)
c
      XCROSS=0.D0
      YCROSS=0.D0

      CALL GCHART(XTICK,XCROSS,YTICK,YCROSS,XNAME,YNAME,IGRID,ICROSS,1)

      END
c
      SUBROUTINE GINIT(IMODE)
c---------------------------------------------------------------------
c
c     Purpose:  switches the screen between graphics and text.
c
c     When the screen is switched to graphics mode
c     GWINDO is called in order to
c     define a window with the corners (0,0) and (1,1) in both
c     screen and world coordinates. By this definition, the screen,
c     window, and world coordinates become identical.
c
c     Inputs:
c     IMODE    > 0 : switch to graphics mode
c              = 0 : switch back to text mode
c
c---------------------------------------------------------------------
      USE pback
      USE pdata2
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      PARAMETER (XLOW=0.D0,XHIGH=1.D0,YLOW=0.D0,YHIGH=1.D0)
      INTEGER OPEN,MAXCOL
      REAL*8 XAUF,YAUF,XCHHI,XCHWID
      COMMON/gxINFO/OPEN,MAXCOL,XAUF,YAUF,XCHHI,XCHWID
cww   common /pswit/ imod
c
      IF (IMODE.GT.0) THEN
        IF (OPEN.NE.1) THEN
c....     open workstation
cww         if(imod.eq.1) then
            iback  = 1
cww         iclear = 0
            call plopen
            open = 1
cww         end if
        END IF

c...    Define window
        CALL GWINDO(XLOW,XHIGH,YLOW,YHIGH,XLOW,XHIGH,YLOW,XHIGH)
      ELSE
c....   close workstation
cwww    if(imod.eq.1) then
           call plclos
cww     end if
        open = 0
      END IF
      return
      end
c
      SUBROUTINE GWINDO(XMIN,XMAX,YMIN,YMAX,XWMIN,XWMAX,YWMIN,YWMAX)
c-----------------------------------------------------------------------
c
c     Purpose: defines world coordinates for a screen window 
c              or the whole screen.
c
c     Inputs:
c      XMIN:   left   margin of window, world coordinates
c      XMAX:   right  margin of window, world coordinates
c      YMIN:   bottom margin of window, world coordinates
c      YMAX:   top    margin of window, world coordinates
c      XWMIN:  left   margin of window, screen coordinates (0-100)
c      XWMAX:  right  margin of window, screen coordinates (0-100)
c      YWMIN:  bottom margin of window, screen coordinates (0-100)
c      YWMAX:  top    margin of window, screen coordinates (0-100)
c
c     COMMON/gxWIND/
c     XLO  =  left   margin of window, world coordinates
c     XHI  =  right  margin of window, world coordinates
c     YLO  =  bottom margin of window, world coordinates
c     YHI  =  top    margin of window, world coordinates
c     XWLO =  left   margin of window, screen coordinates (0-1)
c     XWHI =  right  margin of window, screen coordinates (0-1)
c     YWLO =  bottom margin of window, screen coordinates (0-1)
c     YWHI =  top    margin of window, screen coordinates (0-1)
c
c-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O,Q-Z)
c
      DOUBLE PRECISION XLO,XHI,YLO,YHI,XWLO,XWHI,YWLO,YWHI

      COMMON/gxWIND/XLO,XHI,YLO,YHI,XWLO,XWHI,YWLO,YWHI

      IF ( (XMIN.NE.XMAX)  .AND.  (YMIN.NE.YMAX) ) THEN
        XLO=XMIN
        XHI=XMAX
        YLO=YMIN
        YHI=YMAX
      END IF

      IF ((XWMIN.LT.XWMAX) .AND. (YWMIN.LT.YWMAX) .AND.
     &    (XWMIN.GE.0.D0)  .AND. (XWMAX.LE.1.D0)  .AND.
     &    (YWMIN.GE.0.D0)  .AND. (YWMAX.LE.1.D0)  ) THEN
       XWLO=XWMIN
       XWHI=XWMAX
       YWLO=YWMIN
       YWHI=YWMAX
      END IF
c
      END
c
      SUBROUTINE GCHART
     &       (XTICK,XCROSS,YTICK,YCROSS,XNAME,YNAME,IGRID,ICROSS,IWHERE)
c-----------------------------------------------------------------------
c
c     Purpose: draws a frame, crossed axes and a grid.
c
c     Inputs:
c      XTICK:  separation of tick marks on the x axis,
c              XTICK= 0.D0: no tick marks
c              XTICK=-1.D0: automatic determination of separation
c      XCROSS: x coordinate of the axes intersection
c      YTICK:  separation of tick marks on the y axis,
c              YTICK= 0.D0: no tick marks
c              YTICK=-1.D0: automatic determination of separation
c      YCROSS: y coordinate of the axes intersection
c      XNAME:  labelling of the x axis
c        ' ' : no labelling of the axis and the tick marks
c      YNAME:  labelling of the y axis
c        ' ' : no labelling of the axis and the tick marks
c      IGRID:  IGRID=0 : grid lines not shown
c              IGRID=1 : grid lines shown
c      ICROSS: ICROSS=0 : axes not shown
c              ICROSS=1 : axes shown
c      IWHERE: IWHERE=1 : labelling below and to the left
c              IWHERE=2 : labelling below and to the right
c              IWHERE=3 : labelling above and to the right
c              IWHERE=3 : labelling above and to the left
c
c-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c
      CHARACTER*(*) XNAME,YNAME
c
      CHARACTER*6 FORMS(0:3)
      CHARACTER AXFORM*6,LGFORM*3
      INTEGER gxSKIP,GTRANS
c
      DOUBLE PRECISION XLO,XHI,YLO,YHI,XWLO,XWHI,YWLO,YWHI,gxTICK
      COMMON/gxWIND/XLO,XHI,YLO,YHI,XWLO,XWHI,YWLO,YWHI

      DATA FORMS/'(F6.3)','(F6.2)','(F6.1)','(F6.0)'/
c
      CALL GINFO(ISOPEN,XRES,YRES,NCOLORS,CHWID,CHHIGH)

      CALL pppcol(3)
c
c.... Determine XTICK and YTICK
      XXTICK=XTICK
      IF (XTICK.EQ.-1.D0) XXTICK=gxTICK(XHI,XLO)
      IF (XXTICK.NE.0) THEN
        NXTCKL=(XLO-XCROSS)/XXTICK
        NXTCKH=(XHI-XCROSS)/XXTICK
        NXTCKD=SIGN(1,NXTCKH-NXTCKL)
        NXTCKS=ABS(NXTCKH-NXTCKL)
        IF (NXTCKS.GT.200.OR.XXTICK.LT.0) XXTICK=0
      END IF
c
      YYTICK=YTICK
      IF (YTICK.EQ.-1.D0) YYTICK=gxTICK(YHI,YLO)
      IF (YYTICK.NE.0) THEN
        NYTCKL=(YLO-YCROSS)/YYTICK
        NYTCKH=(YHI-YCROSS)/YYTICK
        NYTCKD=SIGN(1,NYTCKH-NYTCKL)
        NYTCKS=ABS(NYTCKH-NYTCKL)
        IF (NYTCKS.GT.200.OR.YYTICK.LT.0) YYTICK=0
      END IF

c.... Draw grid
      IF (IGRID.NE.0.AND.XXTICK.NE.0) THEN
        DO 10 N=NXTCKL,NXTCKH,NXTCKD
          CALL GLINE(XCROSS+N*XXTICK,YLO,XCROSS+N*XXTICK,YHI,4,3)
10      CONTINUE
      END IF

      IF (IGRID.NE.0.AND.YYTICK.NE.0) THEN
         DO 20 N=NYTCKL,NYTCKH,NYTCKD
           CALL GLINE(XLO,YCROSS+N*YYTICK,XHI,YCROSS+N*YYTICK,4,3)
20       CONTINUE
      END IF

c.... Draw frame
      CALL GLINE(XLO   ,YLO   ,XHI   ,YLO   ,4,1)
      CALL GLINE(XHI   ,YLO   ,XHI   ,YHI   ,4,1)
      CALL GLINE(XHI   ,YHI   ,XLO   ,YHI   ,4,1)
      CALL GLINE(XLO   ,YHI   ,XLO   ,YLO   ,4,1)

c.... Draw axes
      IF (ICROSS.GT.0) THEN
        CALL GLINE(XLO ,YCROSS,XHI   ,YCROSS,1,1)
        CALL GLINE(XCROSS,YLO   ,XCROSS,YHI   ,1,1)
      END IF
 
c.... Label axes
      XHSTEP=CHWID /2
      YHSTEP=CHHIGH/2
cww   XQSTEP=XHSTEP/2
      YQSTEP=YHSTEP/2
c
      IF (XNAME.NE.' ') THEN
        XX=(XWLO+XWHI-CHWID*float(LEN(XNAME)))/2.d0
        IF (IWHERE.LT.3) THEN
          YY=MAX(YWLO-3.d0*YHSTEP,0.D0)
          IF (XXTICK.NE.0) YY=YY-3*YHSTEP
        ELSE
          YY=YWHI+YHSTEP
          IF (XXTICK.NE.0) YY=YY+3*YHSTEP
       END IF
       call pppcol(1)
       CALL GWRITE(1,XX,YY,XNAME,1)
      END IF
c
      IF (YNAME.NE.' ') THEN
       IF (IWHERE.NE.2.AND.IWHERE.NE.3) THEN
         XX=MAX(XWLO-6.d0*CHWID+0.02,0.02D0)
       ELSE
         XX=XWHI+XHSTEP+0.02
       END IF
       YY=YWHI+YHSTEP
       CALL GWRITE(1,XX,YY,YNAME,1)
      END IF

c.... Tick marks on the x axis
      ISTAT=GTRANS(XWLO,YWLO+YHSTEP,XX,YBIG1,-3)
      ISTAT=GTRANS(XWLO,YWHI-YHSTEP,XX,YBIG2,-3)
      YSML1=(YLO+YBIG1)/2
      YSML2=(YHI+YBIG2)/2
      IF (XXTICK.NE.0) THEN
        XT1=XWHI
        XS0=XWLO
        IF (XNAME.EQ.' ') XS0=XT1
        IF (IWHERE.LT.3) THEN
          YEB1=YWLO-3*YHSTEP
        ELSE
          YEB1=YWHI+YHSTEP
          IF (IWHERE.NE.4)XS0=XWLO+MAX(0.d0,FLOAT((LEN(YNAME)-5))*CHWID)
        END IF
        YEB2=YEB1-  YQSTEP
        YEB3=YEB1+  YQSTEP
        AXMAX=MAX(ABS(XLO),ABS(XHI))
        LOGAX=LOG10(AXMAX)
        IF (AXMAX.LT.1) LOGAX=LOGAX-1
        AMFA=1
        IF (LOGAX.LT.0.OR.LOGAX.GT.3.OR.XXTICK.LT.AXMAX*.01D0) THEN
          AMFA=10.D0**(LOGAX)
          WRITE(LGFORM,'(I3)') LOGAX
          ISTART=gxSKIP(LGFORM,1,' ')
          XT1=XT1-(7-ISTART)*CHWID
          IF (XNAME.NE.' ') THEN
            CALL GWRITE(1,XT1    ,YEB2,'*10'           ,1)
            CALL GWRITE(1,XT1+3*CHWID,YEB3,LGFORM(ISTART:3),1)
          END IF
          LOGAX=0
        END IF
        DO 30 N=NXTCKL,NXTCKH,NXTCKD
          XX=XCROSS+N*XXTICK
          AXX=XX/AMFA
          WRITE(AXFORM,FORMS(LOGAX)) AXX
          IEND=gxSKIP(AXFORM,-1,'0')
          IF (AXFORM(IEND:IEND).EQ.'.') IEND=IEND-1
          ISTART=gxSKIP(AXFORM(1:IEND),1,' ')
          ILENG=IEND-ISTART+1
          ISTAT=GTRANS(XX,YLO   ,XS,YS,3)
          XS1=XS- ILENG *XHSTEP
          XT0=XS+(ILENG+1)*XHSTEP
          IF (XS1.GT.XS0.AND.XT0.LT.XT1) THEN
            CALL GWRITE(1,XS1,YEB1,AXFORM(ISTART:IEND),1)
            CALL GLINE(XX,YLO,XX,YBIG1,1,1)
            CALL GLINE(XX,YHI,XX,YBIG2,1,1)
            XS0=XT0
          ELSE
            CALL GLINE(XX,YLO,XX,YSML1,1,1)
            CALL GLINE(XX,YHI,XX,YSML2,1,1)
          END IF
30      CONTINUE
      END IF

c.... Tick marks on the y axis
      ISTAT=GTRANS(XWLO+XHSTEP,YWLO,XBIG1,YY,-3)
      ISTAT=GTRANS(XWHI-XHSTEP,YWLO,XBIG2,YY,-3)
      XSML1=(XLO+XBIG1)/2
      XSML2=(XHI+XBIG2)/2
      IF (YYTICK.NE.0) THEN
        YS0=YWLO
        YT1=YWHI
        IF (YNAME.EQ.' ') YS0=YT1
        YEB1=YWHI -3*YHSTEP
        YEB2=YEB1-  YQSTEP
        YEB3=YEB1+  YQSTEP
        IF (IWHERE.NE.2.AND.IWHERE.NE.3) THEN
          XPS=XWLO-13*XHSTEP
        ELSE
          XPS=XWHI+XHSTEP
        END IF
        AXMAX=MAX(ABS(YLO),ABS(YHI))
        LOGAX=LOG10(AXMAX)
        IF (AXMAX.LT.1) LOGAX=LOGAX-1
        AMFA=1
        IF (LOGAX.LT.0.OR.LOGAX.GT.3.OR.YYTICK.LT.AXMAX*.01D0) THEN
          AMFA=10.D0**(LOGAX)
          WRITE(LGFORM,'(I3)') LOGAX
          ISTART=gxSKIP(LGFORM,1,' ')
          YT1=YEB2
          IF (YNAME.NE.' ') THEN
            CALL GWRITE(1,XPS    ,YEB2,'*10'           ,1)
            CALL GWRITE(1,XPS+3*CHWID,YEB3,LGFORM(ISTART:3),1)
          END IF
          LOGAX=0
        END IF
        DO 40 N=NYTCKL,NYTCKH,NYTCKD
          YY=YCROSS+N*YYTICK
          AYY=YY/AMFA
          WRITE(AXFORM,FORMS(LOGAX)) AYY
          IEND=gxSKIP(AXFORM,-1,'0')
          IF (AXFORM(IEND:IEND).EQ.'.') IEND=IEND-1
          ISTART=gxSKIP(AXFORM(1:IEND),1,' ')
          ILENG=IEND-ISTART+1
          ISTAT=GTRANS(XLO,YY   ,XS,YS,3)
          YS1=YS- YQSTEP
          YT0=YS+ CHHIGH
          IF (YS1.GT.YS0.AND.YT0.LT.YT1) THEN
            IF (IWHERE.NE.2.AND.IWHERE.NE.3) THEN
              XP1=XPS+(6-ILENG)*CHWID
            ELSE
              XP1=XPS
            END IF
            CALL GWRITE(1,XP1,YS1,AXFORM(ISTART:IEND),1)
            CALL GLINE(XLO,YY,XBIG1,YY,1,1)
            CALL GLINE(XHI,YY,XBIG2,YY,1,1)
            YS0=YT0
          ELSE
            CALL GLINE(XLO,YY,XSML1,YY,1,1)
            CALL GLINE(XHI,YY,XSML2,YY,1,1)
          END IF
40     CONTINUE
      END IF
c
      END
c
      DOUBLE PRECISION FUNCTION gxTICK(ZHI,ZLO)
c-----------------------------------------------------------------------
c
c     Purpose: determines separation of tick marks, ~ 10 marks per axis
c
c     Inputs: ZHI,ZLO: coordinates for start and end of axis
c
c-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
      ZTICK=ABS(ZHI-ZLO)
      ZLOG=LOG10(ZTICK)-1
      IZLOG=ZLOG
      IF (ZLOG.LT.0) IZLOG=IZLOG-1
      gxTICK=10.D0**IZLOG
      NTCKS=ZTICK/gxTICK
      IF (NTCKS.GT.60) THEN
        gxTICK=gxTICK*10
      ELSE IF (NTCKS.GT.30) THEN
        gxTICK=gxTICK*5
      ELSE IF (NTCKS.GT.12) THEN
        gxTICK=gxTICK*2
      END IF
c
      END
c
      SUBROUTINE GINFO(ISOPEN,XRES,YRES,NCOL,CHWID,CHHIGH)
c-----------------------------------------------------------------------
c
c     Purpose: give  informations about properties of graphics screen
c
c     Outputs:
c      ISOPEN: ISOPEN = 1 : graphics system opened
c              ISOPEN = 0 : graphics system not yet opened
c     NCOL:    number of foreground colours
c     XRES:    resolution in x direction (screen coordinates)
c     YRES:    resolution in y direction (screen coordinates)
c     CHWID:   width of one character
c              (screen coordinates, including interval space)
c     CHHIGH:  height of one character
c             (screen coordinates, including interval space)
c
c     COMMON/gxINFO/OPEN,MAXCOL,XAUF,YAUF,XCHHI,XCHWID
c     OPEN    = 1 , graphics opened
c             = 0 , graphics not opened
c             < 0 , opening failed
c     MAXCOL  = number of foreground colours
c     XAUF    = resolution in x direction
c     YAUF    = resolution in y direction
c     XCHWID  = character width, including interval space
c     XCHHI   = character height, including interval space
c     COLOR   = colours of the palette, see GINIT
c
c-----------------------------------------------------------------------
      INTEGER ISOPEN,NCOL
      REAL*8 XRES,YRES,CHWID,CHHIGH
c
      INTEGER OPEN,MAXCOL
      REAL*8 XAUF,YAUF,XCHHI,XCHWID
      COMMON/gxINFO/OPEN,MAXCOL,XAUF,YAUF,XCHHI,XCHWID
      CHWID =XCHWID
      CHHIGH=XCHHI
c
      RETURN
      END
c
      SUBROUTINE GLINE(X1,Y1,X2,Y2,ICOLOR,ISTYLE)
c-----------------------------------------------------------------------
c
c     Purpose: draws a straight line between (X1,Y1) and (X2,Y2).
c
c     Inputs:
c     X1:     x value of the first point
c     Y1:     y value of the first point
c     X2:     x value of the second point
c     Y2:     y value of the second point
c
c     ICOLOR: colour index (0-15) , see explanation in GINIT
c
c     ISTYLE: line style
c       ISTYLE=1: solid
c       ISTYLE=2: long dashes
c       ISTYLE=3: dotted
c       ISTYLE=4: chain-dotted (1 dot)
c       ISTYLE=5: short dashes
c       ISTYLE=6: chain-dotted (2 dots)
c-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XY(4)
      INTEGER gxMOD,GTRANS
c
      ISTAT1=GTRANS(X1,Y1,XY(1),XY(2),3)
      ISTAT2=GTRANS(X2,Y2,XY(3),XY(4),3)
      ISTAT=MAX(ISTAT1,ISTAT2)
      IF (ISTAT.NE.0) ISTAT=gxMOD(XY)
      IF (ISTAT.EQ.0) then
        call pppcol(ICOLOR)
        call plotl(XY(1),XY(2),0.0d0,3)
        call plotl(XY(3),XY(4),0.0d0,2)
      end if
      END
c
      INTEGER FUNCTION gxMOD(XY)
c-----------------------------------------------------------------------
c
c     Purpose: shortens the line (XY(1),XY(2)) - (XY(3),XY(4)), 
c              so that it fits entirely within the window.
c
c-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XY(4)
c
      DOUBLE PRECISION XLO,XHI,YLO,YHI,XWLO,XWHI,YWLO,YWHI
      COMMON/gxWIND/XLO,XHI,YLO,YHI,XWLO,XWHI,YWLO,YWHI
      gxMOD=1

c.... Shorten the line so that it fits within the x margins of the window
      DX=XY(3)-XY(1)
      DY=XY(4)-XY(2)
      IF ( DX.EQ.0 ) THEN
        IF ( (XY(1).GT.XWHI) .OR. (XY(1).LT.XWLO) ) RETURN
      ELSE
        AL1=(XWLO-XY(1))/DX
        AL2=(XWHI-XY(1))/DX
        ALO=MAX(AL1,AL2)
        IF (ALO.LT.0) RETURN
        IF (ALO.LT.1) THEN
          XY(4)=XY(2)+ALO*DY
          XY(3)=XY(1)+ALO*DX
        END IF
        ALU=MIN(AL1,AL2)
        IF (ALU.GT.1) RETURN
        IF (ALU.GT.0) THEN
          XY(2)=XY(2)+ALU*DY
          XY(1)=XY(1)+ALU*DX
        END IF
      END IF

c.... Shorten the line so that it fits within the y margins of the window
      DX=XY(3)-XY(1)
      DY=XY(4)-XY(2)
      IF ( DY.EQ.0 ) THEN
       IF ( (XY(2).GT.YWHI) .OR. (XY(2).LT.YWLO) ) RETURN
      ELSE
        AL1=(YWLO-XY(2))/DY
        AL2=(YWHI-XY(2))/DY
        ALO=MAX(AL1,AL2)
        IF (ALO.LT.0) RETURN
        IF (ALO.LT.1) THEN
          XY(4)=XY(2)+ALO*DY
          XY(3)=XY(1)+ALO*DX
        END IF
        ALU=MIN(AL1,AL2)
        IF (ALU.GT.1) RETURN
        IF (ALU.GT.0) THEN
          XY(2)=XY(2)+ALU*DY
          XY(1)=XY(1)+ALU*DX
        END IF
      END IF
c
      gxMOD=0
c
      END
c
      SUBROUTINE GWRITE(IAREA,XREL,YREL,NOTE,ICOLOR)
c-----------------------------------------------------------------------
c
c     Purpose: displays a row of text on the screen.
c
c     Inputs:
c       IAREA:   coordinate system for XREL and YREL
c       IAREA=1: screen coordinates (0-1)
c       IAREA=2: window coordinates (0-1)
c       IAREA=3: world coordinates
c       XREL:    x coordinate of the left edge of the text
c       YREL:    y coordinate of the baseline of the text
c       NOTE:    text to be output
c       ICOLOR:  colour index (0-15), see explanation in GINIT
c
c-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*(*) NOTE
      CHARACTER*1 NOTE1(80)
      DIMENSION XY(2)
      INTEGER GTRANS
c
      IF (IAREA.GT.0.AND.IAREA.LE.3) THEN
        ISTAT=GTRANS(XREL,YREL,XY(1),XY(2),IAREA)
        nstrng=len(note)
        call dplot(1.25d0*XY(1),1.25d0*XY(2),3,0)
        do 100 i=1,nstrng
100     note1(i)=note(i:i)
        call tplot(NOTE1,nstrng)
      END IF
      END
c
      INTEGER FUNCTION GTRANS(XIN,YIN,XOUT,YOUT,IDIR)
c-----------------------------------------------------------------------
c
c     Purpose:  GTRANS carries out a coordinate transformation.
c
c     Inputs:
c     XIN :  old x coordinate
c     YIN :  old x coordinate
c     IDIR : direction of transformation
c            IDIR = -3 : screen coordinates -> world coordinates
c            IDIR = -2 : screen coordinates -> window cooordinates
c            IDIR = -1, 0, 1 : identity transformation
c             screen coordinates -> screen coordinates
c            IDIR =  2 : window coordinates -> screen coordinates
c            IDIR =  3 : world coordinates  -> screen coordinates
c
c     Outputs:
c     XOUT :  new x coordinate
c     YOUT :  new y coordinate
c
c     
c     GTRANS : GTRANS = 0 : coordinates lie in permissible range
c              GTRANS = 1 : coordinates lie OUTSIDE permissible range
c              GTRANS = 2 : wrong value for IDIR
c
c-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (EPS=1.D-10)
      LOGICAL BAD
c
      DOUBLE PRECISION XLO,XHI,YLO,YHI,XWLO,XWHI,YWLO,YWHI
      COMMON/gxWIND/XLO,XHI,YLO,YHI,XWLO,XWHI,YWLO,YWHI
*     For the significance of the COMMON quantities, see GWINDO
c
      ATRAN(S,SMIN,SMAX,TMIN,TMAX)=
     &  TMIN+(TMAX-TMIN)*(S-SMIN)/(SMAX-SMIN)
      BAD(X,X0,X1) = (X.LT.X0-EPS) .OR. (X.GT.X1+EPS)
c
      GTRANS=2
      IF (ABS(IDIR).GT.3) RETURN
c
      GTRANS=0
      GOTO (103,102,1,1,1,2,3),IDIR+4
c
1     XOUT=XIN
      YOUT=YIN
      IF (BAD(YOUT,0.d0,1.d0).OR.BAD(XOUT,0.d0,1.d0)) GTRANS=1
      RETURN
c
2     XOUT=ATRAN(XIN,0.d0,1.d0,XWLO,XWHI)
      YOUT=ATRAN(YIN,0.d0,1.d0,YWLO,YWHI)
      IF (BAD(YOUT,YWLO,YWHI).OR.BAD(XOUT,XWLO,XWHI)) GTRANS=1
      RETURN
c
3     XOUT=ATRAN(XIN,XLO,XHI,XWLO,XWHI)
      YOUT=ATRAN(YIN,YLO,YHI,YWLO,YWHI)
      IF (BAD(YOUT,YWLO,YWHI).OR.BAD(XOUT,XWLO,XWHI)) GTRANS=1
      RETURN
c
102   IF (BAD(XIN,XWLO,XWHI).OR.BAD(YIN,YWLO,YWHI)) GTRANS=1
      XOUT=ATRAN(XIN,XWLO,XWHI,0.d0,1.d0)
      YOUT=ATRAN(YIN,YWLO,YWHI,0.d0,1.d0)
      RETURN
c
103   IF (BAD(XIN,XWLO,XWHI).OR.BAD(YIN,YWLO,YWHI)) GTRANS=1
      XOUT=ATRAN(XIN,XWLO,XWHI,XLO,XHI)
      YOUT=ATRAN(YIN,YWLO,YWHI,YLO,YHI)
      RETURN
c
      END
c
      SUBROUTINE GDRAW(NUMBER,X,Y)
c-----------------------------------------------------------------------
c
c     Purpose: draws a straight line to the current point.
c              GDRAW may be used an arbitrary number of times in 
c              succession, in order to draw a curve. To fix the starting 
c              point of the curve, subroutine GMOVE1 must be called 
c              before the first call of GDRAW.
c
c     Inputs:
c
c     NUMBER: identification number of the curve
c       (1 < NUMBER < 20)
c     X:      x value of the current point
c     Y:      y value of the current point
c
c-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c
      PARAMETER (NMOVE=20)
      DOUBLE PRECISION XX(NMOVE),YY(NMOVE)
      INTEGER PCOLOR(NMOVE), PSTYLE(NMOVE)
      COMMON /gxPLOT/ XX, YY, PCOLOR, PSTYLE
c
      IF (NUMBER.GT.0.AND.NUMBER.LE.NMOVE) THEN
       CALL GLINE(XX(NUMBER),YY(NUMBER),X,Y,
     &      PCOLOR(NUMBER),PSTYLE(NUMBER))
       XX(NUMBER)=X
       YY(NUMBER)=Y
      END IF
c
      END
c
      SUBROUTINE GMOVE1(NUMBER,X,Y,ICOLOR,ISTYLE)
c-----------------------------------------------------------------------
c
c     Purpose: defines the starting point for a curve to be drawn. 
c              Up to 20 curves may be drawn simultaneously, each being
c              identified by the parameter NUMBER.
c              The actual drawing of the curves is carried out by
c              subroutine GDRAW.
c
c     Inputs:
c     NUMBER: Identification number of the curve (1 < NUMBER < 20)
c     X:      x value of the starting point
c     Y:      y-value of the starting point
c     ICOLOR: colour index (0-15), see explanation in GINIT
c     ISTYLE: line type, see explanation in GLINE
c
c     COMMON /gxPLOT/ XX, YY, PCOLOR, PSTYLE
c     XX(I) =    current x coordinate for curve no. I
c     YY(I) =    current y coordinate for curve no. I
c     PCOLOR(I) =  current colour index for curve no. I
c     PSTYLE(I) =  current line type  for curve no. I
c
c-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c
      PARAMETER (NMOVE=20)
      DOUBLE PRECISION XX(NMOVE),YY(NMOVE)
      INTEGER PCOLOR(NMOVE), PSTYLE(NMOVE)
      COMMON /gxPLOT/ XX, YY, PCOLOR, PSTYLE
c
      call pppcol(icolor)
      IF (NUMBER.GT.0.AND.NUMBER.LE.NMOVE) THEN
        PSTYLE(NUMBER)=ISTYLE
        PCOLOR(NUMBER)=ICOLOR
        XX(NUMBER)=X
        YY(NUMBER)=Y
      END IF
c
      END
c
      SUBROUTINE GCLOSE
c-----------------------------------------------------------------------
c
c     Purpose: asks the user to press a key and closes the graphics.
c
c-----------------------------------------------------------------------
c
      CALL GINIT(-1)
      END
c
      INTEGER FUNCTION gxSKIP(STR,IDIR,CHA)
c-----------------------------------------------------------------------
c
c     Purpose: searches for the first character different from CHA in 
c              string STR
c       IDIR ò 0 : Search from start of string to end
c       IDIR < 0 : Search from end of string to start
c       function value: position found
c                       If all characters are equal CHA, then
c                       gxSKIP = 1 ,       if IDIR < 0
c                       gxSKIP = LEN(STR), if IDIR >=0
c
c-----------------------------------------------------------------------
c
      CHARACTER*(*) STR
      CHARACTER*1 CHA
      LOGICAL GO
c
      ILEN=LEN(STR)
      GO=.TRUE.
      IF (IDIR.GE.0) THEN
        I1=1
        I2=ILEN
      ELSE
        I1=ILEN
        I2=1
      END IF
      ID=SIGN(1,IDIR)
c
      DO 10 I=I1,I2,ID
        IF (GO) gxSKIP=I
        IF (STR(I:I).NE.CHA) GO=.FALSE.
10    CONTINUE
      RETURN
      END
c
      SUBROUTINE GMARK(X,Y,ICOLOR,ISTYLE)
c-----------------------------------------------------------------------
c
c     Purpose: GMARK draws a mark.
c
c     Inputs:
c     X:      x value of the mark
c     Y:      y value of the mark
c     ICOLOR: colour index (0-15), see explanation in GINIT
c     ISTYLE: type of mark
c       ISTYLE=0: small dot
c       ISTYLE=1: large dot
c       ISTYLE=2: plus sign
c     [ ISTYLE=3: asterisk  ]
c       ISTYLE=3: triangle
c       ISTYLE=4: square
c       ISTYLE=5: cross
c       ISTYLE=6: lozenge
c
c-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c
      COMMON/gxWIND/XLO,XHI,YLO,YHI,XWLO,XWHI,YWLO,YWHI
      INTEGER GTRANS
c.... plot mark if point is in window
      if(xlo.le.x.and.x.le.xhi.and.ylo.le.y.and.y.le.yhi) then
C
      ISTAT=GTRANS(X,Y,xx,yy,3)
      dx1 = .005d0
      dx12= .0025d0
      call pppcol(icolor)
      zz = 0.d0
      if (ISTYLE.eq.0) then
        call plotl(xx-dx12, yy+dx12, zz, 3)
        call plotl(xx-dx12, yy+dx12, zz, 1)
        call plotl(xx-dx12, yy-dx12, zz, 2)
        call plotl(xx+dx12, yy-dx12, zz, 2)
        call plotl(xx+dx12, yy+dx12, zz, 2)
        call plotl(xx-dx12, yy+dx12, zz, 2)
        call clpan
      elseif (ISTYLE.eq.1) then
        call plotl(xx-dx1 , yy+dx1 , zz, 3)
        call plotl(xx-dx1 , yy+dx1 , zz, 1)
        call plotl(xx-dx1 , yy-dx1 , zz, 2)
        call plotl(xx+dx1 , yy-dx1 , zz, 2)
        call plotl(xx+dx1 , yy+dx1 , zz, 2)
        call plotl(xx-dx1 , yy+dx1 , zz, 2)
        call clpan
      elseif (ISTYLE.eq.2) then
        call plotl(xx-dx1 , yy     , zz, 3)
        call plotl(xx+dx1 , yy     , zz, 2)
        call plotl(xx   , yy-dx1 , zz, 3)
        call plotl(xx   , yy+dx1 , zz, 2)
      elseif (ISTYLE.eq.3) then
        call plotl(xx-dx1 , yy-dx1*0.66 , zz, 3)
        call plotl(xx-dx1 , yy-dx1*0.66 , zz, 1)
        call plotl(xx   , yy+dx1*1.33 , zz, 2)
        call plotl(xx+dx1 , yy-dx1*0.66 , zz, 2)
        call plotl(xx-dx1 , yy-dx1*0.66 , zz, 2)
        call clpan
      elseif (ISTYLE.eq.4) then
        call plotl(xx-dx1 , yy+dx1 , zz, 3)
        call plotl(xx-dx1 , yy-dx1 , zz, 2)
        call plotl(xx+dx1 , yy-dx1 , zz, 2)
        call plotl(xx+dx1 , yy+dx1 , zz, 2)
        call plotl(xx-dx1 , yy+dx1 , zz, 2)
      elseif (ISTYLE.eq.5) then
        call plotl(xx-dx1 , yy+dx1 , zz, 3)
        call plotl(xx+dx1 , yy-dx1 , zz, 2)
        call plotl(xx-dx1 , yy-dx1 , zz, 3)
        call plotl(xx+dx1 , yy+dx1 , zz, 2)
      elseif (ISTYLE.eq.6) then
        call plotl(xx   , yy+dx1 , zz, 3)
        call plotl(xx   , yy+dx1 , zz, 1)
        call plotl(xx+dx1 , yy     , zz, 2)
        call plotl(xx   , yy-dx1 , zz, 2)
        call plotl(xx-dx1 , yy     , zz, 2)
        call plotl(xx   , yy+dx1 , zz, 2)
        call clpan
      endif
      else
c.... do not plot mark 
      endif
      end


      subroutine modscal(isw)
c-----------------------------------------------------------------------
c.... modify scale values for TPLO used in dplot,plotl                 |
c     W. Wagner  IBS UKA              12/99                            |
c-----------------------------------------------------------------------
      USE pdata1
      USE pdata2
      USE ppers
      implicit double precision(a-h,o-z)
      logical isos
c.... for correct plotting
      dimension s0s(2),sxs(2),dxs(2)
      save scales,scalegs,s0s,sxs,dxs,facts,isos,kperss

      goto(1,2) isw
c....   save original scale values 
1       scales  = scale
        scalegs = scaleg
        s0s(1)  = s0(1)
        s0s(2)  = s0(2)
        sxs(1)  = sx(1)
        sxs(2)  = sx(2)
        dxs(1)  = dx(1) 
        dxs(2)  = dx(2)
        facts   = fact
        isos    = iso 
        kperss  = kpers
c
c....   set values for correct scaling in tplo
        scale  = .77d0/1.3d0
        scaleg = .77d0/1.3d0
        s0(1)  = 0.d0
        s0(2)  = 0.d0
        sx(1)  = 0.d0
        sx(2)  = 0.d0
        dx(1)  = 0.d0
        dx(2)  = 0.d0
        fact   = 1.d0
c...    iso
        iso    = .false. 
c...    pers
        kpers  = 0
c
        return

c....   set values for correct scaling in tplo
2       scale  = scales       
        scaleg = scalegs    
        s0(1)  = s0s(1)
        s0(2)  = s0s(2)
        sx(1)  = sxs(1)
        sx(2)  = sxs(2)
        dx(1)  = dxs(1)
        dx(2)  = dxs(2)
        fact   = facts 
c...    iso
        iso    = isos 
c...    pers
        kpers  = kperss
c
c....   wipe for next plot
        iclear = 0

        return
      end