      subroutine elmlib(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,iel,isw)
c----------------------------------------------------------------------
c      Purpose: Element library driver routine

c      Inputs:
c         d(*)    - Material parameters
c         u(*)    - Element solution parameters
c         x(*)    - Element nodal coordinates
c         ix(*)   - Element nodal numbers
c         t(*)    - Element temperatures
c         h1(*)   - History array h1(nhmax)
c         h2(*)   - History array h2(nhmax)
c         h3(*)   - History array h3(nh3max)
c         i       - Number dof/node           (ndf)
c         j       - Spatial dimension of mesh (ndm)
c         k       - Size of element arrays    (nst)
c         iel     - Element type number
c         isw     - Switch    
c
c      Outputs:
c         d(*)    - Material parameters (isw = 1 only)
c         s(*,*)  - Element array
c         p(*)    - Element vector
c         h1(*)   - History array h1(nhmax)
c         h2(*)   - History array h2(nhmax)
c         h3(*)   - History array h3(nh3max)
c----------------------------------------------------------------------
      USE iofile
      implicit double precision (a-h,o-z)
      dimension p(k),s(k,k),d(*),u(*),x(*),ix(*),t(*),h1(*),h2(*),h3(*)
      if(iel.le.0.or.iel.gt.100) go to 400
      if(isw.ge.3) then 
        call pzero(s,k*k)
        call pzero(p,k)
      end if
      go to ( 1, 2, 3, 4, 5, 6, 7, 8, 9,10,
     1	      11,12,13,14,15,16,17,18,19,20,
     2	      21,22,23,24,25,26,27,28,29,30,
     3	      31,32,33,34,35,36,37,38,39,40,
     4	      41,42,43,44,45,46,47,48,49,50,
     5	      51,52,53,54,55,56,57,58,59,60,
     6	      61,62,63,64,65,66,67,68,69,70,
     7	      71,72,73,74,75,76,77,78,79,80,
     8	      81,82,83,84,85,86,87,88,89,90,
     9	      91,92,93,94,95,96,97,98,99,100),iel

1     call elmt01(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200
2     call elmt02(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200
3     call elmt03(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200
4     call elmt04(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200
5     call elmt05(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200
6     call elmt06(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200
7     call elmt07(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200
8     call elmt08(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200
9     call elmt09(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200
10    call elmt10(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200
11    call elmt11(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200
12    call elmt12(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200
13    call elmt13(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200
14    call elmt14(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200
15    call elmt15(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200
16    call elmt16(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
17    call elmt17(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
18    call elmt18(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
19    call elmt19(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
20    call elmt20(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
21    call elmt21(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
22    call elmt22(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
23    call elmt23(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
24    call elmt24(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
25    call elmt25(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
26    call elmt26(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
27    call elmt27(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
28    call elmt28(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
29    call elmt29(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
30    call elmt30(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200
31    call elmt31(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
32    call elmt32(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
33    call elmt33(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
34    call elmt34(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
35    call elmt35(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
36    call elmt36(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
37    call elmt37(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
38    call elmt38(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
39    call elmt39(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
40    call elmt40(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
41    call elmt41(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
42    call elmt42(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
43    call elmt43(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
44    call elmt44(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
45    call elmt45(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200
46    call elmt46(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
47    call elmt47(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
48    call elmt48(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
49    call elmt49(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
50    call elmt50(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
51    call elmt51(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
52    call elmt52(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
53    call elmt53(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
54    call elmt54(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
55    call elmt55(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
56    call elmt56(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
57    call elmt57(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
58    call elmt58(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
59    call elmt59(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
60    call elmt60(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200
61    call elmt61(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
62    call elmt62(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
63    call elmt63(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
64    call elmt64(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
65    call elmt65(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
66    call elmt66(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
67    call elmt67(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
68    call elmt68(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
69    call elmt69(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
70    call elmt70(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
71    call elmt71(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
72    call elmt72(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
73    call elmt73(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
74    call elmt74(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
75    call elmt75(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200
76    call elmt76(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
77    call elmt77(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
78    call elmt78(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
79    call elmt79(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
80    call elmt80(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
81    call elmt81(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
82    call elmt82(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
83    call elmt83(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
84    call elmt84(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
85    call elmt85(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
86    call elmt86(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
87    call elmt87(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
88    call elmt88(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
89    call elmt89(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
90    call elmt90(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200
91    call elmt91(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
92    call elmt92(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
93    call elmt93(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
94    call elmt94(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
95    call elmt95(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
96    call elmt96(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
97    call elmt97(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
98    call elmt98(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
99    call elmt99(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
100   call elmt100(d,u,x,ix,t,s,p,h1,h2,h3,i,j,k,isw)
      go to 200                           
200   return                      
400   if(ior.gt.0) write(iow,4000) iel
      if(ior.lt.0) write(  *,4000) iel
      stop                                
4000  format('  **error** elementh1,h2,h3, type number',i3,' input')
      end                                 
c                                 
c      subroutine elmt01                  
c      write(*,2000)              
c2000  format(1x,'**** Warning dummy subroutine elmt01 called ****')
c      end
c
c      subroutine elmt02
c      write(*,2000)
c2000  format(1x,'**** Warning dummy subroutine elmt02 called ****')
c      end
c
c      subroutine elmt03
c      write(*,2000)
c2000  format(1x,'**** Warning dummy subroutine elmt03 called ****')
c      end
c
c      subroutine elmt04
c      write(*,2000)
c2000  format(1x,'**** Warning dummy subroutine elmt04 called ****')
c      end
c
c      subroutine elmt05
c      write(*,2000)
c2000  format(1x,'**** Warning dummy subroutine elmt05 called ****')
c      end
c
c      subroutine elmt06
c      write(*,2000)
c2000  format(1x,'**** Warning dummy subroutine elmt06 called ****')
c      end
c
c      subroutine elmt07
c      write(*,2000)
c2000  format(1x,'**** Warning dummy subroutine elmt07 called ****')
c      end
c
c      subroutine elmt08
c      write(*,2000)
c2000  format(1x,'**** Warning dummy subroutine elmt08 called ****')
c      end
c
c      subroutine elmt09
c      write(*,2000)
c2000  format(1x,'**** Warning dummy subroutine elmt09 called ****')
c      end
c
c      subroutine elmt10
c      write(*,2000)
c2000  format(1x,'**** Warning dummy subroutine elmt10 called ****')
c      end
c
c      subroutine elmt11
c      write(*,2000)
c2000  format(1x,'**** Warning dummy subroutine elmt11 called ****')
c      end
c
c      subroutine elmt12
c      write(*,2000)
c2000  format(1x,'**** Warning dummy subroutine elmt12 called ****')
c      end
c
c      subroutine elmt13
c      write(*,2000)
c2000  format(1x,'**** Warning dummy subroutine elmt13 called ****')
c      end
c
c      subroutine elmt14
c      write(*,2000)
c2000  format(1x,'**** Warning dummy subroutine elmt14 called ****')
c      end
c
c      subroutine elmt15
c      write(*,2000)
c2000  format(1x,'**** Warning dummy subroutine elmt15 called ****')
c      end
c
c      subroutine elmt16
c      write(*,2000)
c2000  format(1x,'**** Warning dummy subroutine elmt16 called ****')
c      end
c
c      subroutine elmt17
c      write(*,2000)
c2000  format(1x,'**** Warning dummy subroutine elmt17 called ****')
c      end
c
c      subroutine elmt18
c      write(*,2000)
c2000  format(1x,'**** Warning dummy subroutine elmt18 called ****')
c      end
c
c      subroutine elmt19
c      write(*,2000)
c2000  format(1x,'**** Warning dummy subroutine elmt19 called ****')
c      end
c
c      subroutine elmt20
c      write(*,2000)
c2000  format(1x,'**** Warning dummy subroutine elmt20 called ****')
c      end
c
c      subroutine elmt21
c      write(*,2000)
c2000  format(1x,'**** Warning dummy subroutine elmt21 called ****')
c      end
c
c      subroutine elmt22
c      write(*,2000)
c2000  format(1x,'**** Warning dummy subroutine elmt22 called ****')
c      end
c
c      subroutine elmt23
c      write(*,2000)
c2000  format(1x,'**** Warning dummy subroutine elmt23 called ****')
c      end
c
c      subroutine elmt24
c      write(*,2000)
c2000  format(1x,'**** Warning dummy subroutine elmt24 called ****')
c      end
c
c      subroutine elmt25
c      write(*,2000)
c2000  format(1x,'**** Warning dummy subroutine elmt25 called ****')
c      end
c
c      subroutine elmt26
c      write(*,2000)
c2000  format(1x,'**** Warning dummy subroutine elmt26 called ****')
c      end
c
      subroutine elmt27
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt27 called ****')
      end
c
      subroutine elmt28
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt28 called ****')
      end
c
c      subroutine elmt29
c      write(*,2000)
c2000  format(1x,'**** Warning dummy subroutine elmt29 called ****')
c      end
c
c      subroutine elmt30
c      write(*,2000)
c2000  format(1x,'**** Warning dummy subroutine elmt30 called ****')
c      end
c
c      subroutine elmt31
c      write(*,2000)
c2000  format(1x,'**** Warning dummy subroutine elmt31 called ****')
c      end
c
c      subroutine elmt32
c      write(*,2000)
c2000  format(1x,'**** Warning dummy subroutine elmt32 called ****')
c      end
c
c      subroutine elmt33
c      write(*,2000)
c2000  format(1x,'**** Warning dummy subroutine elmt33 called ****')
c      end
c
c      subroutine elmt34
c      write(*,2000)
c2000  format(1x,'**** Warning dummy subroutine elmt34 called ****')
c      end
c
      subroutine elmt35
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt35 called ****')
      end
c
      subroutine elmt36
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt36 called ****')
      end
c
      subroutine elmt37
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt37 called ****')
      end
c
      subroutine elmt38
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt38 called ****')
      end
c
      subroutine elmt39
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt39 called ****')
      end
c
      subroutine elmt40
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt40 called ****')
      end
c
      subroutine elmt41
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt41 called ****')
      end
c
      subroutine elmt42
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt42 called ****')
      end
c
      subroutine elmt43
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt43 called ****')
      end
c
      subroutine elmt44
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt44 called ****')
      end
c
c      subroutine elmt45
c      write(*,2000)
c2000  format(1x,'**** Warning dummy subroutine elmt45 called ****')
c      end
c
      subroutine elmt46
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt46 called ****')
      end
c
      subroutine elmt47
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt47 called ****')
      end
c
      subroutine elmt48
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt48 called ****')
      end
c
      subroutine elmt49
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt49 called ****')
      end
c
      subroutine elmt50
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt50 called ****')
      end
c
      subroutine elmt51
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt51 called ****')
      end
c
      subroutine elmt52
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt52 called ****')
      end
c
      subroutine elmt53
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt53 called ****')
      end
c
      subroutine elmt54
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt54 called ****')
      end
c
      subroutine elmt55
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt55 called ****')
      end
c
      subroutine elmt56
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt56 called ****')
      end
c
      subroutine elmt57
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt57 called ****')
      end
c
      subroutine elmt58
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt58 called ****')
      end
c
      subroutine elmt59
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt59 called ****')
      end
c
      subroutine elmt60
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt60 called ****')
      end
c
c      subroutine elmt61
c      write(*,2000)
c2000  format(1x,'**** Warning dummy subroutine elmt61 called ****')
c      end
c
c      subroutine elmt62
c      write(*,2000)
c2000  format(1x,'**** Warning dummy subroutine elmt62 called ****')
c      end
c
c      subroutine elmt63
c      write(*,2000)
c2000  format(1x,'**** Warning dummy subroutine elmt63 called ****')
c      end
c
      subroutine elmt64
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt64 called ****')
      end
c
      subroutine elmt65
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt65 called ****')
      end
c
      subroutine elmt66
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt66 called ****')
      end
c
      subroutine elmt67
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt67 called ****')
      end
c
      subroutine elmt68
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt68 called ****')
      end
c
      subroutine elmt69
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt69 called ****')
      end
c
      subroutine elmt70
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt70 called ****')
      end
c
      subroutine elmt71
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt71 called ****')
      end
c
      subroutine elmt72
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt72 called ****')
      end
c
      subroutine elmt73
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt73 called ****')
      end
c
      subroutine elmt74
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt74 called ****')
      end
c
      subroutine elmt75
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt75 called ****')
      end
c
      subroutine elmt76
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt76 called ****')
      end
c
      subroutine elmt77
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt77 called ****')
      end
c
      subroutine elmt78
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt78 called ****')
      end
c
      subroutine elmt79
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt79 called ****')
      end
c
      subroutine elmt80
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt80 called ****')
      end
c
      subroutine elmt81
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt81 called ****')
      end
c
      subroutine elmt82
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt82 called ****')
      end
c
      subroutine elmt83
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt83 called ****')
      end
c
      subroutine elmt84
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt84 called ****')
      end
c
      subroutine elmt85
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt85 called ****')
      end
c
      subroutine elmt86
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt86 called ****')
      end
c
      subroutine elmt87
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt87 called ****')
      end
c
      subroutine elmt88
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt88 called ****')
      end
c
      subroutine elmt89
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt89 called ****')
      end
c
      subroutine elmt90
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt90 called ****')
      end
c
      subroutine elmt91
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt91 called ****')
      end
c
      subroutine elmt92
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt92 called ****')
      end
c
      subroutine elmt93
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt93 called ****')
      end
c
      subroutine elmt94
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt94 called ****')
      end
c
      subroutine elmt95
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt95 called ****')
      end
c
      subroutine elmt96
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt96 called ****')
      end
c
      subroutine elmt97
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt97 called ****')
      end
c
      subroutine elmt98
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt98 called ****')
      end
c
      subroutine elmt99
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt99 called ****')
      end
c
      subroutine elmt100
      write(*,2000)
2000  format(1x,'**** Warning dummy subroutine elmt100 called ****')
      end
c
