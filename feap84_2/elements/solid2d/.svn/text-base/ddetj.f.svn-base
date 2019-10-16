c$Id:$
      function   ddetj(j,theta)

c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Purpose: Series expansion for
c              ((1+theta)**1/3 - (1+j)**1/3)/(1+j)**1/3

c     Input:
c       j      - Value to expand

c     Output:
c       ddetj  - ((1+theta)**1/3 - (1+j)**1/3)/(1+j)**1/3
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      real*8     ddetj, j, theta
      real*8     c1,c2,c3,c4,c5,c6

      save

      data       c1 / 0.33333333333333331d0 /
      data       c2 /-0.11111111111111110d0 /
      data       c3 / 6.17283950617283916d-02 /
      data       c4 /-4.11522633744855967d-02 /
      data       c5 / 3.01783264746227700d-02 /
      data       c6 /-2.17954580094497780d-02 /

c     Compute difference by series expansion

      ddetj = c1*(theta    - j   )
     &      + c2*(theta**2 - j**2)
     &      + c3*(theta**3 - j**3)
     &      + c4*(theta**4 - j**4)
     &      + c5*(theta**5 - j**5)
     &      + c6*(theta**6 - j**6)

      end
