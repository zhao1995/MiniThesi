      subroutine matelib3d(h1,h2,nh,d,md,Eps,Sig,Cmat,nsig,ntyp,plout,
     +        xgp,tgp,dvp,detf,skfy,skfz,ngp,lgp,lay1gp,lay2gp,imat,isw)
c-----------------------------------------------------------------------
c
c      Purpose: set values for different 3D material laws
c
c      Inputs:
c         h1(nh)   - history array h1
c         h2(nh)   - history array h2
c         d(md)    - local d-array
c         Eps(nsig)- strains
c         nsig     - number of strain/stress components
c         ntyp     - 1=3D, 2=shell, 3=shell ext., 4=beam, 5=beam ext.
c         xgp(3)   - coordinates of Gauss point
c         tgp      - temperature load at Gauss point
c         dvp      - det J*W at Gauss point
c         detf     - det F at Gauss point
c         skfy     - sqrt of shear correct. factor k_y or kappa if only 1
c         skfz     - sqrt of shear correct. factor k_z
c         ngp      - element number
c         lgp      - Gauss point number
c         lay1gp   - layer number
c         lay2gp   - Gauss point number at layer
c         imat     - material type
c         isw      - solution option from element
c
c      Outputs:
c         md       - number of used data for control of d-array
c         nh       - length of history array at Gauss-Point
c         Sig(nsig)- stresses
c         Cmat(nsig,nsig)- material matrix
c         plout(10)- plot data
c
c      Scratch:
c
c-----------------------------------------------------------------------
c
c      Implemented materials
c        1 = linear elastic isotropic
c        2 = linear elastic orthotropic
c        3 = linear elastic transversal isotropic
c        4 = small elasto-(visco-)plastic strains isotropic E = Ee + Ep
c        5 = finite elastoplastic strains F = Fe Fp
c        6 = finite elastic strains, Ogden
c        7 = embedded cell approach Aboudi(Gardner)
c        8 = FE^2 Test
c        9 = small strain isotropic damage
c       10 = concrete model Schütt
c       11 = small strain visco-elastic
c       12 = linear piezoelectric
c       13 = future: dielectric elastomers
c       14 = finite elastic strains, Blatz-Ko
c       15 = linear elastic transversal isotropic with damage
c
c        material formulations assume the following order of components
c        11,22,33,12,13,23
c
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension h1(*),h2(*),d(*),xgp(3)
      dimension Eps(nsig),Sig(nsig),Cmat(nsig,nsig),plout(10)
c
      sig   = 0
      cmat  = 0
      plout = 0


c...  material models
c
      if(imat.eq.1) then
c       linear elastic isotropic
        call mate3d01(h1,h2,nh,d,md,Eps,Sig,Cmat,nsig,ntyp,plout,
     +             xgp,tgp,dvp,detf,skfy,skfz,ngp,lgp,lay1gp,lay2gp,isw)
c
      else if(imat.eq.2) then
c       linear elastic orthotropic
        call mate3d02(h1,h2,nh,d,md,Eps,Sig,Cmat,nsig,ntyp,plout,
     +             xgp,tgp,dvp,detf,skfy,skfz,ngp,lgp,lay1gp,lay2gp,isw)
c
      else if(imat.eq.3) then
c       linear elastic transversal isotropic
        call mate3d03(h1,h2,nh,d,md,Eps,Sig,Cmat,nsig,ntyp,plout,
     +             xgp,tgp,dvp,detf,skfy,skfz,ngp,lgp,lay1gp,lay2gp,isw)
c
      else if(imat.eq.4) then
c       small elasto-(visco-)plastic strains E = Ee + Ep
        call mate3d04(h1,h2,nh,d,md,Eps,Sig,Cmat,nsig,ntyp,plout,
     +             xgp,tgp,dvp,detf,skfy,skfz,ngp,lgp,lay1gp,lay2gp,isw)
c
      else if(imat.eq.5) then
c       finite elastoplastic strains F = Fe * Fp
        call mate3d05(h1,h2,nh,d,md,Eps,Sig,Cmat,nsig,ntyp,plout,
     +             xgp,tgp,dvp,detf,skfy,skfz,ngp,lgp,lay1gp,lay2gp,isw)
c
      else if(imat.eq.6) then
c       finite elastic strains, Ogden
        call mate3d06(h1,h2,nh,d,md,Eps,Sig,Cmat,nsig,ntyp,plout,
     +             xgp,tgp,dvp,detf,skfy,skfz,ngp,lgp,lay1gp,lay2gp,isw)
c
      else if(imat.eq.7) then
c       embedded cell approach Aboudi(Gardner)
        call mate3d07(h1,h2,nh,d,md,Eps,Sig,Cmat,nsig,ntyp,plout,
     +             xgp,tgp,dvp,detf,skfy,skfz,ngp,lgp,lay1gp,lay2gp,isw)
c
      else if(imat.eq.8) then
c       FE^2
        call mate3d08(h1,h2,nh,d,md,Eps,Sig,Cmat,nsig,ntyp,plout,
     +             xgp,tgp,dvp,detf,skfy,skfz,ngp,lgp,lay1gp,lay2gp,isw)
c
      else if(imat.eq.9) then
c       small strain isotropic damage
        call mate3d09(h1,h2,nh,d,md,Eps,Sig,Cmat,nsig,ntyp,plout,
     +             xgp,tgp,dvp,detf,skfy,skfz,ngp,lgp,lay1gp,lay2gp,isw)
c
      else if(imat.eq.10) then
c       concrete model Schütt
        call mate3d10(h1,h2,nh,d,md,Eps,Sig,Cmat,nsig,ntyp,plout,
     +             xgp,tgp,dvp,detf,skfy,skfz,ngp,lgp,lay1gp,lay2gp,isw)
c
      else if(imat.eq.11) then
c       small strain visco-elastic
        call mate3d11(h1,h2,nh,d,md,Eps,Sig,Cmat,nsig,ntyp,plout,
     +             xgp,tgp,dvp,detf,skfy,skfz,ngp,lgp,lay1gp,lay2gp,isw)
c
      else if(imat.eq.12) then
c       linear piezoelectric
        call mate3d12(h1,h2,nh,d,md,Eps,Sig,Cmat,nsig,ntyp,plout,
     +             xgp,tgp,dvp,detf,skfy,skfz,ngp,lgp,lay1gp,lay2gp,isw)
c
      else if(imat.eq.13) then
c       dielectric elastomers
c       call mate3d13(h1,h2,nh,d,md,Eps,Sig,Cmat,nsig,ntyp,plout,
c    +             xgp,tgp,dvp,detf,skfy,skfz,ngp,lgp,lay1gp,lay2gp,isw)
c
c
      else if(imat.eq.14) then
c       finite elastic strains, Blatz-Ko
        call mate3d14(h1,h2,nh,d,md,Eps,Sig,Cmat,nsig,ntyp,plout,
     +             xgp,tgp,dvp,detf,skfy,skfz,ngp,lgp,lay1gp,lay2gp,isw)

      else if(imat.eq.15) then
c       linear elastic transversal isotropic with damage
        call mate3d15(h1,h2,nh,d,md,Eps,Sig,Cmat,nsig,ntyp,plout,
     +             xgp,tgp,dvp,detf,skfy,skfz,ngp,lgp,lay1gp,lay2gp,isw)

      else if(imat.eq.16) then
c       linear elastic transversal isotropic with damage
        call mate3d16(h1,h2,nh,d,md,Eps,Sig,Cmat,nsig,ntyp,plout,
     +             xgp,tgp,dvp,detf,skfy,skfz,ngp,lgp,lay1gp,lay2gp,isw)

      else if(imat.eq.17) then
c       linear elastic transversal isotropic with damage
        call mate3d17(h1,h2,nh,d,md,Eps,Sig,Cmat,nsig,ntyp,plout,
     +             xgp,tgp,dvp,detf,skfy,skfz,ngp,lgp,lay1gp,lay2gp,isw)

      end if
c
      return
      end
c
