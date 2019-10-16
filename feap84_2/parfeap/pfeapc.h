
      Vec             rhs, sol, xvec
      common /pfeapc/ rhs, sol, xvec

      Vec             yvec, zvec
      common /pfeapc/ yvec, zvec

      Vec             Mdiag, Msqrt
      common /pfeapc/ Mdiag, Msqrt

      Mat             Kmat, Mmat, Pmat
      common /pfeapc/ Kmat, Mmat, Pmat

      KSP             kspsol
      common /pfeapc/ kspsol
