
                             ! CONTACT KEYS FOR HISTORY VECTOR

      integer         p1,p3,lh1,lh3,nset,tlch1,tlch3,eslave,emastr,eoffs
      common /c_keyh/
     &  p1(c_lp1),           ! Pointer for the variables    in  CH1
     &  p3(c_lp3),           ! Pointer for the variables    in  CH3
     &  lh1,                 ! Length of ch1 (dynamic)
     &  lh3,                 ! Length of ch3 (dynamic)
     &  nset,                ! # of history set required for pairs
     &  tlch1,               ! Total length of CH1 & CH2
     &  tlch3,               ! Total length of CH3
     &  eslave,              ! # of Explicit SLAVE  nodes
     &  emastr,              ! # of Explicit MASTeR nodes
     &  eoffs                ! # of Explicit OFFSet

