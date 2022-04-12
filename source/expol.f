c
c
c     #########################################################
c     ##  COPYRIGHT (C) 2021 by Moses Chung & Jay W. Ponder  ##
c     ##                   All Rights Reserved               ##
c     #########################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  module expol  --  exch polarization for current structure   ##        
c     ##                                                              ##
c     ##################################################################
c
c
c     nexpol    total number of exch polarization sites in the system
c     prepep    exchange polarization constant prefactor at each site
c     dmppep    exchange polarization damping alpha at each site
c     scrtyp    type of screening
c
c
      module expol
      implicit none
      integer nexpol
      real*8, allocatable :: prepep(:)
      real*8, allocatable :: dmppep(:)
      character*3 scrtyp
      save
      end