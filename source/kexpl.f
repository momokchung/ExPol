c
c
c     #########################################################
c     ##  COPYRIGHT (C) 2021 by Moses Chung & Jay W. Ponder  ##
c     ##                   All Rights Reserved               ##
c     #########################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  module kexpl  --  exch polarization forcefield parameters   ##
c     ##                                                              ##
c     ##################################################################
c
c
c     pepk       exchange polarization spring constant
c     peppre     exchange polarization constant prefactor
c     pepdmp     exchange polarization damping alpha
c
c
      module kexpl
      implicit none
      real*8, allocatable :: pepk(:)
      real*8, allocatable :: peppre(:)
      real*8, allocatable :: pepdmp(:)
      save
      end