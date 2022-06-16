c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2021 by Moses Chung, Zhi Wang & Jay Ponder  ##
c     ##                    All Rights Reserved                     ##
c     ################################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine expolar1  --  ExchPol variable polarizability  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "expolar1" calculates the variable polarizability force components
c     due to exchange polarization
c
c     literature reference:
c
c
      subroutine expolar1
      use atoms
      use bound
      use cell
      use chgpot
      use couple
      use deriv
      use expol
      use mpole
      use polar
      use polgrp
      use polpot
      use shunt
      use units
      implicit none
      integer i,j,k,m
      integer ii,kk
      real*8 xr,yr,zr
      real*8 r,r2
      real*8 f
      real*8 sizi,sizk,sizik
      real*8 alphai,alphak
      real*8 springi,springk
      real*8 s2,ds2
      real*8 s2i, s2k
      real*8 ds2i, ds2k
      real*8 uixl,ukxl
      real*8 uiyl,ukyl
      real*8 uizl,ukzl
      real*8 frcxil,frcxkl
      real*8 frcyil,frcykl
      real*8 frczil,frczkl
      real*8 frctxi,frctxk
      real*8 frctyi,frctyk
      real*8 frctzi,frctzk
      real*8 frcil(3)
      real*8 frckl(3)
      real*8 frcxi,frcyi,frczi
      real*8 frcxk,frcyk,frczk
      real*8 tqxil,tqyil
      real*8 tqxkl,tqykl
      real*8 ai(3,3)
      real*8 ak(3,3)
      real*8, allocatable :: pscale(:)
      character*6 mode
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
c
c     set conversion factor, cutoff and switching coefficients
c
      f = 0.5d0 * electric / dielec
      mode = 'MPOLE'
      call switch (mode)
c
c     set array needed to scale atom and group interactions
c
      do i = 1, n
         pscale(i) = 1.0d0
      end do
c
c     find the exchange polarization gradient
c
      do ii = 1, npole-1
         i = ipole(ii)
         springi = kpep(ii)/polarity(ii)
         sizi = prepep(ii)
         alphai = dmppep(ii)
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            pscale(i12(j,i)) = p2scale
            do k = 1, np11(i)
               if (i12(j,i) .eq. ip11(k,i))
     &               pscale(i12(j,i)) = p2iscale
            end do
         end do
         do j = 1, n13(i)
            pscale(i13(j,i)) = p3scale
            do k = 1, np11(i)
               if (i13(j,i) .eq. ip11(k,i))
     &               pscale(i13(j,i)) = p3iscale
            end do
         end do
         do j = 1, n14(i)
            pscale(i14(j,i)) = p4scale
            do k = 1, np11(i)
               if (i14(j,i) .eq. ip11(k,i))
     &               pscale(i14(j,i)) = p4iscale
            end do
         end do
         do j = 1, n15(i)
            pscale(i15(j,i)) = p5scale
            do k = 1, np11(i)
               if (i15(j,i) .eq. ip11(k,i))
     &               pscale(i15(j,i)) = p5iscale
            end do
         end do
c
c     evaluate all sites within the cutoff distance
c
         do kk = ii+1, npole
            k = ipole(kk)
            xr = x(k) - x(i)
            yr = y(k) - y(i)
            zr = z(k) - z(i)
            if (use_bounds)  call image (xr,yr,zr)
            r2 = xr*xr + yr* yr + zr*zr
            if (r2 .le. off2) then
               r = sqrt(r2)
               springk = kpep(kk)/polarity(kk)
               sizk = prepep(kk)
               alphak = dmppep(kk)
               sizik = sizi*sizk
               call dampep (r,sizik,alphai,alphak,s2,ds2)
               s2i = springi*s2*pscale(k)
               s2k = springk*s2*pscale(k)
               ds2i = springi*ds2*pscale(k)
               ds2k = springk*ds2*pscale(k)
               call exrotmat(i,k,ai,ak)
               uixl = 0.0d0
               ukxl = 0.0d0
               uiyl = 0.0d0
               ukyl = 0.0d0
               uizl = 0.0d0
               ukzl = 0.0d0
               do j = 1, 3
                  uixl = uixl + uind(j,ii)*ai(1,j)
                  ukxl = ukxl - uind(j,kk)*ak(1,j)
                  uiyl = uiyl + uind(j,ii)*ai(2,j)
                  ukyl = ukyl - uind(j,kk)*ak(2,j)
                  uizl = uizl + uind(j,ii)*ai(3,j)
                  ukzl = ukzl - uind(j,kk)*ak(3,j)
               end do
              frcil(3) = uizl**2 * ds2i
              frckl(3) = ukzl**2 * ds2k
c
c     compute torque in local frame
c
               tqxil =  2.0d0 * uiyl * uizl * s2i
               tqyil = -2.0d0 * uixl * uizl * s2i
               tqxkl =  2.0d0 * ukyl * ukzl * s2k
               tqykl = -2.0d0 * ukxl * ukzl * s2k
c
c     convert torque to forces
c
               frcil(1) = -tqyil / r
               frcil(2) = tqxil / r
               frckl(1) = -tqykl / r
               frckl(2) = tqxkl / r
c
c     rotate force to global frame
c
               frcxi = f*ai(3,1)*frcil(3)
               frcyi = f*ai(3,2)*frcil(3)
               frczi = f*ai(3,3)*frcil(3)
               frcxk = f*ak(3,1)*frckl(3)
               frcyk = f*ak(3,2)*frckl(3)
               frczk = f*ak(3,3)*frckl(3)

               frctxi = f*(ai(1,1)*frcil(1)+ai(2,1)*frcil(2))
               frctyi = f*(ai(1,2)*frcil(1)+ai(2,2)*frcil(2))
               frctzi = f*(ai(1,3)*frcil(1)+ai(2,3)*frcil(2))
               frctxk = f*(ak(1,1)*frckl(1)+ak(2,1)*frckl(2))
               frctyk = f*(ak(1,2)*frckl(1)+ak(2,2)*frckl(2))
               frctzk = f*(ak(1,3)*frckl(1)+ak(2,3)*frckl(2))
c
c     increment force-based gradient on the interaction sites
c
               dep(1,i) = dep(1,i) - frcxi + frcxk - frctxi + frctxk
               dep(2,i) = dep(2,i) - frcyi + frcyk - frctyi + frctyk
               dep(3,i) = dep(3,i) - frczi + frczk - frctzi + frctzk
               dep(1,k) = dep(1,k) - frcxk + frcxi - frctxk + frctxi
               dep(2,k) = dep(2,k) - frcyk + frcyi - frctyk + frctyi
               dep(3,k) = dep(3,k) - frczk + frczi - frctzk + frctzi
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            pscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            pscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            pscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            pscale(i15(j,i)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (pscale)
      return
      end

c
c
c     #######################################################
c     ##                                                   ##
c     ##  subroutine exrotmat  --  compute rotation matrix ##
c     ##                                                   ##
c     #######################################################
c
c
c     "exrotmat" computes rotation matrix
c
c
      subroutine exrotmat (i,k,ai,ak)
      use atoms
      use math
      use mpole
      use polpot
      implicit none
      integer i,ii,k,kk
      integer j,l,m,o
      real*8 r,dot
      real*8 eps,angle
      real*8 xi,yi,zi
      real*8 xk,yk,zk
      real*8 dx,dy,dz
      real*8 ai(3,3)
      real*8 ak(3,3)
c
c
c     use the identity matrix as the default rotation matrix
c
      ai(1,1) = 1.0d0
      ai(2,1) = 0.0d0
      ai(3,1) = 0.0d0
      ai(1,2) = 0.0d0
      ai(2,2) = 1.0d0
      ai(3,2) = 0.0d0
      ai(1,3) = 0.0d0
      ai(2,3) = 0.0d0
      ai(3,3) = 1.0d0
c
c     get coordinates and frame definition for the multipole site
c
      xi = x(i)
      yi = y(i)
      zi = z(i)
      xk = x(k)
      yk = y(k)
      zk = z(k)
c
c     get the rotation matrix elements
c
      dx = xk - xi
      dy = yk - yi
      dz = zk - zi
      r = sqrt(dx*dx + dy*dy + dz*dz)
      ai(3,1) = dx / r
      ai(3,2) = dy / r
      ai(3,3) = dz / r
      dx = 1.0d0
      dy = 0.0d0
      dz = 0.0d0
      dot = ai(3,1)
      eps = 0.707d0
      if (abs(dot) .gt. eps) then
         dx = 0.0d0
         dy = 1.0d0
         dot = ai(3,2)
      end if
      dx = dx - dot*ai(3,1)
      dy = dy - dot*ai(3,2)
      dz = dz - dot*ai(3,3)
      r = sqrt(dx*dx + dy*dy + dz*dz)
      ai(1,1) = dx / r
      ai(1,2) = dy / r
      ai(1,3) = dz / r
      ai(2,1) = ai(1,3)*ai(3,2) - ai(1,2)*ai(3,3)
      ai(2,2) = ai(1,1)*ai(3,3) - ai(1,3)*ai(3,1)
      ai(2,3) = ai(1,2)*ai(3,1) - ai(1,1)*ai(3,2)
      ak = ai
      ak(2,1) = -ak(2,1)
      ak(2,2) = -ak(2,2)
      ak(2,3) = -ak(2,3)
      ak(3,1) = -ak(3,1)
      ak(3,2) = -ak(3,2)
      ak(3,3) = -ak(3,3)
      return
      end