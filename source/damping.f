c
c
c     ############################################################
c     ##  COPYRIGHT (C) 2018 by Joshua Rackers & Jay W. Ponder  ##
c     ##                   All Rights Reserved                  ## 
c     ############################################################ 
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine dampewald  --  find Ewald damping coefficients  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "dampewald" finds coefficients for error function damping used
c     for Ewald real space interactions
c
c
      subroutine dampewald (rorder,r,r2,scale,dmpe)
      use ewald
      use math
      implicit none
      integer i,niter
      integer rorder
      real*8 r,r2,scale
      real*8 bfac,erfc
      real*8 aesq2,afac
      real*8 expterm,ra
      real*8 bn(0:5)
      real*8 dmpe(*)
      external erfc
c
c
c     initialize the Ewald damping factor coefficients
c
      do i = 1, rorder
         dmpe(i) = scale
      end do
c     
c     compute the successive Ewald damping factors
c
      ra = aewald * r
      bn(0) = erfc(ra) / r
      dmpe(1) = scale * bn(0)
      expterm = exp(-ra*ra)
      aesq2 = 2.0d0 * aewald * aewald
      afac = 0.0d0
      if (aewald .gt. 0.0d0)  afac = 1.0d0 / (rootpi*aewald)
      niter = (rorder-1) / 2
      do i = 1, niter
         bfac = dble(2*i-1)
         afac = aesq2 * afac
         bn(i) = (bfac*bn(i-1)+afac*expterm) / r2
         dmpe(2*i+1) = scale * bn(i)
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine dampthole  --  find Thole damping coefficients  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "dampthole" finds coefficients for the original Thole damping
c     function used by AMOEBA or for the alternate direct polarization
c     damping used by AMOEBA+
c
c     literature reference:
c
c     B. T. Thole, "Molecular Polarizabilities Calculated with a
c     Modified Dipole Interaction", Chemical Physics, 59, 341-350 (1981)
c
c
      subroutine dampthole (i,k,rorder,r,dmpik)
      use polar
      use polpot
      implicit none
      integer i,j,k
      integer rorder
      real*8 r,damp
      real*8 damp2
      real*8 damp3
      real*8 expdamp
      real*8 pgamma
      real*8 dmpik(*)
c
c
c     initialize the Thole damping factors to a value of one
c
      do j = 1, rorder
         dmpik(j) = 1.0d0
      end do
c
c     use alternate Thole model for AMOEBA+ direct polarization
c
      damp = pdamp(i) * pdamp(k)
      if (use_dirdamp) then
         pgamma = min(dirdamp(i),dirdamp(k))
         if (pgamma .eq. 0.0d0)  pgamma = max(dirdamp(i),dirdamp(k))
         if (damp.ne.0.0d0 .and. pgamma.ne.0.0d0) then
            damp = pgamma * (r/damp)**(1.5d0)
            if (damp .lt. 50.0d0) then
               expdamp = exp(-damp)
               dmpik(3) = 1.0d0 - expdamp
               dmpik(5) = 1.0d0 - expdamp*(1.0d0+0.5d0*damp)
               if (rorder .ge. 7) then
                  damp2 = damp * damp
                  dmpik(7) = 1.0d0 - expdamp*(1.0d0+0.65d0*damp
     &                                  +0.15d0*damp2)
               end if
            end if
         end if
c
c     use original AMOEBA Thole polarization damping factors
c
      else
         pgamma = min(thole(i),thole(k))
         if (pgamma .eq. 0.0d0)  pgamma = max(thole(i),thole(k))
         if (damp.ne.0.0d0 .and. pgamma.ne.0.0d0) then
            damp = pgamma * (r/damp)**3
            if (damp .lt. 50.0d0) then
               expdamp = exp(-damp)
               dmpik(3) = 1.0d0 - expdamp
               dmpik(5) = 1.0d0 - expdamp*(1.0d0+damp)
               if (rorder .ge. 7) then
                  damp2 = damp * damp
                  dmpik(7) = 1.0d0 - expdamp*(1.0d0+damp+0.6d0*damp2)
                  if (rorder .ge. 9) then
                     damp3 = damp * damp2
                     dmpik(9) = 1.0d0 - expdamp*(1.0d0+damp
     &                                     +(18.0d0/35.0d0)*damp2
     &                                     +(9.0d0/35.0d0)*damp3)
                  end if
               end if
            end if
         end if
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine dampthole2  --  original Thole damping values  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "dampthole2" finds coefficients for the original Thole damping
c     function used by AMOEBA and for mutual polarization by AMOEBA+
c
c     literature reference:
c
c     B. T. Thole, "Molecular Polarizabilities Calculated with a
c     Modified Dipole Interaction", Chemical Physics, 59, 341-350 (1981)
c
c
      subroutine dampthole2 (i,k,rorder,r,dmpik)
      use polar
      implicit none
      integer i,j,k
      integer rorder
      real*8 r,damp
      real*8 damp2
      real*8 damp3
      real*8 expdamp
      real*8 pgamma
      real*8 dmpik(*)
c
c
c     initialize the Thole damping factors to a value of one
c
      do j = 1, rorder
         dmpik(j) = 1.0d0
      end do
c
c     assign original Thole polarization model damping factors
c
      damp = pdamp(i) * pdamp(k)
      pgamma = min(thole(i),thole(k))
      if (pgamma .eq. 0.0d0)  pgamma = max(thole(i),thole(k))
      if (damp.ne.0.0d0 .and. pgamma.ne.0.0d0) then
         damp = pgamma * (r/damp)**3
         if (damp .lt. 50.0d0) then
            expdamp = exp(-damp)
            dmpik(3) = 1.0d0 - expdamp
            dmpik(5) = 1.0d0 - expdamp*(1.0d0+damp)
            if (rorder .ge. 7) then
               damp2 = damp * damp
               dmpik(7) = 1.0d0 - expdamp*(1.0d0+damp+0.6d0*damp2)
               if (rorder .ge. 9) then
                  damp3 = damp * damp2
                  dmpik(9) = 1.0d0 - expdamp*(1.0d0+damp
     &                                  +(18.0d0/35.0d0)*damp2
     &                                  +(9.0d0/35.0d0)*damp3)
               end if
            end if
         end if
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine damppole  --  penetration damping coefficents  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "damppole" finds coefficients for two alternative Gordon charge
c     penetration damping function
c
c     literature references:
c
c     L. V. Slipchenko and M. S. Gordon, "Electrostatic Energy in the
c     Effective Fragment Potential Method: Theory and Application to
c     the Benzene Dimer", Journal of Computational Chemistry, 28,
c     276-291 (2007)  [Gordon f1 and f2 models]
c
c     J. A. Rackers, Q. Wang, C. Liu, J.-P. Piquemal, P. Ren and
c     J. W. Ponder, "An Optimized Charge Penetration Model for Use with
c     the AMOEBA Force Field", Physical Chemistry Chemical Physics, 19,
c     276-291 (2017)
c
c
      subroutine damppole (r,rorder,alphai,alphak,dmpi,dmpk,dmpik)
      use mplpot
      implicit none
      integer rorder
      real*8 termi,termk
      real*8 termi2,termk2
      real*8 alphai,alphak
      real*8 alphai2,alphak2
      real*8 r,eps,diff
      real*8 expi,expk
      real*8 dampi,dampk
      real*8 dampi2,dampi3
      real*8 dampi4,dampi5
      real*8 dampi6,dampi7
      real*8 dampi8
      real*8 dampk2,dampk3
      real*8 dampk4,dampk5
      real*8 dampk6
      real*8 dmpi(*)
      real*8 dmpk(*)
      real*8 dmpik(*)
c
c
c     compute tolerance and exponential damping factors
c
      eps = 0.001d0
      diff = abs(alphai-alphak)
      dampi = alphai * r
      dampk = alphak * r
      expi = exp(-dampi)
      expk = exp(-dampk)
c
c     core-valence charge penetration damping for Gordon f1
c
      if (pentyp .eq. 'GORDON1') then
         dampi2 = dampi * dampi
         dampi3 = dampi * dampi2
         dampi4 = dampi2 * dampi2
         dampi5 = dampi2 * dampi3
         dmpi(1) = 1.0d0 - (1.0d0 + 0.5d0*dampi)*expi
         dmpi(3) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi2)*expi
         dmpi(5) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi2
     &                + dampi3/6.0d0)*expi
         dmpi(7) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi2
     &                + dampi3/6.0d0 + dampi4/30.0d0)*expi
         dmpi(9) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi2
     &                + dampi3/6.0d0 + 4.0d0*dampi4/105.0d0
     &                + dampi5/210.0d0)*expi
         if (diff .lt. eps) then
            dmpk(1) = dmpi(1)
            dmpk(3) = dmpi(3)
            dmpk(5) = dmpi(5)
            dmpk(7) = dmpi(7)
            dmpk(9) = dmpi(9)
         else
            dampk2 = dampk * dampk
            dampk3 = dampk * dampk2
            dampk4 = dampk2 * dampk2
            dampk5 = dampk2 * dampk3
            dmpk(1) = 1.0d0 - (1.0d0 + 0.5d0*dampk)*expk
            dmpk(3) = 1.0d0 - (1.0d0 + dampk + 0.5d0*dampk2)*expk
            dmpk(5) = 1.0d0 - (1.0d0 + dampk + 0.5d0*dampk2
     &                   + dampk3/6.0d0)*expk
            dmpk(7) = 1.0d0 - (1.0d0 + dampk + 0.5d0*dampk2
     &                   + dampk3/6.0d0 + dampk4/30.0d0)*expk
            dmpk(9) = 1.0d0 - (1.0d0 + dampk + 0.5d0*dampk2
     &                   + dampk3/6.0d0 + 4.0d0*dampk4/105.0d0
     &                   + dampk5/210.0d0)*expk
         end if
c
c     valence-valence charge penetration damping for Gordon f1
c
         if (diff .lt. eps) then
            dampi6 = dampi3 * dampi3
            dampi7 = dampi3 * dampi4
            dmpik(1) = 1.0d0 - (1.0d0 + 11.0d0*dampi/16.0d0
     &                    + 3.0d0*dampi2/16.0d0
     &                    + dampi3/48.0d0)*expi
            dmpik(3) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi2
     &                    + 7.0d0*dampi3/48.0d0
     &                    + dampi4/48.0d0)*expi
            dmpik(5) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi2
     &                    + dampi3/6.0d0 + dampi4/24.0d0
     &                    + dampi5/144.0d0)*expi
            dmpik(7) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi2
     &                    + dampi3/6.0d0 + dampi4/24.0d0
     &                    + dampi5/120.0d0 + dampi6/720.0d0)*expi
            dmpik(9) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi2
     &                    + dampi3/6.0d0 + dampi4/24.0d0
     &                    + dampi5/120.0d0 + dampi6/720.0d0
     &                    + dampi7/5040.0d0)*expi
            if (rorder .ge. 11) then
               dampi8 = dampi4 * dampi4
               dmpik(11) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi2
     &                        + dampi3/6.0d0 + dampi4/24.0d0
     &                        + dampi5/120.0d0 + dampi6/720.0d0
     &                        + dampi7/5040.0d0 + dampi8/45360.0d0)*expi
            end if
         else
            alphai2 = alphai * alphai
            alphak2 = alphak * alphak
            termi = alphak2 / (alphak2-alphai2)
            termk = alphai2 / (alphai2-alphak2)
            termi2 = termi * termi
            termk2 = termk * termk
            dmpik(1) = 1.0d0 - termi2*(1.0d0 + 2.0d0*termk
     &                    + 0.5d0*dampi)*expi
     &                 - termk2*(1.0d0 + 2.0d0*termi
     &                      + 0.5d0*dampk)*expk
            dmpik(3) = 1.0d0 - termi2*(1.0d0+dampi+0.5d0*dampi2)*expi
     &                    - termk2*(1.0d0+dampk+0.5d0*dampk2)*expk
     &                    - 2.0d0*termi2*termk*(1.0d0+dampi)*expi
     &                    - 2.0d0*termk2*termi*(1.0d0+dampk)*expk
            dmpik(5) = 1.0d0 - termi2*(1.0d0 + dampi + 0.5d0*dampi2
     &                    + dampi3/6.0d0)*expi
     &                 - termk2*(1.0d0 + dampk + 0.5d0*dampk2
     &                      + dampk3/6.0d0)*expk
     &                 - 2.0d0*termi2*termk
     &                      *(1.0d0 + dampi + dampi2/3.0d0)*expi
     &                 - 2.0d0*termk2*termi
     &                      *(1.0d0 + dampk + dampk2/3.0d0)*expk
            dmpik(7) = 1.0d0 - termi2*(1.0d0 + dampi + 0.5d0*dampi2
     &                    + dampi3/6.0d0 + dampi4/30.0d0)*expi
     &                 - termk2*(1.0d0 + dampk + 0.5d0*dampk2
     &                      + dampk3/6.0d0 + dampk4/30.0d0)*expk
     &                 - 2.0d0*termi2*termk*(1.0d0 + dampi
     &                      + 2.0d0*dampi2/5.0d0 + dampi3/15.0d0)*expi
     &                 - 2.0d0*termk2*termi*(1.0d0 + dampk
     &                      + 2.0d0*dampk2/5.0d0 + dampk3/15.0d0)*expk
            dmpik(9) = 1.0d0 - termi2*(1.0d0 + dampi + 0.5d0*dampi2
     &                    + dampi3/6.0d0 + 4.0d0*dampi4/105.0d0
     &                    + dampi5/210.0d0)*expi
     &                 - termk2*(1.0d0 + dampk + 0.5d0*dampk2
     &                      + dampk3/6.0d0 + 4.0d0*dampk4/105.0d0
     &                      + dampk5/210.0d0)*expk
     &                 - 2.0d0*termi2*termk*(1.0d0 + dampi
     &                      + 3.0d0*dampi2/7.0d0
     &                      + 2.0d0*dampi3/21.0d0
     &                      + dampi4/105.0d0)*expi 
     &                 - 2.0d0*termk2*termi*(1.0d0 + dampk
     &                      + 3.0d0*dampk2/7.0d0
     &                      + 2.0d0*dampk3/21.0d0
     &                      + dampk4/105.0d0)*expk
            if (rorder .ge. 11) then
               dampi6 = dampi3 * dampi3
               dampk6 = dampk3 * dampk3
               dmpik(11) = 1.0d0 - termi2*(1.0d0 + dampi
     &                        + 0.5d0*dampi2 + dampi3/6.0d0
     &                        + 5.0d0*dampi4/126.0d0
     &                        + 2.0d0*dampi5/315.0d0
     &                        + dampi6/1890.0d0)*expi
     &                     - termk2*(1.0d0 + dampk
     &                          + 0.5d0*dampk2 + dampk3/6.0d0
     &                          + 5.0d0*dampk4/126.0d0
     &                          + 2.0d0*dampk5/315.0d0
     &                          + dampk6/1890.0d0)*expk
     &                     - 2.0d0*termi2*termk*(1.0d0 + dampi
     &                          + 4.0d0*dampi2/9.0d0 + dampi3/9.0d0
     &                          + dampi4/63.0d0 + dampi5/945.0d0)*expi
     &                     - 2.0d0*termk2*termi*(1.0d0 + dampk
     &                          + 4.0d0*dampk2/9.0d0 + dampk3/9.0d0
     &                          + dampk4/63.0d0 + dampk5/945.0d0)*expk 
            end if
         end if
c
c     core-valence charge penetration damping for Gordon f2
c
      else if (pentyp .eq. 'GORDON2') then
         dampi2 = dampi * dampi
         dampi3 = dampi * dampi2
         dmpi(1) = 1.0d0 - expi
         dmpi(3) = 1.0d0 - (1.0d0 + dampi)*expi
         dmpi(5) = 1.0d0 - (1.0d0 + dampi + dampi2/3.0d0)*expi
         dmpi(7) = 1.0d0 - (1.0d0 + dampi + 0.4d0*dampi2
     &                + dampi3/15.0d0)*expi
         if (diff .lt. eps) then
            dmpk(1) = dmpi(1)
            dmpk(3) = dmpi(3)
            dmpk(5) = dmpi(5)
            dmpk(7) = dmpi(7)
         else
            dampk2 = dampk * dampk
            dampk3 = dampk * dampk2
            dmpk(1) = 1.0d0 - expk
            dmpk(3) = 1.0d0 - (1.0d0 + dampk)*expk
            dmpk(5) = 1.0d0 - (1.0d0 + dampk + dampk2/3.0d0)*expk
            dmpk(7) = 1.0d0 - (1.0d0 + dampk + 0.4d0*dampk2
     &                   + dampk3/15.0d0)*expk
         end if
c
c     valence-valence charge penetration damping for Gordon f2
c
         dampi4 = dampi2 * dampi2
         dampi5 = dampi2 * dampi3
         if (diff .lt. eps) then
            dampi6 = dampi3 * dampi3
            dmpik(1) = 1.0d0 - (1.0d0 + 0.5d0*dampi)*expi
            dmpik(3) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi2)*expi
            dmpik(5) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi2
     &                    + dampi3/6.0d0)*expi
            dmpik(7) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi2
     &                    + dampi3/6.0d0 + dampi4/30.0d0)*expi
            dmpik(9) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi2
     &                    + dampi3/6.0d0 + 4.0d0*dampi4/105.0d0
     &                    + dampi5/210.0d0)*expi
            if (rorder .ge. 11) then
               dmpik(11) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi2
     &                        + dampi3/6.0d0 + 5.0d0*dampi4/126.0d0
     &                        + 2.0d0*dampi5/315.0d0
     &                        + dampi6/1890.0d0)*expi
            end if
         else
            dampk4 = dampk2 * dampk2
            dampk5 = dampk2 * dampk3
            alphai2 = alphai * alphai
            alphak2 = alphak * alphak
            termi = alphak2 / (alphak2-alphai2)
            termk = alphai2 / (alphai2-alphak2)
            dmpik(1) = 1.0d0 - termi*expi - termk*expk
            dmpik(3) = 1.0d0 - termi*(1.0d0 + dampi)*expi
     &                    - termk*(1.0d0 + dampk)*expk
            dmpik(5) = 1.0d0 - termi*(1.0d0 + dampi + dampi2/3.0d0)*expi
     &                    - termk*(1.0d0 + dampk + dampk2/3.0d0)*expk
            dmpik(7) = 1.0d0 - termi*(1.0d0 + dampi + 0.4d0*dampi2
     &                    + dampi3/15.0d0)*expi
     &                 - termk*(1.0d0 + dampk + 0.4d0*dampk2
     &                      + dampk3/15.0d0)*expk
            dmpik(9) = 1.0d0 - termi*(1.0d0 + dampi + 3.0d0*dampi2/7.0d0
     &                    + 2.0d0*dampi3/21.0d0 + dampi4/105.0d0)*expi
     &                 - termk*(1.0d0 + dampk + 3.0d0*dampk2/7.0d0
     &                      + 2.0d0*dampk3/21.0d0 + dampk4/105.0d0)*expk
            if (rorder .ge. 11) then
               dmpik(11) = 1.0d0 - termi*(1.0d0 + dampi
     &                        + 4.0d0*dampi2/9.0d0 + dampi3/9.0d0
     &                        + dampi4/63.0d0 + dampi5/945.0d0)*expi
     &                     - termk*(1.0d0 + dampk
     &                          + 4.0d0*dampk2/9.0d0 + dampk3/9.0d0
     &                          + dampk4/63.0d0 + dampk5/945.0d0)*expk
            end if
         end if
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine dampdir  --  direct field damping coefficents  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "dampdir" finds coefficients for two alternative Gordon direct
c     field damping functions
c
c
      subroutine dampdir (r,alphai,alphak,dmpi,dmpk)
      use mplpot
      implicit none
      real*8 alphai,alphak
      real*8 r,eps,diff
      real*8 expi,expk
      real*8 dampi,dampk
      real*8 dampi2,dampk2
      real*8 dampi3,dampk3
      real*8 dampi4,dampk4
      real*8 dmpi(*)
      real*8 dmpk(*)
c
c
c     compute tolerance and exponential damping factors
c
      eps = 0.001d0
      diff = abs(alphai-alphak)
      dampi = alphai * r
      dampk = alphak * r
      expi = exp(-dampi)
      expk = exp(-dampk)
c
c     core-valence charge penetration damping for Gordon f1
c
      if (pentyp .eq. 'GORDON1') then
         dampi2 = dampi * dampi
         dampi3 = dampi * dampi2
         dampi4 = dampi2 * dampi2
         dmpi(3) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi2)*expi
         dmpi(5) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi2 
     &                + dampi3/6.0d0)*expi
         dmpi(7) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi2
     &                + dampi3/6.0d0 + dampi4/30.0d0)*expi
         if (diff .lt. eps) then
            dmpk(3) = dmpi(3)
            dmpk(5) = dmpi(5)
            dmpk(7) = dmpi(7)
         else
            dampk2 = dampk * dampk
            dampk3 = dampk * dampk2
            dampk4 = dampk2 * dampk2
            dmpk(3) = 1.0d0 - (1.0d0 + dampk + 0.5d0*dampk2)*expk
            dmpk(5) = 1.0d0 - (1.0d0 + dampk + 0.5d0*dampk2
     &                   + dampk3/6.0d0)*expk
            dmpk(7) = 1.0d0 - (1.0d0 + dampk + 0.5d0*dampk2
     &                   + dampk3/6.0d0 + dampk4/30.0d0)*expk
         end if
c
c     core-valence charge penetration damping for Gordon f2
c
      else if (pentyp .eq. 'GORDON2') then
         dampi2 = dampi * dampi
         dampi3 = dampi * dampi2
         dmpi(3) = 1.0d0 - (1.0d0 + dampi)*expi
         dmpi(5) = 1.0d0 - (1.0d0 + dampi + dampi2/3.0d0)*expi
         dmpi(7) = 1.0d0 - (1.0d0 + dampi + 0.4d0*dampi2
     &                + dampi3/15.0d0)*expi
         if (diff .lt. eps) then
            dmpk(3) = dmpi(3)
            dmpk(5) = dmpi(5)
            dmpk(7) = dmpi(7)
         else
            dampk2 = dampk * dampk
            dampk3 = dampk * dampk2
            dmpk(3) = 1.0d0 - (1.0d0 + dampk)*expk
            dmpk(5) = 1.0d0 - (1.0d0 + dampk + dampk2/3.0d0)*expk
            dmpk(7) = 1.0d0 - (1.0d0 + dampk + 0.4d0*dampk2
     &                   + dampk3/15.0d0)*expk
         end if
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine dampmut  --  mutual field damping coefficents  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "dampmut" finds coefficients for two alternative Gordon mutual
c     field damping functions
c
c
      subroutine dampmut (r,alphai,alphak,dmpik)
      use mplpot
      implicit none
      real*8 termi,termk
      real*8 termi2,termk2
      real*8 alphai,alphak
      real*8 alphai2,alphak2
      real*8 r,eps,diff
      real*8 expi,expk
      real*8 dampi,dampk
      real*8 dampi2,dampi3
      real*8 dampi4,dampi5
      real*8 dampk2,dampk3
      real*8 dmpik(*)
c
c
c     compute tolerance and exponential damping factors
c
      eps = 0.001d0
      diff = abs(alphai-alphak)
      dampi = alphai * r
      dampk = alphak * r
      expi = exp(-dampi)
      expk = exp(-dampk)
c
c     valence-valence charge penetration damping for Gordon f1
c
      if (pentyp .eq. 'GORDON1') then
         dampi2 = dampi * dampi
         dampi3 = dampi * dampi2
         if (diff .lt. eps) then
            dampi4 = dampi2 * dampi2
            dampi5 = dampi2 * dampi3
            dmpik(3) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi2
     &                    + 7.0d0*dampi3/48.0d0
     &                    + dampi4/48.0d0)*expi
            dmpik(5) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi2
     &                    + dampi3/6.0d0 + dampi4/24.0d0
     &                    + dampi5/144.0d0)*expi
         else
            dampk2 = dampk * dampk
            dampk3 = dampk * dampk2
            alphai2 = alphai * alphai
            alphak2 = alphak * alphak
            termi = alphak2 / (alphak2-alphai2)
            termk = alphai2 / (alphai2-alphak2)
            termi2 = termi * termi
            termk2 = termk * termk
            dmpik(3) = 1.0d0 - termi2*(1.0d0+dampi+0.5d0*dampi2)*expi
     &                    - termk2*(1.0d0+dampk+0.5d0*dampk2)*expk
     &                    - 2.0d0*termi2*termk*(1.0d0+dampi)*expi
     &                    - 2.0d0*termk2*termi*(1.0d0+dampk)*expk
            dmpik(5) = 1.0d0 - termi2*(1.0d0+dampi+0.5d0*dampi2
     &                            +dampi3/6.0d0)*expi
     &                    - termk2*(1.0d0+dampk+0.5d0*dampk2
     &                         +dampk3/6.00)*expk
     &                    - 2.0d0*termi2*termk
     &                         *(1.0+dampi+dampi2/3.0d0)*expi
     &                    - 2.0d0*termk2*termi
     &                         *(1.0+dampk+dampk2/3.0d0)*expk
         end if
c
c     valence-valence charge penetration damping for Gordon f2
c
      else if (pentyp .eq. 'GORDON2') then
         dampi2 = dampi * dampi
         if (diff .lt. eps) then
            dampi3 = dampi * dampi2
            dmpik(3) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi2)*expi
            dmpik(5) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi2
     &                    + dampi3/6.0d0)*expi
         else
            dampk2 = dampk * dampk
            alphai2 = alphai * alphai
            alphak2 = alphak * alphak
            termi = alphak2 / (alphak2-alphai2)
            termk = alphai2 / (alphai2-alphak2)
            dmpik(3) = 1.0d0 - termi*(1.0d0 + dampi)*expi
     &                    - termk*(1.0d0 + dampk)*expk
            dmpik(5) = 1.0d0 - termi*(1.0d0 + dampi + dampi2/3.0d0)*expi
     &                    - termk*(1.0d0 + dampk + dampk2/3.0d0)*expk
         end if
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine damppot  --  electrostatic potential damping  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "damppot" finds coefficients for two alternative Gordon charge
c     penetration damping functions for the electrostatic potential
c
c
      subroutine damppot (r,alphak,dmpk)
      use mplpot
      implicit none
      real*8 r,alphak
      real*8 expk,dampk
      real*8 dampk2,dampk3
      real*8 dmpk(*)
c
c
c     compute common exponential factors for damping
c
      dampk = alphak * r
      expk = exp(-dampk)
c
c     core-valence charge penetration damping for Gordon f1
c
      if (pentyp .eq. 'GORDON1') then
         dampk2 = dampk * dampk
         dampk3 = dampk * dampk2
         dmpk(1) = 1.0d0 - (1.0d0 + 0.5d0*dampk)*expk
         dmpk(3) = 1.0d0 - (1.0d0 + dampk + 0.5d0*dampk2)*expk
         dmpk(5) = 1.0d0 - (1.0d0 + dampk + 0.5d0*dampk2
     &                + dampk3/6.0d0)*expk
c
c     core-valence charge penetration damping for Gordon f2
c
      else if (pentyp .eq. 'GORDON2') then
         dampk2 = dampk * dampk
         dmpk(1) = 1.0d0 - expk
         dmpk(3) = 1.0d0 - (1.0d0 + dampk)*expk
         dmpk(5) = 1.0d0 - (1.0d0 + dampk + dampk2/3.0d0)*expk
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine damprep  --  Pauli exchange repulsion damping  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "damprep" finds coefficients for the Pauli repulsion damping
c     function used by HIPPO
c
c     literature reference:
c
c     J. A. Rackers and J. W. Ponder, "Classical Pauli Repulsion: An
c     Anisotropic, Atomic Multipole Model", Journal of Chemical Physics,
c     150, 084104 (2019)
c
c
      subroutine damprep (r,r2,rr1,rr3,rr5,rr7,rr9,rr11,
     &                       rorder,dmpi,dmpk,dmpik)
      implicit none
      integer rorder
      real*8 r,r2,r3,r4
      real*8 r5,r6,r7,r8
      real*8 rr1,rr3,rr5
      real*8 rr7,rr9,rr11
      real*8 s,ds,d2s
      real*8 d3s,d4s,d5s
      real*8 dmpi,dmpk
      real*8 dmpi2,dmpk2
      real*8 dmpi22,dmpi23
      real*8 dmpi24,dmpi25
      real*8 dmpi26,dmpi27
      real*8 dmpk22,dmpk23
      real*8 dmpk24,dmpk25
      real*8 dmpk26
      real*8 eps,diff
      real*8 expi,expk
      real*8 dampi,dampk
      real*8 pre,term,tmp
      real*8 dmpik(*)
c
c
c     compute tolerance value for damping exponents
c
      eps = 0.001d0
      diff = abs(dmpi-dmpk)
c
c     treat the case where alpha damping exponents are equal
c
      if (diff .lt. eps) then
         r3 = r2 * r
         r4 = r3 * r
         r5 = r4 * r
         r6 = r5 * r
         dmpi2 = 0.5d0 * dmpi
         dampi = dmpi2 * r
         expi = exp(-dampi)
         dmpi22 = dmpi2 * dmpi2
         dmpi23 = dmpi22 * dmpi2
         dmpi24 = dmpi23 * dmpi2
         dmpi25 = dmpi24 * dmpi2
         pre = 2.0d0
         s = (r + dmpi2*r2 + dmpi22*r3/3.0d0) * expi
         ds = (dmpi22*r3 + dmpi23*r4) * expi / 3.0d0
         d2s = dmpi24 * expi * r5 / 9.0d0
         d3s = dmpi25 * expi * r6 / 45.0d0
         if (rorder .ge. 9) then
            r7 = r6 * r
            dmpi26 = dmpi25 * dmpi2
            d4s = (dmpi25*r6 + dmpi26*r7) * expi / 315.0d0
            if (rorder .ge. 11) then
               r8 = r7 * r
               dmpi27 = dmpi2 * dmpi26
               d5s = (dmpi25*r6 + dmpi26*r7 + dmpi27*r8/3.0d0)
     &                   * expi / 945.0d0
            end if
         end if
c
c     treat the case where alpha damping exponents are unequal
c
      else
         r3 = r2 * r
         r4 = r3 * r
         dmpi2 = 0.5d0 * dmpi
         dmpk2 = 0.5d0 * dmpk
         dampi = dmpi2 * r
         dampk = dmpk2 * r
         expi = exp(-dampi)
         expk = exp(-dampk)
         dmpi22 = dmpi2 * dmpi2
         dmpi23 = dmpi22 * dmpi2
         dmpi24 = dmpi23 * dmpi2
         dmpk22 = dmpk2 * dmpk2
         dmpk23 = dmpk22 * dmpk2
         dmpk24 = dmpk23 * dmpk2
         term = dmpi22 - dmpk22
         pre = 128.0d0 * dmpi23 * dmpk23 / term**4
         tmp = 4.0d0 * dmpi2 * dmpk2 / term
         s = (dampi-tmp)*expk + (dampk+tmp)*expi
         ds = (dmpi2*dmpk2*r2 - 4.0d0*dmpi2*dmpk22*r/term
     &            - 4.0d0*dmpi2*dmpk2/term) * expk
     &      + (dmpi2*dmpk2*r2 + 4.0d0*dmpi22*dmpk2*r/term
     &            + 4.0d0*dmpi2*dmpk2/term) * expi
         d2s = (dmpi2*dmpk2*r2/3.0d0
     &             + dmpi2*dmpk22*r3/3.0d0
     &             - (4.0d0/3.0d0)*dmpi2*dmpk23*r2/term
     &             - 4.0d0*dmpi2*dmpk22*r/term
     &             - 4.0d0*dmpi2*dmpk2/term) * expk
     &       + (dmpi2*dmpk2*r2/3.0d0
     &             + dmpi22*dmpk2*r3/3.0d0
     &             + (4.0d0/3.0d0)*dmpi23*dmpk2*r2/term
     &             + 4.0d0*dmpi22*dmpk2*r/term
     &             + 4.0d0*dmpi2*dmpk2/term) * expi
         d3s = (dmpi2*dmpk23*r4/15.0d0
     &             + dmpi2*dmpk22*r3/5.0d0
     &             + dmpi2*dmpk2*r2/5.0d0
     &             - (4.0d0/15.0d0)*dmpi2*dmpk24*r3/term
     &             - (8.0d0/5.0d0)*dmpi2*dmpk23*r2/term
     &             - 4.0d0*dmpi2*dmpk22*r/term
     &             - 4.0d0/term*dmpi2*dmpk2) * expk
     &       + (dmpi23*dmpk2*r4/15.0d0 
     &             + dmpi22*dmpk2*r3/5.0d0
     &             + dmpi2*dmpk2*r2/5.0d0 
     &             + (4.0d0/15.0d0)*dmpi24*dmpk2*r3/term
     &             + (8.0d0/5.0d0)*dmpi23*dmpk2*r2/term
     &             + 4.0d0*dmpi22*dmpk2*r/term
     &             + 4.0d0/term*dmpi2*dmpk2) * expi
         if (rorder .ge. 9) then
            r5 = r4 * r
            dmpi25 = dmpi24 * dmpi2
            dmpk25 = dmpk24 * dmpk2
            d4s = (dmpi2*dmpk24*r5/105.0d0
     &                + (2.0d0/35.0d0)*dmpi2*dmpk23*r4
     &                + dmpi2*dmpk22*r3/7.0d0
     &                + dmpi2*dmpk2*r2/7.0d0
     &                - (4.0d0/105.0d0)*dmpi2*dmpk25*r4/term
     &                - (8.0d0/21.0d0)*dmpi2*dmpk24*r3/term
     &                - (12.0d0/7.0d0)*dmpi2*dmpk23*r2/term
     &                - 4.0d0*dmpi2*dmpk22*r/term
     &                - 4.0d0*dmpi2*dmpk2/term) * expk
     &          + (dmpi24*dmpk2*r5/105.0d0
     &                + (2.0d0/35.0d0)*dmpi23*dmpk2*r4
     &                + dmpi22*dmpk2*r3/7.0d0
     &                + dmpi2*dmpk2*r2/7.0d0
     &                + (4.0d0/105.0d0)*dmpi25*dmpk2*r4/term
     &                + (8.0d0/21.0d0)*dmpi24*dmpk2*r3/term
     &                + (12.0d0/7.0d0)*dmpi23*dmpk2*r2/term
     &                + 4.0d0*dmpi22*dmpk2*r/term
     &                + 4.0d0*dmpi2*dmpk2/term) * expi
            if (rorder .ge. 11) then
               r6 = r5 * r
               dmpi26 = dmpi25 * dmpi2
               dmpk26 = dmpk25 * dmpk2
               d5s = (dmpi2*dmpk25*r6/945.0d0
     &                   + (2.0d0/189.0d0)*dmpi2*dmpk24*r5
     &                   + dmpi2*dmpk23*r4/21.0d0
     &                   + dmpi2*dmpk22*r3/9.0d0
     &                   + dmpi2*dmpk2*r2/9.0d0
     &                   - (4.0d0/945.0d0)*dmpi2*dmpk26*r5/term
     &                   - (4.0d0/63.0d0)*dmpi2*dmpk25*r4/term
     &                   - (4.0d0/9.0d0)*dmpi2*dmpk24*r3/term
     &                   - (16.0d0/9.0d0)*dmpi2*dmpk23*r2/term
     &                   - 4.0d0*dmpi2*dmpk22*r/term
     &                   - 4.0d0*dmpi2*dmpk2/term) * expk
     &             + (dmpi25*dmpk2*r6/945.0d0
     &                   + (2.0d0/189.0d0)*dmpi24*dmpk2*r5
     &                   + dmpi23*dmpk2*r4/21.0d0
     &                   + dmpi22*dmpk2*r3/9.0d0
     &                   + dmpi2*dmpk2*r2/9.0d0
     &                   + (4.0d0/945.0d0)*dmpi26*dmpk2*r5/term
     &                   + (4.0d0/63.0d0)*dmpi25*dmpk2*r4/term
     &                   + (4.0d0/9.0d0)*dmpi24*dmpk2*r3/term
     &                   + (16.0d0/9.0d0)*dmpi23*dmpk2*r2/term
     &                   + 4.0d0*dmpi22*dmpk2*r/term
     &                   + 4.0d0*dmpi2*dmpk2/term) * expi
            end if
         end if
      end if
c
c     convert partial derivatives into full derivatives
c
      s = s * rr1
      ds = ds * rr3
      d2s = d2s * rr5
      d3s = d3s * rr7
      dmpik(1) = 0.5d0 * pre * s * s
      dmpik(3) = pre * s * ds
      dmpik(5) = pre * (s*d2s + ds*ds)
      dmpik(7) = pre * (s*d3s + 3.0d0*ds*d2s)
      if (rorder .ge. 9) then
         d4s = d4s * rr9
         dmpik(9) = pre * (s*d4s + 4.0d0*ds*d3s + 3.0d0*d2s*d2s)
         if (rorder .ge. 11) then
            d5s = d5s * rr11
            dmpik(11) = pre * (s*d5s + 5.0d0*ds*d4s + 10.0d0*d2s*d3s)
         end if
      end if
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine dampep  --  exchange polarization damping  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "dampep" finds overlap for exchange polarization damping
c     function used by HIPPO
c
c     literature reference:
c
c
      subroutine dampep (r,preik,alphai,alphak,s2,ds2)
      use expol
      implicit none
      real*8 r
      real*8 s,s2,ds2
      real*8 alphai,alphak,alphaik
      real*8 dmpi2,dmpk2
      real*8 dmpi22,dmpk22
      real*8 dmpik2,dampik,dampik2
      real*8 eps,diff
      real*8 expi,expk,expik
      real*8 dampi,dampk,dampi2
      real*8 pre,term,preik
c
c
      if (scrtyp .eq. 'S2U') then
         alphaik = sqrt(alphai * alphak)
         dmpik2 = 0.5d0 * alphaik
         dampik = dmpik2 * r
         dampik2 = dampik * dampik
         expik = exp(-dampik)
         s =(1+dampik+dampik2/3.0d0)*expik
         s2 = s*s
         ds2 = s * (-alphaik/3.0d0)*(dampik+dampik2)*expik
      else if (scrtyp .eq. 'S2 ') then
c
c     compute tolerance value for damping exponents
c
         eps = 0.001d0
         diff = abs(alphai-alphak)
c
c     treat the case where alpha damping exponents are equal
c
         if (diff .lt. eps) then
            dmpi2 = 0.5d0 * alphai
            dampi = dmpi2 * r
            dampi2 = dampi * dampi
            expi = exp(-dampi)
            s = (1+dampi+dampi2/3.0d0)*expi
            ds2 = s * (-alphai/3.0d0)*(dampi+dampi2)*expi
c
c     treat the case where alpha damping exponents are unequal
c
         else
            dmpi2 = 0.5d0 * alphai
            dmpk2 = 0.5d0 * alphak
            dampi = dmpi2 * r
            dampk = dmpk2 * r
            expi = exp(-dampi)
            expk = exp(-dampk)
            dmpi22 = dmpi2 * dmpi2
            dmpk22 = dmpk2 * dmpk2
            term = dmpi22 - dmpk22
            pre = sqrt(alphai**3 * alphak**3) / (r * term**3)
            s = pre*(dmpi2*(r*term - 4*dmpk2) * expk
     &            + dmpk2*(r*term + 4*dmpi2) * expi)
            ds2 = 2.0d0*s*pre*dmpi2*dmpk2 *
     &       ((4.0d0/r-(r*term-4.0d0*dmpk2))*expk -
     &       ((4.0d0/r+(r*term+4.0d0*dmpi2))*expi))
         end if
         s2 = s*s
      else if (scrtyp .eq. 'G  ') then
         alphaik = sqrt(alphai * alphak)
         s2 = exp(-alphaik/10.0d0 * r**2)
         ds2 = (-alphaik/5.0d0)*r*s2
      end if
      s2 = preik*s2
      ds2 = preik*ds2
      return
      end
