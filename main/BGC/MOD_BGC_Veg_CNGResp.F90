#include <define.h>
#ifdef BGC
MODULE MOD_BGC_Veg_CNGResp

!-----------------------------------------------------------------------------------------------------
! !DESCRIPTION:
! This module calculate growth respiration rate.
!
! !REFERENCES:
! Atkin, O.K., Bahar, N.H., Bloomfield, K.J., Griffin, K.L., Heskel, M.A., Huntingford, C., de la Torre, A.M.
! and Turnbull, M.H., 2017. Leaf respiration in terrestrial biosphere models. Plant respiration: metabolic
! fluxes and carbon balance, pp.107-142.
!
! !ORIGINAL:
! The Community Land Model version 5.0 (CLM5)
!
! !REVISION:
! Xingjie Lu, 2021, revised the CLM5 code to be compatible with CoLM code structure.

   USE MOD_Precision
   USE MOD_Namelist, only: DEF_USE_SOILACIDIFICATION, DEF_SOILACID_FACTOR_ER, &
       DEF_SOILACID_pH_opt, DEF_SOILACID_pH_sens
   USE MOD_Const_PFT, only: &
       grperc, grpnow, woody, rootfr_p

   USE MOD_Vars_PFTimeInvariants, only: pftclass
   USE MOD_BGC_Vars_TimeVariables, only: pH_vr

   USE MOD_BGC_Vars_1DPFTFluxes, only: &
       cpool_to_leafc_p     , cpool_to_leafc_storage_p     , leafc_xfer_to_leafc_p          , &
       cpool_to_frootc_p    , cpool_to_frootc_storage_p    , frootc_xfer_to_frootc_p        , &
       cpool_to_livestemc_p , cpool_to_livestemc_storage_p , livestemc_xfer_to_livestemc_p  , &
       cpool_to_deadstemc_p , cpool_to_deadstemc_storage_p , deadstemc_xfer_to_deadstemc_p  , &
       cpool_to_livecrootc_p, cpool_to_livecrootc_storage_p, livecrootc_xfer_to_livecrootc_p, &
       cpool_to_deadcrootc_p, cpool_to_deadcrootc_storage_p, deadcrootc_xfer_to_deadcrootc_p, &
       cpool_to_grainc_p    , cpool_to_grainc_storage_p    , grainc_xfer_to_grainc_p        , &
       cpool_leaf_gr_p      , cpool_leaf_storage_gr_p      , transfer_leaf_gr_p             , &
       cpool_froot_gr_p     , cpool_froot_storage_gr_p     , transfer_froot_gr_p            , &
       cpool_livestem_gr_p  , cpool_livestem_storage_gr_p  , transfer_livestem_gr_p         , &
       cpool_deadstem_gr_p  , cpool_deadstem_storage_gr_p  , transfer_deadstem_gr_p         , &
       cpool_livecroot_gr_p , cpool_livecroot_storage_gr_p , transfer_livecroot_gr_p        , &
       cpool_deadcroot_gr_p , cpool_deadcroot_storage_gr_p , transfer_deadcroot_gr_p        , &
       cpool_grain_gr_p     , cpool_grain_storage_gr_p     , transfer_grain_gr_p

   IMPLICIT NONE

   PUBLIC CNGResp

CONTAINS

   SUBROUTINE CNGResp(i, ps, pe, nl_soil, npcropmin)

   integer ,intent(in) :: i         ! patch index
   integer ,intent(in) :: ps        ! start pft index
   integer ,intent(in) :: pe        ! end pft index
   integer ,intent(in) :: nl_soil   ! number of total soil layers
   integer ,intent(in) :: npcropmin ! first crop pft index

   ! !LOCAL VARIABLES:
   real(r8):: respfact_leaf
   real(r8):: respfact_froot
   real(r8):: respfact_livecroot
   real(r8):: respfact_livestem
   real(r8):: respfact_leaf_storage
   real(r8):: respfact_froot_storage
   real(r8):: respfact_livecroot_storage
   real(r8):: respfact_livestem_storage
   integer :: ivt, m, j
   real(r8) :: soil_pH_root_p, root_weight_sum, pH_factor_p, pH_opt, pH_sens


      DO m = ps, pe
         ivt = pftclass(m)
         IF(DEF_USE_SOILACIDIFICATION .and. DEF_SOILACID_FACTOR_ER) THEN
            pH_opt = DEF_SOILACID_pH_opt    ! optimal pH for soil microbial respiration
            pH_sens = DEF_SOILACID_pH_sens  ! sensitivity parameter (larger = less sensitive to pH changes)
            root_weight_sum = 0.0_r8
            soil_pH_root_p = 0.0_r8
            
            DO j = 1, nl_soil
               IF (rootfr_p(j,ivt) > 0.0_r8) THEN
                  soil_pH_root_p = soil_pH_root_p + pH_vr(j,i) * rootfr_p(j,ivt)
                  root_weight_sum = root_weight_sum + rootfr_p(j,ivt)
               ENDIF
            ENDDO

            IF (root_weight_sum > 0.0_r8) THEN
               soil_pH_root_p = soil_pH_root_p / root_weight_sum
            ELSE
               soil_pH_root_p = 6.5_r8
            ENDIF
            
            ! Gaussian function: h(pH) = exp(-((pH - pH_opt) / pH_sens)^2)
            pH_factor_p = EXP(-((soil_pH_root_p - pH_opt) / pH_sens)**2)
            write(*,*) 'Soil RR_G : pH =', soil_pH_root_p, ' pH_factor_p =', pH_factor_p
         ELSE
            pH_factor_p = 1.0_r8
         ENDIF
         respfact_leaf              = 1.0_r8
         respfact_froot             = 1.0_r8
         respfact_livecroot         = 1.0_r8
         respfact_livestem          = 1.0_r8
         respfact_livecroot         = 1.0_r8
         respfact_livestem          = 1.0_r8
         respfact_leaf_storage      = 1.0_r8
         respfact_froot_storage     = 1.0_r8
         respfact_livecroot_storage = 1.0_r8
         respfact_livestem_storage  = 1.0_r8
         respfact_livecroot_storage = 1.0_r8
         respfact_livestem_storage  = 1.0_r8

         IF (ivt >= npcropmin) THEN ! skip 2 generic crops
            cpool_livestem_gr_p         (m) = cpool_to_livestemc_p           (m) * grperc(ivt) * respfact_livestem

            cpool_livestem_storage_gr_p (m) = cpool_to_livestemc_storage_p   (m) * grperc(ivt) * grpnow(ivt) * respfact_livestem_storage

            transfer_livestem_gr_p      (m) = livestemc_xfer_to_livestemc_p  (m) * grperc(ivt) * (1._r8 - grpnow(ivt)) * respfact_livestem_storage

            cpool_grain_gr_p            (m) = cpool_to_grainc_p              (m) * grperc(ivt)

            cpool_grain_storage_gr_p    (m) = cpool_to_grainc_storage_p      (m) * grperc(ivt) * grpnow(ivt)

            transfer_grain_gr_p         (m) = grainc_xfer_to_grainc_p        (m) * grperc(ivt) * (1._r8 - grpnow(ivt))
         ENDIF

        ! leaf and fine root growth respiration
         cpool_leaf_gr_p                (m) = cpool_to_leafc_p               (m) * grperc(ivt) * respfact_leaf

         cpool_leaf_storage_gr_p        (m) = cpool_to_leafc_storage_p       (m) * grperc(ivt) * grpnow(ivt) * respfact_leaf_storage

         transfer_leaf_gr_p             (m) = leafc_xfer_to_leafc_p          (m) * grperc(ivt) * (1._r8 - grpnow(ivt)) * respfact_leaf_storage

         cpool_froot_gr_p               (m) = cpool_to_frootc_p              (m) * grperc(ivt) * respfact_froot * pH_factor_p

         cpool_froot_storage_gr_p       (m) = cpool_to_frootc_storage_p      (m) * grperc(ivt) * grpnow(ivt) * respfact_froot_storage * pH_factor_p

         transfer_froot_gr_p            (m) = frootc_xfer_to_frootc_p        (m) * grperc(ivt) * (1._r8 - grpnow(ivt)) * respfact_froot_storage * pH_factor_p

         IF (woody(ivt) == 1._r8) THEN
            cpool_livestem_gr_p         (m) = cpool_to_livestemc_p           (m) * grperc(ivt) * respfact_livestem

            cpool_livestem_storage_gr_p (m) = cpool_to_livestemc_storage_p   (m) * grperc(ivt) * grpnow(ivt) * respfact_livestem_storage

            transfer_livestem_gr_p      (m) = livestemc_xfer_to_livestemc_p  (m) * grperc(ivt) * (1._r8 - grpnow(ivt)) * respfact_livestem_storage

            cpool_deadstem_gr_p         (m) = cpool_to_deadstemc_p           (m) * grperc(ivt)

            cpool_deadstem_storage_gr_p (m) = cpool_to_deadstemc_storage_p   (m) * grperc(ivt) * grpnow(ivt)

            transfer_deadstem_gr_p      (m) = deadstemc_xfer_to_deadstemc_p  (m) * grperc(ivt) * (1._r8 - grpnow(ivt))

            cpool_livecroot_gr_p        (m) = cpool_to_livecrootc_p          (m) * grperc(ivt) * respfact_livecroot * pH_factor_p

            cpool_livecroot_storage_gr_p(m) = cpool_to_livecrootc_storage_p  (m) * grperc(ivt) * grpnow(ivt) * respfact_livecroot_storage * pH_factor_p

            transfer_livecroot_gr_p     (m) = livecrootc_xfer_to_livecrootc_p(m) * grperc(ivt) * (1._r8 - grpnow(ivt)) * respfact_livecroot_storage * pH_factor_p

            cpool_deadcroot_gr_p        (m) = cpool_to_deadcrootc_p          (m) * grperc(ivt) * pH_factor_p

            cpool_deadcroot_storage_gr_p(m) = cpool_to_deadcrootc_storage_p  (m) * grperc(ivt) * grpnow(ivt) * pH_factor_p

            transfer_deadcroot_gr_p     (m) = deadcrootc_xfer_to_deadcrootc_p(m) * grperc(ivt) * (1._r8 - grpnow(ivt)) * pH_factor_p
         ENDIF
      ENDDO

   END SUBROUTINE CNGResp

END MODULE MOD_BGC_Veg_CNGResp
#endif
