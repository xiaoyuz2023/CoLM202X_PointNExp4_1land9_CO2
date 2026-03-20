#include <define.h>
#ifdef BGC

MODULE MOD_BGC_Soil_BiogeochemNStateUpdate1

!---------------------------------------------------------------------------------------------
! !DESCRIPTION:
! Updates soil mineral nitrogen pool sizes. The dynamics of soil mineral nitrogen pool is
! simulated according to fertilisation, nitrogen deposition, biological fixation, plant uptake,
! mineralisation and immobilisation in this module. IF nitrification is activated, nitrate nitrogen
! has a separated pool against ammonium nitrogen pool. Accumulated nitrogen transfer
! network is also recorded for semi-analytic spinup.
!
! !ORIGINAL:
! The Community Land Model version 5.0 (CLM5.0)
!
! !REVISION:
! Xingjie Lu, 2022, 1) modify original CLM5 to be compatible with CoLM code structure.
!                   2) Record accumulated nitrogen transfer network for semi-analytic spinup

   USE MOD_Precision
   USE MOD_Namelist, only: DEF_USE_SASU, DEF_USE_DiagMatrix, DEF_USE_NITRIF, DEF_USE_CNSOYFIXN, &
       DEF_SOILACID_neu_const, DEF_SOILACID_acid_const, DEF_SOILACID_alk_const, &
       DEF_SOILACID_k_ref_co2, DEF_SOILACID_z_ref_co2, DEF_SOILACID_exchange_calib
   USE MOD_Const_Physical, only : denice, denh2o, tfrz
   USE MOD_Vars_TimeVariables, only : wliq_soisno, wice_soisno, t_soisno
   USE MOD_Vars_TimeInvariants, only: porsl, BD_all, wfc, bsw, OM_density, vf_clay, vf_sand, wf_clay, wf_sand
   USE MOD_Vars_1DForcing, only : forc_pco2m
   USE MOD_Const_PFT, only: rootfr_p
   USE MOD_Vars_PFTimeInvariants, only: pftclass, pftfrac
   USE MOD_Vars_Global, only: spval

   USE MOD_BGC_Vars_TimeInvariants, only: &
     ! bgc constants
       i_met_lit, i_cel_lit, i_lig_lit, i_cwd, i_soil1, i_soil2, i_soil3
   USE MOD_BGC_Vars_TimeInvariants, only: &
       receiver_pool, donor_pool, nitrif_n2o_loss_frac, organic_max, &
       d_con_g21_CO2, d_con_g22_CO2, d_con_w21_CO2, d_con_w22_CO2, d_con_w23_CO2, sf, sf_no3

   USE MOD_BGC_Vars_TimeVariables, only: &
       ! Mineral nitrogen pools (inout)
       sminn_vr                   , smin_nh4_vr                , smin_no3_vr                , &
       ndep_prof                  , nfixation_prof             , depth_scalar               , &
       pCO2_hr                    , pH_vr                      , pH_init_vr                 , h_conc_vr                  , &
       ca_conc_vr                 , Bc_conc_vr                 , Bct_conc_vr                , ANC_conc_vr                , &
       al3_conc_vr                , al2_conc_vr                , al_conc_vr                , oh_conc_vr                 , &
       DIC_conc_vr                , h2co3_conc_vr              , hco3_conc_vr               , co3_conc_vr                , & 
       h_conc_m1                  , oh_conc_m1                 , ca_conc_m1                 , Bc_conc_m1                 , &
       DIC_conc_m1                , &
       nh4_conc_vr                , no3_conc_vr                , h_from_N                   , E_Bc                       , &
       rate_cal                   , rate_CO2                   , rate_Bcex                  , rate_Hex                   , &
       pCO2_soil                  , calcite_mass               , calcite_content            , CEC                        , &
       AKX_met_to_soil1_n_vr_acc  , AKX_cel_to_soil1_n_vr_acc  , AKX_lig_to_soil2_n_vr_acc  , AKX_soil1_to_soil2_n_vr_acc, &
       AKX_cwd_to_cel_n_vr_acc    , AKX_cwd_to_lig_n_vr_acc    , AKX_soil1_to_soil3_n_vr_acc, AKX_soil2_to_soil1_n_vr_acc, &
       AKX_soil2_to_soil3_n_vr_acc, AKX_soil3_to_soil1_n_vr_acc, &
       AKX_met_exit_n_vr_acc      , AKX_cel_exit_n_vr_acc      , AKX_lig_exit_n_vr_acc      , AKX_cwd_exit_n_vr_acc      , &
       AKX_soil1_exit_n_vr_acc    , AKX_soil2_exit_n_vr_acc    , AKX_soil3_exit_n_vr_acc

   USE MOD_BGC_Vars_1DFluxes, only: &
       ! Decomposition fluxes variables (inout)
              decomp_npools_sourcesink, decomp_ntransfer_vr      , decomp_sminn_flux_vr     , sminn_to_denit_decomp_vr, &
              gross_nmin_vr           , actual_immob_nh4_vr      , actual_immob_no3_vr      , &
              sminn_to_plant_vr       , smin_nh4_to_plant_vr     , smin_no3_to_plant_vr     , supplement_to_sminn_vr, &
              sminn_to_plant_fun_vr   , sminn_to_plant_fun_nh4_vr, sminn_to_plant_fun_no3_vr, &
              sminn_to_denit_excess_vr, f_nit_vr                 , f_denit_vr               , soyfixn_to_sminn, &
              ndep_to_sminn           , ffix_to_sminn            , nfix_to_sminn            , fert_to_sminn , &
              soil_er                 , soil_er_vr               ,decomp_hr                  , ar            ,&
              h_content_cons_vr       , h_content_prod_vr       
   USE MOD_SPMD_Task

   IMPLICIT NONE

   PUBLIC SoilBiogeochemNStateUpdate1
   PUBLIC soil_acidification_buffering

CONTAINS

   SUBROUTINE SoilBiogeochemNStateUpdate1(i,deltim,nl_soil,ndecomp_transitions,dz_soi)

   integer ,intent(in) :: i                  ! patch idnex
   real(r8),intent(in) :: deltim             ! time step in seconds
   integer ,intent(in) :: nl_soil            ! number of total soil layers
   integer ,intent(in) :: ndecomp_transitions! number of total transfers between soil and litter pools in the decomposition
   real(r8),intent(in) :: dz_soi(1:nl_soil)  ! thicknesses of each soil layer (m)

   integer j,k
   real(r8):: sminflux,minerflux

      IF(.not. DEF_USE_NITRIF)THEN
         DO j = 1, nl_soil
            ! N deposition and fixation
            sminn_vr(j,i) = sminn_vr(j,i) + ndep_to_sminn(i)*deltim * ndep_prof(j,i)
            sminn_vr(j,i) = sminn_vr(j,i) + nfix_to_sminn(i)*deltim * nfixation_prof(j,i)
         ENDDO
      ELSE
         DO j = 1, nl_soil
            ! N deposition and fixation (put all into NH4 pool)
            smin_nh4_vr(j,i) = smin_nh4_vr(j,i) + ndep_to_sminn(i)*deltim * ndep_prof(j,i)
            smin_nh4_vr(j,i) = smin_nh4_vr(j,i) + nfix_to_sminn(i)*deltim * nfixation_prof(j,i)
         ENDDO
      ENDIF

      ! repeating N dep and fixation for crops
#ifdef CROP
      IF(.not. DEF_USE_NITRIF)THEN
         DO j = 1, nl_soil
            ! column loop
            ! N deposition and fixation
            sminn_vr(j,i) = sminn_vr(j,i) &
                          + fert_to_sminn(i) * deltim * ndep_prof(j,i)
         ENDDO
         IF(DEF_USE_CNSOYFIXN)THEN
            DO j = 1, nl_soil
               sminn_vr(j,i) = sminn_vr(j,i) &
                             + soyfixn_to_sminn(i) * deltim * nfixation_prof(j,i)
            ENDDO
         ENDIF
      ELSE
         DO j = 1, nl_soil
            ! N deposition and fixation (put all into NH4 pool)
            smin_nh4_vr(j,i) = smin_nh4_vr(j,i) &
                             + fert_to_sminn(i) * deltim * ndep_prof(j,i)
         ENDDO
         IF(DEF_USE_CNSOYFIXN)THEN
            DO j = 1, nl_soil
               smin_nh4_vr(j,i) = smin_nh4_vr(j,i) &
                                + soyfixn_to_sminn(i) * deltim * nfixation_prof(j,i)
            ENDDO
         ENDIF
      ENDIF
#endif

    ! decomposition fluxes
      DO k = 1, ndecomp_transitions
         DO j = 1, nl_soil
            decomp_npools_sourcesink(j,donor_pool(k),i) = &
                      decomp_npools_sourcesink(j,donor_pool(k),i) - &
                      decomp_ntransfer_vr(j,k,i) * deltim
         ENDDO
      ENDDO


      DO k = 1, ndecomp_transitions
         IF ( receiver_pool(k) /= 0 ) THEN  ! skip terminal transitions
            DO j = 1, nl_soil
               decomp_npools_sourcesink(j,receiver_pool(k),i) = &
                      decomp_npools_sourcesink(j,receiver_pool(k),i) + &
                         (decomp_ntransfer_vr(j,k,i) + &
                          decomp_sminn_flux_vr(j,k,i)) * deltim
            ENDDO
         ELSE  ! terminal transitions
            DO j = 1, nl_soil
               decomp_npools_sourcesink(j,donor_pool(k),i) = &
                      decomp_npools_sourcesink(j,donor_pool(k),i) - &
                         decomp_sminn_flux_vr(j,k,i) * deltim
            ENDDO
         ENDIF
      ENDDO

      IF(DEF_USE_SASU .or. DEF_USE_DiagMatrix)THEN
         DO j = 1, nl_soil
            AKX_met_to_soil1_n_vr_acc  (j,i) = AKX_met_to_soil1_n_vr_acc  (j,i) + (decomp_ntransfer_vr(j, 1,i) + decomp_sminn_flux_vr(j, 1,i)) * deltim
            AKX_cel_to_soil1_n_vr_acc  (j,i) = AKX_cel_to_soil1_n_vr_acc  (j,i) + (decomp_ntransfer_vr(j, 2,i) + decomp_sminn_flux_vr(j, 2,i)) * deltim
            AKX_lig_to_soil2_n_vr_acc  (j,i) = AKX_lig_to_soil2_n_vr_acc  (j,i) + (decomp_ntransfer_vr(j, 3,i) + decomp_sminn_flux_vr(j, 3,i)) * deltim
            AKX_soil1_to_soil2_n_vr_acc(j,i) = AKX_soil1_to_soil2_n_vr_acc(j,i) + (decomp_ntransfer_vr(j, 4,i) + decomp_sminn_flux_vr(j, 4,i)) * deltim
            AKX_cwd_to_cel_n_vr_acc    (j,i) = AKX_cwd_to_cel_n_vr_acc    (j,i) + (decomp_ntransfer_vr(j, 5,i) + decomp_sminn_flux_vr(j, 5,i)) * deltim
            AKX_cwd_to_lig_n_vr_acc    (j,i) = AKX_cwd_to_lig_n_vr_acc    (j,i) + (decomp_ntransfer_vr(j, 6,i) + decomp_sminn_flux_vr(j, 6,i)) * deltim
            AKX_soil1_to_soil3_n_vr_acc(j,i) = AKX_soil1_to_soil3_n_vr_acc(j,i) + (decomp_ntransfer_vr(j, 7,i) + decomp_sminn_flux_vr(j, 7,i)) * deltim
            AKX_soil2_to_soil1_n_vr_acc(j,i) = AKX_soil2_to_soil1_n_vr_acc(j,i) + (decomp_ntransfer_vr(j, 8,i) + decomp_sminn_flux_vr(j, 8,i)) * deltim
            AKX_soil2_to_soil3_n_vr_acc(j,i) = AKX_soil2_to_soil3_n_vr_acc(j,i) + (decomp_ntransfer_vr(j, 9,i) + decomp_sminn_flux_vr(j, 9,i)) * deltim
            AKX_soil3_to_soil1_n_vr_acc(j,i) = AKX_soil3_to_soil1_n_vr_acc(j,i) + (decomp_ntransfer_vr(j,10,i) + decomp_sminn_flux_vr(j,10,i)) * deltim

            AKX_met_exit_n_vr_acc      (j,i) = AKX_met_exit_n_vr_acc      (j,i) + decomp_ntransfer_vr(j, 1,i) * deltim
            AKX_cel_exit_n_vr_acc      (j,i) = AKX_cel_exit_n_vr_acc      (j,i) + decomp_ntransfer_vr(j, 2,i) * deltim
            AKX_lig_exit_n_vr_acc      (j,i) = AKX_lig_exit_n_vr_acc      (j,i) + decomp_ntransfer_vr(j, 3,i) * deltim
            AKX_soil1_exit_n_vr_acc    (j,i) = AKX_soil1_exit_n_vr_acc    (j,i) + decomp_ntransfer_vr(j, 4,i) * deltim
            AKX_cwd_exit_n_vr_acc      (j,i) = AKX_cwd_exit_n_vr_acc      (j,i) + decomp_ntransfer_vr(j, 5,i) * deltim
            AKX_cwd_exit_n_vr_acc      (j,i) = AKX_cwd_exit_n_vr_acc      (j,i) + decomp_ntransfer_vr(j, 6,i) * deltim
            AKX_soil1_exit_n_vr_acc    (j,i) = AKX_soil1_exit_n_vr_acc    (j,i) + decomp_ntransfer_vr(j, 7,i) * deltim
            AKX_soil2_exit_n_vr_acc    (j,i) = AKX_soil2_exit_n_vr_acc    (j,i) + decomp_ntransfer_vr(j, 8,i) * deltim
            AKX_soil2_exit_n_vr_acc    (j,i) = AKX_soil2_exit_n_vr_acc    (j,i) + decomp_ntransfer_vr(j, 9,i) * deltim
            AKX_soil3_exit_n_vr_acc    (j,i) = AKX_soil3_exit_n_vr_acc    (j,i) + decomp_ntransfer_vr(j,10,i) * deltim
         ENDDO
      ENDIF

      IF(.not. DEF_USE_NITRIF)THEN

           !--------------------------------------------------------
           !-------------    NITRIF_DENITRIF OFF -------------------
           !--------------------------------------------------------

           ! immobilization/mineralization in litter-to-SOM and SOM-to-SOM fluxes and denitrification fluxes
         DO k = 1, ndecomp_transitions
            IF ( receiver_pool(k) /= 0 ) THEN  ! skip terminal transitions
               DO j = 1, nl_soil
                  sminn_vr(j,i)  = sminn_vr(j,i) - &
                                  (sminn_to_denit_decomp_vr(j,k,i) + &
                                  decomp_sminn_flux_vr(j,k,i))* deltim
               ENDDO
            ELSE
               DO j = 1, nl_soil
                  sminn_vr(j,i)  = sminn_vr(j,i) - &
                                  sminn_to_denit_decomp_vr(j,k,i)* deltim

                  sminn_vr(j,i)  = sminn_vr(j,i) + &
                                  decomp_sminn_flux_vr(j,k,i)* deltim

               ENDDO
            ENDIF
         ENDDO


         DO j = 1, nl_soil
            ! "bulk denitrification"
            sminn_vr(j,i) = sminn_vr(j,i) - sminn_to_denit_excess_vr(j,i) * deltim

            ! total plant uptake from mineral N
            sminn_vr(j,i) = sminn_vr(j,i) - sminn_to_plant_vr(j,i)*deltim
            ! flux that prevents N limitation (when Carbon_only is set)
            sminn_vr(j,i) = sminn_vr(j,i) + supplement_to_sminn_vr(j,i)*deltim
         ENDDO

      ELSE

           !--------------------------------------------------------
           !-------------    NITRIF_DENITRIF ON --------------------
           !--------------------------------------------------------

         DO j = 1, nl_soil

            ! mineralization fluxes (divert a fraction of this stream to nitrification flux, add the rest to NH4 pool)
            smin_nh4_vr(j,i) = smin_nh4_vr(j,i) + gross_nmin_vr(j,i)*deltim

            ! immobilization fluxes
            smin_nh4_vr(j,i) = smin_nh4_vr(j,i) - actual_immob_nh4_vr(j,i)*deltim

            smin_no3_vr(j,i) = smin_no3_vr(j,i) - actual_immob_no3_vr(j,i)*deltim

            ! plant uptake fluxes
            smin_nh4_vr(j,i) = smin_nh4_vr(j,i) - smin_nh4_to_plant_vr(j,i)*deltim

            smin_no3_vr(j,i) = smin_no3_vr(j,i) - smin_no3_to_plant_vr(j,i)*deltim


            ! Account for nitrification fluxes
            smin_nh4_vr(j,i) = smin_nh4_vr(j,i) - f_nit_vr(j,i) * deltim

            smin_no3_vr(j,i) = smin_no3_vr(j,i) + f_nit_vr(j,i) * deltim &
                             * (1._r8 - nitrif_n2o_loss_frac)

            ! Account for denitrification fluxes
            smin_no3_vr(j,i) = smin_no3_vr(j,i) - f_denit_vr(j,i) * deltim

            ! flux that prevents N limitation (when Carbon_only is set; put all into NH4)
            smin_nh4_vr(j,i) = smin_nh4_vr(j,i) + supplement_to_sminn_vr(j,i)*deltim

            ! update diagnostic total
            sminn_vr(j,i) = smin_nh4_vr(j,i) + smin_no3_vr(j,i)

         ENDDO
      ENDIF

   END SUBROUTINE SoilBiogeochemNStateUpdate1

   SUBROUTINE soil_acidification_buffering (i,ps,pe,deltim,nl_soil, dz_soi, z_soi)
   
   IMPLICIT NONE

   integer ,intent(in) :: i                  ! patch index
   integer ,intent(in) :: ps                 ! start pft index
   integer ,intent(in) :: pe                 ! end pft index
   real(r8),intent(in) :: deltim             ! time step in seconds
   integer ,intent(in) :: nl_soil            ! number of total soil layers
   real(r8),intent(in) :: dz_soi(1:nl_soil)  ! thicknesses of each soil layer (m)
   real(r8),intent(in) :: z_soi (1:nl_soil)   ! depth of each soil layer (m)

   real(r8) :: eff_porosity, vol_ice, vol_liq(1:nl_soil)
   integer  :: j, m, newton_iter, max_iter = 50
   logical  :: converged, do_chemistry(1:nl_soil)

   real(r8) :: factor_co2(1:nl_soil), factor_cec(1:nl_soil), factor_cal(1:nl_soil)
   real(r8) :: M_Ca = 0.04_r8, M_CaCO3 = 0.1_r8  ! kg/mol
   real(r8) :: ca_in_caco3 = 0.4_r8 ! soluble Ca in CaCO3
   real(r8) :: ss_vliq_min = 5e-2_r8, threshold_calcite_mass = 1e-14_r8, sigmoid_slope = 0.5_r8
   real(r8) :: TOL_R_SHALLOW = 1.0E-7_r8, TOL_D_SHALLOW = 1.0E-3_r8
   real(r8) :: TOL_R_DEEP = 1.0E-8_r8, TOL_D_DEEP = 1.0E-3_r8
   real(r8) :: gN_to_mol = 1.0_r8 / 14.0_r8  ! convert gN to moles of N
   real(r8) :: rgas = 8.314467591                 ! universal gas constant (J mol-1 K-1)

   ! --- CO2 exchange ---
   real(r8) :: CO2_exchange(1:nl_soil), h2co3_soil(1:nl_soil)
   real(r8) :: rate_CO2_explicit, CO2_exchange_amount, DIC_after_CO2
   
   ! --- Cation exchange ---
   real(r8) :: rate_EBc, base_iter, base_final, base_start
   real(r8) :: Bc_ex_flux, Bc_ex_iter, Bc_ex_start, H_ex_flux, H_ex_iter, H_ex_start, H_ex_final, Bct_ex_flux, Bct_ex_iter, Bct_ex_start, Bct_ex_final
   real(r8) :: dEBc_dh, dEBc_dBc, dH_ex_flux_dBct, dH_ex_flux_dh, dBc_ex_flux_dBc, dBc_ex_flux_dh, dBc_ex_flux_dca, dBc_ex_flux_dBct, dca_ex_flux_dBc, dca_ex_flux_dh, dca_ex_flux_dca
   real(r8) :: dBct_ex_flux_dh, dBct_ex_flux_dBct, dEBc_dBct
   real(r8) :: CEC_total, ca_ex_flux, weight_cec, K_HBc(1:nl_soil)

   ! --- Calcite buffer ---
   real(r8) :: rate_CaCO3 = 0.0_r8, k_calcite = 0.0_r8, Q_calcite = 0.0_r8, saturation_ratio = 0.0_r8, k0_cal = 0.0_r8, k_H_cal = 0.0_r8, k_OH_cal = 0.0_r8
   real(r8) :: K_H(1:nl_soil), K1(1:nl_soil), K2(1:nl_soil), Kw(1:nl_soil), log_KH, log_K1, log_K2, log_Kw
   real(r8) :: Ksp_calcite(1:nl_soil), log_Ksp_calcite, weight_cal, A_calcite(1:nl_soil), depth_factor(1:nl_soil)
   real(r8) :: init_calcite_content(10) = [5.0_r8,5.0_r8,5.0_r8,5.0_r8,5.0_r8,5.0_r8,8.0_r8,8.0_r8,8.0_r8,15.0_r8] ! US-Ne1 CaCO3 content (%)
   
   ! --- Gibbsite buffer ---
   real(r8) :: Kal_1(1:nl_soil), Kal_2(1:nl_soil), Kal_3(1:nl_soil), log_Kal_1, log_Kal_2, log_Kal_3
   real(r8) :: Ksp_gibbsite(1:nl_soil)
   real(r8) :: pH_iter, al_weight, k_sigmoid = 5.0_r8, al_alk_iter_raw, dal_weight_dh, dal_alk_dh_raw, pH_start

   ! newton iteration
   real(r8) :: StateVector(4), Residuals(4), Jacobian(4,4), DeltaX(4)
   real(r8) :: dF1_dh, dF2_dh, dF3_dh,dF4_dh, dh2co3_dh, dhco3_dh, dco3_dh, doh_dh, drate_cal_dh, drate_EBc_dh, dal3_dh, dal2_dh, dal_dh, dal_alk_dh
   real(r8) :: dF1_dlogH, dF2_dlogH, dF3_dlogH, dF4_dlogH, dF1_dlogDIC, dF2_dlogDIC, dF3_dlogDIC, dF4_dlogDIC
   real(r8) :: dF1_dca, dF2_dca, dF3_dca, dF4_dca, drate_cal_dca, dF1_dlogca, dF2_dlogca, dF3_dlogca, dF4_dlogca
   real(r8) :: dF1_dDIC, dF2_dDIC, dF3_dDIC, dF4_dDIC, dh2co3_dDIC, dhco3_dDIC, dco3_dDIC, drate_cal_dDIC
   real(r8) :: dF1_dBc, dF2_dBc, dF3_dBc, dF4_dBc, dF1_dlogBc, dF2_dlogBc, dF3_dlogBc, dF4_dlogBc
   real(r8) :: ANC_start, h_start, oh_start, ca_start, DIC_start, Bc_start, Bct_start, co3_start, h2co3_start, hco3_start, A_start, al3_start, al2_start, al_start, al_alk_start  ! mol/l
   real(r8) :: ANC_final, h_final, oh_final, ca_final, DIC_final, Bc_final, Bct_final, co3_final, h2co3_final, hco3_final, A_final, al3_final, al2_final, al_final  ! mol/l
   real(r8) :: ANC_iter, h_iter, oh_iter, ca_iter, DIC_iter, Bc_iter, Bct_iter, co3_iter, h2co3_iter, hco3_iter, A_iter, al3_iter, al2_iter, al_iter, al_alk_iter  ! mol/l
   real(r8) :: logH_iter, logDIC_iter, logca_iter, logBc_iter, dh_dlogH, dca_dlogCa, dBc_dlogBc, dDIC_dlogDIC
   real(r8) :: residual_norm, residual_norm_sq, delta_norm, delta_norm_sq, prev_residual_norm
   real(r8) :: nr_tol_residual, nr_tol_delta

   ! Time step reduction mechanism
   real(r8) :: step_scale, min_dt_factor = 0.1_r8, dt_reduction_factor = 0.5_r8
   real(r8) :: cum_rate_CO2, cum_rate_cal, cum_rate_Bcex, cum_rate_Hex
   integer :: stagnation_count, min_step_count, residual_div_count
   logical :: success
   real(r8) :: backup_E_Bc, backup_calcite_mass, backup_calcite_content
   real(r8) :: convert = 1.957e-3_r8, neu_const, acid_const, alk_const
   real(r8) :: k_ref_co2  ! CO2 gas-liquid exchange rate constant at surface [s-1] （Johnson, K.S. 1982)
   real(r8) :: z_ref_co2  ! depth decay scale for CO2 gas-liquid exchange [m]

      ! assign tunable parameters from namelist (no recompile needed)
      neu_const  = DEF_SOILACID_neu_const
      acid_const = DEF_SOILACID_acid_const
      alk_const  = DEF_SOILACID_alk_const
      k_ref_co2  = DEF_SOILACID_k_ref_co2
      z_ref_co2  = DEF_SOILACID_z_ref_co2

      DO j = 1, nl_soil
         factor_co2(j) = 1.0_r8 - EXP(-k_ref_co2 * EXP(-z_soi(j) / z_ref_co2) * deltim)
         ! write(*,*) 'layer =', j, ' k_ref_co2 =', k_ref_co2, ' factor_co2 =', factor_co2(j)
         factor_cec(j) = 1.0_r8
         factor_cal(j) = 1.0_r8
         IF(calcite_mass(j,i) <= 0.0_r8) THEN
            calcite_mass(j,i) = MAX(0.0_r8, (20.0_r8/100.0_r8) * BD_all(j,i) * dz_soi(j))
         ENDIF
         IF(t_soisno(j,i) > tfrz - 10.0_r8) THEN
            CALL calc_equilibrium_constants(t_soisno(j,i), K_H(j), K1(j), K2(j), Ksp_calcite(j), Kw(j), &
                                             Kal_1(j), Kal_2(j), Kal_3(j), Ksp_gibbsite(j))
         ELSE
            K_H(j) = 1.0e-2_r8; K1(j) = 1.0e-7_r8; K2(j) = 1.0e-11_r8
            Ksp_calcite(j) = 1.0e-9_r8; Kw(j) = 1.0e-14_r8
            Kal_1(j) = 1.0e-9_r8; Kal_2(j) = 1.0e-18_r8; Kal_3(j) = 1.0e-27_r8
            Ksp_gibbsite(j) = 1.3e-33_r8
         ENDIF
         vol_ice = min(porsl(j,i), wice_soisno(j,i)/(dz_soi(j)*denice))
         eff_porosity = max(0.01, porsl(j,i)-vol_ice)
         vol_liq(j) = min(eff_porosity, wliq_soisno(j,i)/(dz_soi(j)*denh2o))
         do_chemistry(j) = (vol_liq(j) >= ss_vliq_min .and. t_soisno(j,i) > tfrz - 10.0_r8)
         depth_factor(j) = exp(-z_soi(j)/5.0_r8)
         A_calcite(j) = (8_r8*wf_clay(j,i) + 0.3_r8*wf_sand(j,i) + 2.2_r8*(1_r8-wf_clay(j,i)-wf_sand(j,i))) * BD_all(j,i)/(max(vol_liq(j), 1e-20_r8)*1000_r8)*0.08 
         K_HBc(j) = 10**(3.6_r8*wf_clay(j,i) + 3.2_r8*(1_r8-wf_clay(j,i)))
      ENDDO

      ! Nitrogen cycle
      DO j = 1, nl_soil
         IF(do_chemistry(j)) THEN 
            ! H+ production and consumption from N-cycle（mol/m3/s）
            h_content_prod_vr(j,i) = 0.0_r8
            h_content_cons_vr(j,i) = 0.0_r8

            ! Mineralization and dissolution (RNH2 + H+ + H2O -> ROH + NH4+)
            h_content_cons_vr(j,i) = h_content_cons_vr(j,i) + gross_nmin_vr(j,i) * gN_to_mol 
            
            ! NO3- uptake (ROH + NO3- + H+ + 2CH2O -> RNH2 + 2CO2 + 2H2O)
            h_content_cons_vr(j,i) =  h_content_cons_vr(j,i) + smin_no3_to_plant_vr(j,i) * gN_to_mol 
            
            ! Incomplete nitration (2NO2- + 6H+ + 4e- -> N2O + 3H2O)
            h_content_cons_vr(j,i) = h_content_cons_vr(j,i)  + 3.0_r8 * f_nit_vr(j,i) * nitrif_n2o_loss_frac * gN_to_mol 

            ! Denitrification (5CH20 + 4NO3- + 4H+ -> 2N2 + 5CO2 + 7H2O)
            h_content_cons_vr(j,i) =  h_content_cons_vr(j,i) + f_denit_vr(j,i) * gN_to_mol
               
            ! NH4+ immobilization (ROH + NH4+ -> RNH2 + H+ + H2O)
            h_content_prod_vr(j,i) = h_content_prod_vr(j,i) + actual_immob_nh4_vr(j,i) * gN_to_mol 
               
            ! NH4+ uptake (NH4+ + ROH -> RNH2 + H+ + H2O)
            h_content_prod_vr(j,i) =  h_content_prod_vr(j,i) + smin_nh4_to_plant_vr(j,i) * gN_to_mol 
            
            ! Ammoxidation (NH4+ + 1.5O2 -> NO2- + 2H+ + H2O) + (NO2- + 0.5O2 -> NO3-)
            h_content_prod_vr(j,i) = h_content_prod_vr(j,i)  + 2.0_r8 * f_nit_vr(j,i) * gN_to_mol 

            ! H+ net production from N-cycle 
            h_from_N(j,i) = (h_content_prod_vr(j,i) - h_content_cons_vr(j,i)) / (max(vol_liq(j), 1e-20_r8) * 1000.0_r8)  ! mol/m3-soil/s → mol/l-water/s
         ENDIF

      ENDDO

      CALL Gas_liquid_CO2_exchange_M2(i, ps, pe, nl_soil, z_soi, dz_soi, vol_liq, K1, K2, K_H, convert) ! for pCO2_soil
      ! CALL Gas_liquid_CO2_exchange_M1(i, ps, pe, nl_soil, z_soi, dz_soi, vol_liq, K1, K2, K_H, convert) ! for pCO2_soil

      DO j = 1, nl_soil  
         IF (.not.do_chemistry(j)) CYCLE
         backup_E_Bc = E_Bc(j,i)
         backup_calcite_mass = calcite_mass(j,i)
         backup_calcite_content = calcite_content(j,i)
         
         IF (pH_init_vr(j,i) <= 6_r8) THEN
            weight_cal = 0.0_r8
         ELSE
            weight_cal = 1.0_r8 / (1.0_r8 + EXP(-(pH_init_vr(j,i) - 7.0_r8) / sigmoid_slope))
            weight_cal = MIN(MAX(weight_cal, 0.0_r8), 1.0_r8)
         ENDIF
         weight_cec = 1.0_r8 - weight_cal

         IF (j <= 5) THEN
            nr_tol_residual = TOL_R_SHALLOW
            nr_tol_delta    = TOL_D_SHALLOW
         ELSE
            nr_tol_residual = TOL_R_DEEP
            nr_tol_delta    = TOL_D_DEEP
         ENDIF

         h2co3_start = (h_conc_vr(j,i)**2)/(h_conc_vr(j,i)**2 + K1(j)*h_conc_vr(j,i) + K1(j)*K2(j)) * DIC_conc_vr(j,i)
         h2co3_soil(j) = pCO2_soil(j,i)/1e6_r8 * K_H(j) ! mol/l
         IF (vol_liq(j) >= ss_vliq_min) THEN
            rate_CO2_explicit = factor_co2(j) * (h2co3_soil(j) - h2co3_start)/deltim   ! mol/l/s
         ELSE
            rate_CO2_explicit = 0.0_r8
         ENDIF
         CO2_exchange_amount = rate_CO2_explicit * deltim
         DIC_after_CO2 = DIC_conc_vr(j,i) + CO2_exchange_amount

         ! --- The rate constant of calcium carbonate mineralization A: Xu, 2012 ---
         k0_cal = neu_const * exp((-2.35e4_r8/rgas)*(1.0_r8/t_soisno(j,i) - 1.0_r8/298.15))
         k_H_cal = acid_const * exp((-1.44e4_r8/rgas)*(1.0_r8/t_soisno(j,i) - 1.0_r8/298.15))
         k_OH_cal = alk_const * exp((-3.54e4_r8/rgas)*(1.0_r8/t_soisno(j,i) - 1.0_r8/298.15))

         IF (wliq_soisno(j,i) > 1e-10_r8) THEN
            CEC_total = CEC(j,i) * 0.01_r8 * BD_all(j,i) * dz_soi(j)/wliq_soisno(j,i)  ! mol(+)/l
         ELSE
            CEC_total = 0.0_r8
         ENDIF

         ! --- Multivariable Newton-Laverson method ---
         h_start = MIN(MAX(h_conc_vr(j,i), 1.0e-14_r8), 1.0_r8)
         oh_start = MIN(MAX(oh_conc_vr(j,i), 1.0e-14_r8), 1.0_r8)
         DIC_start = MAX(DIC_after_CO2, 1.0e-10_r8)
         ca_start = MAX(ca_conc_vr(j,i), 1.0e-10_r8)
         Bc_start = MAX(Bc_conc_vr(j,i), 1.0e-10_r8)  ! K+Mg
         
         IF(pH_init_vr(j,i) <= 6_r8) THEN
            Bct_start = ca_start + Bc_start
            IF (Bct_start > 1.0e-10_r8) THEN
               base_start = SQRT(MAX(Bct_start, 1.0e-20_r8)) / &
                           (SQRT(MAX(Bct_start, 1.0e-20_r8)) + K_HBc(j) * h_start)
               base_start = MIN(MAX(base_start, 1.0e-14_r8), 1.0_r8)
            ELSE
               base_start = MIN(MAX(backup_E_Bc/100_r8, 1.0e-14_r8), 1.0_r8)
            ENDIF
         ELSE
            base_start = MIN(MAX(backup_E_Bc/100_r8, 1.0e-14_r8), 1.0_r8)
         ENDIF

         Bct_ex_start = base_start/2 * CEC_total ! mol/l
         H_ex_start = (1-base_start) * CEC_total
         
         ! Fast reaction of carbonic acid and aluminum system
         A_start = h_start**2 + K1(j)*h_start + K1(j)*K2(j)
         h2co3_start = (h_start**2/A_start) * DIC_start
         hco3_start = (K1(j)*h_start/A_start) * DIC_start
         co3_start = (K1(j)*K2(j)/A_start) * DIC_start

         al3_start = MIN(Ksp_gibbsite(j) / (oh_start**3), 1.0_r8)
         al2_start = al3_start * Kal_2(j)/(h_start**2)
         al_start = al3_start * Kal_1(j)/h_start
         al_alk_start = 3.0_r8*al3_start + 2.0_r8*al_start + al2_start

         ! IF (-log10(h_start) > 5.5_r8) THEN
         !    al_alk_start = 0.0_r8  ! 中性和碱性条件下忽略铝
         ! ENDIF
         !=========================================================================
         pH_start     = -LOG10(h_start)
         al_weight    = 1.0_r8 / (1.0_r8 + EXP(k_sigmoid * (pH_start - 5.5_r8)))
         al_alk_start = al_alk_start * al_weight
         !=========================================================================

         ! ANC = [HCO3-] + 2[CO3-] + [OH-] - [H+] - 3[Al3+] - 2[AlOH2+] - [Al(OH)2+]
         ANC_start = hco3_start + 2.0_r8*co3_start + oh_start - h_start - al_alk_start
         
         StateVector(1) = LOG10(MAX(h_start, 1.0e-14_r8))
         StateVector(2) = LOG10(MAX(ca_start, 1.0e-10_r8))
         StateVector(3) = LOG10(MAX(DIC_start, 1.0e-10_r8))
         StateVector(4) = LOG10(MAX(Bc_start, 1.0e-10_r8))

         converged = .false.
         prev_residual_norm = 1e20_r8  ! 初始化为很大的值
         stagnation_count = 0  ! 停滞检测计数器
         min_step_count = 0  ! 最小步长计数器
         residual_div_count = 0  ! 残差发散计数器

         DO newton_iter = 1, max_iter

            logH_iter = MIN(MAX(StateVector(1), -14.0_r8), 0.0_r8)
            h_iter = MIN(MAX(10.0_r8 ** logH_iter, 1.0e-14_r8), 1.0_r8)
            oh_iter = MIN(MAX(Kw(j) / h_iter, 1.0e-14_r8), 1.0_r8)

            logCa_iter = MIN(MAX(StateVector(2), -10.0_r8), 2.0_r8)
            ca_iter = MIN(MAX(10.0_r8 ** logCa_iter, 1.0e-10_r8), 100.0_r8)
            logDIC_iter = MIN(MAX(StateVector(3), -10.0_r8), 2.0_r8)
            DIC_iter = MIN(MAX(10.0_r8 ** logDIC_iter, 1.0e-10_r8), 100.0_r8)
            logBc_iter = MIN(MAX(StateVector(4), -10.0_r8), 2.0_r8)
            Bc_iter = MIN(MAX(10.0_r8 ** logBc_iter, 1.0e-10_r8), 100.0_r8)

            A_iter = h_iter**2 + K1(j)*h_iter + K1(j)*K2(j)
            h2co3_iter = (h_iter**2/A_iter) * DIC_iter
            hco3_iter = (K1(j)*h_iter / A_iter) * DIC_iter
            co3_iter = (K1(j)*K2(j) / A_iter) * DIC_iter

            al3_iter = MIN(Ksp_gibbsite(j) / (oh_iter**3), 1.0_r8)
            al2_iter = al3_iter * Kal_2(j)/(h_iter**2)
            al_iter = al3_iter * Kal_1(j)/h_iter
            al_alk_iter_raw = 3.0_r8*al3_iter + 2.0_r8*al_iter + al2_iter
                 
            ! ! Al缓冲 硬切换
            ! IF (-log10(h_iter) > 5.5_r8) THEN
            !    al_alk_iter = 0.0_r8
            ! ENDIF

            !=======================================================================
            ! Al缓冲 sigmoid连续过渡
            pH_iter   = -LOG10(h_iter)
            al_weight = 1.0_r8 / (1.0_r8 + EXP(k_sigmoid * (pH_iter - 5.5_r8)))
            al_alk_iter = al_alk_iter_raw * al_weight   ! 连续过渡替代硬开关
            !=======================================================================

            ANC_iter = hco3_iter + 2.0_r8*co3_iter + oh_iter - h_iter - al_alk_iter

            IF (calcite_mass(j,i) > threshold_calcite_mass) THEN
               Q_calcite = ca_iter * co3_iter
               saturation_ratio = Q_calcite / Ksp_calcite(j)
               k_calcite = k0_cal + k_H_cal * h_iter + k_OH_cal * oh_iter
               rate_CaCO3 = weight_cal * factor_cal(j) * k_calcite * A_calcite(j) * (1.0_r8 - saturation_ratio)
            ELSE
               rate_CaCO3 = 0.0_r8
            ENDIF
            
            Bct_iter = ca_iter + Bc_iter
            IF (Bct_iter > 1.0e-10_r8 .AND. ca_iter > 1.0e-10_r8 .AND. Bc_iter > 1.0e-10_r8) THEN
                base_iter = SQRT(MAX(Bct_iter, 1.0e-20_r8)) / (SQRT(MAX(Bct_iter, 1.0e-20_r8)) + K_HBc(j) * h_iter)
                Bct_ex_iter = base_iter/2 * CEC_total
                H_ex_iter = (1-base_iter) * CEC_total
                H_ex_flux = weight_cec * factor_cec(j) * (H_ex_iter - H_ex_start) / deltim 
                Bct_ex_flux = weight_cec * factor_cec(j) * (Bct_ex_iter - Bct_ex_start) / deltim 
                ca_ex_flux = Bct_ex_flux * (ca_iter / Bct_iter)
                Bc_ex_flux = Bct_ex_flux * (Bc_iter / Bct_iter)
            ELSE
                base_iter = base_start
                H_ex_flux = 0.0_r8
                Bct_ex_flux = 0.0_r8
                ca_ex_flux = 0.0_r8
                Bc_ex_flux = 0.0_r8
            ENDIF
                        
            Residuals(1) = (ANC_iter - ANC_start)/deltim  - 2.0_r8*rate_CaCO3 - H_ex_flux + h_from_N(j,i)
            Residuals(2) = (ca_iter - ca_start)/deltim  - rate_CaCO3 + ca_ex_flux
            Residuals(3) = (DIC_iter - DIC_start)/deltim  - rate_CaCO3
            Residuals(4) = (Bc_iter - Bc_start)/deltim  + Bc_ex_flux

            ! Jacobian matrix calculation
            dh_dlogH = h_iter * log(10.0_r8)
            doh_dh = - Kw(j)/(h_iter**2.0_r8)
            dca_dlogCa = ca_iter * log(10.0_r8)
            dDIC_dlogDIC = DIC_iter * log(10.0_r8)
            dBc_dlogBc = Bc_iter * log(10.0_r8)

            dh2co3_dh = DIC_iter * (2.0_r8*h_iter/A_iter - (h_iter**2.0_r8) * (2.0_r8*h_iter + K1(j))/(A_iter**2.0_r8))
            dh2co3_dDIC = h_iter**2.0_r8/A_iter
            dhco3_dh = DIC_iter * (K1(j)/A_iter -  K1(j)*h_iter*(2.0_r8*h_iter + K1(j))/(A_iter**2.0_r8))
            dhco3_dDIC = K1(j)*h_iter/A_iter
            dco3_dh = -DIC_iter * (K1(j)*K2(j)*(2.0_r8*h_iter + K1(j))/(A_iter**2.0_r8))
            dco3_dDIC = K1(j)*K2(j)/A_iter

            ! IF (-log10(h_iter) > 5.5_r8) THEN
            !    dal_alk_dh = 0.0_r8
            ! ELSE
            !    dal3_dh = - 3.0_r8 * Ksp_gibbsite(j)/(oh_iter**4.0_r8) * doh_dh
            !    dal2_dh = dal3_dh * Kal_2(j)/(h_iter**2.0_r8)- 2.0_r8 * al3_iter * Kal_2(j)/(h_iter**3)
            !    dal_dh = dal3_dh * Kal_1(j)/h_iter - al3_iter * Kal_1(j)/(h_iter**2)
            !    dal_alk_dh = 3.0_r8 * dal3_dh + 2.0_r8 * dal_dh + dal2_dh
            ! ENDIF

            !=======================================================================
            dal3_dh = - 3.0_r8 * Ksp_gibbsite(j)/(oh_iter**4.0_r8) * doh_dh
            dal2_dh = dal3_dh * Kal_2(j)/(h_iter**2.0_r8)- 2.0_r8 * al3_iter * Kal_2(j)/(h_iter**3)
            dal_dh = dal3_dh * Kal_1(j)/h_iter - al3_iter * Kal_1(j)/(h_iter**2)
            dal_alk_dh_raw = 3.0_r8 * dal3_dh + 2.0_r8 * dal_dh + dal2_dh

            dal_weight_dh = al_weight * (1.0_r8 - al_weight) * k_sigmoid / (h_iter * LOG(10.0_r8))
            dal_alk_dh = dal_alk_dh_raw * al_weight + al_alk_iter_raw * dal_weight_dh
            !=======================================================================

            drate_cal_dh = factor_cal(j) * A_calcite(j) * ((k_H_cal + k_OH_cal * doh_dh) * (1.0_r8 - (ca_iter*co3_iter)/Ksp_calcite(j)) + k_calcite * (-ca_iter/Ksp_calcite(j))*dco3_dh)
            drate_cal_dca = factor_cal(j) * (k_calcite * A_calcite(j) * (-co3_iter/Ksp_calcite(j)))
            drate_cal_dDIC = factor_cal(j) * (k_calcite * A_calcite(j) * (-ca_iter/Ksp_calcite(j))*dco3_dDIC)
   
            IF (Bct_iter > 1.0e-10_r8 .AND. ca_iter > 1.0e-10_r8 .AND. Bc_iter > 1.0e-10_r8) THEN
               dEBc_dh = (-sqrt(MAX(Bct_iter, 1.0e-20_r8))*K_HBc(j))/((sqrt(MAX(Bct_iter, 1.0e-20_r8))+K_HBc(j)*h_iter)**2)
               dEBc_dBct = (K_HBc(j)*h_iter)/(2*sqrt(MAX(Bct_iter, 1.0e-20_r8))*(sqrt(MAX(Bct_iter, 1.0e-20_r8))+K_HBc(j)*h_iter)**2)
               
               dBct_ex_flux_dh = factor_cec(j) * dEBc_dh * CEC_total / (2.0_r8 * deltim )
               dBct_ex_flux_dBct = factor_cec(j) *  dEBc_dBct * CEC_total / (2.0_r8 * deltim )

               dH_ex_flux_dh = - factor_cec(j) * dEBc_dh * CEC_total / deltim 
               dH_ex_flux_dBct = - factor_cec(j) * dEBc_dBct * CEC_total / deltim 
               dca_ex_flux_dca = dBct_ex_flux_dBct * (ca_iter/Bct_iter) + Bct_ex_flux * (Bc_iter / (Bct_iter**2.0_r8))
               dca_ex_flux_dBc = dBct_ex_flux_dBct * (ca_iter/Bct_iter) + Bct_ex_flux * (-ca_iter / (Bct_iter**2.0_r8))
               dca_ex_flux_dh = dBct_ex_flux_dh * (ca_iter/Bct_iter)
               dBc_ex_flux_dca = dBct_ex_flux_dBct * (Bc_iter/Bct_iter) + Bct_ex_flux * (-Bc_iter / (Bct_iter**2.0_r8))
               dBc_ex_flux_dBc = dBct_ex_flux_dBct * (Bc_iter/Bct_iter) + Bct_ex_flux * (ca_iter / (Bct_iter**2.0_r8))
               dBc_ex_flux_dh = dBct_ex_flux_dh * (Bc_iter/Bct_iter)
            ELSE
               dH_ex_flux_dh = 0.0_r8
               dH_ex_flux_dBct = 0.0_r8
               dca_ex_flux_dca = 0.0_r8
               dca_ex_flux_dBc = 0.0_r8
               dca_ex_flux_dh = 0.0_r8
               dBc_ex_flux_dca = 0.0_r8
               dBc_ex_flux_dBc = 0.0_r8
               dBc_ex_flux_dh = 0.0_r8
            ENDIF

            ! Equation 1
            dF1_dh = (dhco3_dh + 2.0_r8*dco3_dh + doh_dh - 1.0_r8 - dal_alk_dh)/deltim  - 2.0_r8*weight_cal*drate_cal_dh - weight_cec*dH_ex_flux_dh
            dF1_dlogH = dF1_dh * dh_dlogH
            dF1_dca = -2.0_r8*weight_cal*drate_cal_dca - weight_cec*dH_ex_flux_dBct
            dF1_dlogca = dF1_dca * dca_dlogCa
            dF1_dDIC = (dhco3_dDIC + 2.0_r8*dco3_dDIC)/deltim  - 2.0_r8*weight_cal*drate_cal_dDIC
            dF1_dlogDIC = dF1_dDIC * dDIC_dlogDIC
            dF1_dBc = - weight_cec*dH_ex_flux_dBct
            dF1_dlogBc = dF1_dBc * dBc_dlogBc

            ! Equation 2
            dF2_dh = - weight_cal*drate_cal_dh + weight_cec*dca_ex_flux_dh
            dF2_dlogH = dF2_dh * dh_dlogH
            dF2_dca = 1.0_r8/deltim  - weight_cal*drate_cal_dca + weight_cec*dca_ex_flux_dca
            dF2_dlogca = dF2_dca * dca_dlogCa
            dF2_dDIC = - weight_cal*drate_cal_dDIC
            dF2_dlogDIC = dF2_dDIC * dDIC_dlogDIC
            dF2_dBc = weight_cec*dca_ex_flux_dBc
            dF2_dlogBc = dF2_dBc * dBc_dlogBc

            ! Equation 3
            dF3_dh = - weight_cal*drate_cal_dh
            dF3_dlogH = dF3_dh * dh_dlogH
            dF3_dca = - weight_cal*drate_cal_dca
            dF3_dlogca = dF3_dca * dca_dlogCa
            dF3_dDIC = 1.0_r8/deltim  - weight_cal*drate_cal_dDIC
            dF3_dlogDIC = dF3_dDIC * dDIC_dlogDIC
            dF3_dBc = 0.0_r8
            dF3_dlogBc = 0.0_r8

            ! Equation 4
            dF4_dh = weight_cec*dBc_ex_flux_dh
            dF4_dlogH = dF4_dh * dh_dlogH
            dF4_dca = weight_cec*dBc_ex_flux_dca
            dF4_dlogca = dF4_dca * dca_dlogCa
            dF4_dDIC = 0.0_r8
            dF4_dlogDIC = dF4_dDIC * dDIC_dlogDIC
            dF4_dBc = 1.0_r8/deltim + weight_cec*dBc_ex_flux_dBc
            dF4_dlogBc = dF4_dBc * dBc_dlogBc
            
            Jacobian(1,:) = [dF1_dlogH, dF1_dlogca, dF1_dlogDIC, dF1_dlogBc]
            Jacobian(2,:) = [dF2_dlogH, dF2_dlogca, dF2_dlogDIC, dF2_dlogBc]
            Jacobian(3,:) = [dF3_dlogH, dF3_dlogca, dF3_dlogDIC, dF3_dlogBc]
            Jacobian(4,:) = [dF4_dlogH, dF4_dlogca, dF4_dlogDIC, dF4_dlogBc]
            
            ! J·Δ = -Residuals 线性方程求解器
            CALL SolveLinearSystem(4, Jacobian, Residuals, DeltaX, success)
            IF (.NOT. success) THEN
                residual_norm_sq = SUM(Residuals**2)
                residual_norm = SQRT(MAX(residual_norm_sq, 1.0e-20_r8))
                 
                IF (residual_norm < nr_tol_residual) THEN
                  converged = .true.
                  EXIT
                ELSEIF (residual_norm < nr_tol_residual * 10.0_r8) THEN
                  converged = .true.
                  EXIT
                ELSE
                  converged = .false.
                  EXIT
                ENDIF
            ENDIF

            residual_norm_sq = SUM(Residuals**2)
            delta_norm_sq = SUM(DeltaX**2)

            ! --- 自适应牛顿阻尼 - 基于残差减少的线搜索 ---
            residual_norm = SQRT(MAX(residual_norm_sq, 1.0e-20_r8))
            step_scale = 1.0_r8
            IF (newton_iter > 1) THEN
                CALL AdaptiveDamping(StateVector, DeltaX, Residuals, residual_norm, nr_tol_residual, &
                                    step_scale, j, i, deltim , K1(j), K2(j), Kw(j), Ksp_calcite(j), &
                                    Kal_1(j), Kal_2(j), Kal_3(j), Ksp_gibbsite(j), CEC_total, &
                                    weight_cal, weight_cec, factor_cal(j), factor_cec(j), &
                                    k0_cal, k_H_cal, k_OH_cal, A_calcite(j), K_HBc(j), factor_co2(j), &
                                    h2co3_soil(j), ANC_start, ca_start, DIC_start, Bc_start, &
                                    H_ex_start, Bct_ex_start)
            ELSE
               IF (j <= 5) THEN
                  step_scale = 0.5_r8  ! 表面层更保守
               ELSE
                  step_scale = 0.8_r8  ! 深层稍微激进
               ENDIF
            ENDIF

            StateVector = StateVector + step_scale * DeltaX

            StateVector(1) = MIN(MAX(StateVector(1), -14.0_r8), 0.0_r8)
            StateVector(2) = MIN(MAX(StateVector(2), -10.0_r8), 2.0_r8)
            StateVector(3) = MIN(MAX(StateVector(3), -10.0_r8), 2.0_r8)
            StateVector(4) = MIN(MAX(StateVector(4), -10.0_r8), 2.0_r8)
            
            ! 检测连续最小步长（表明数值困难）
            IF (step_scale <= 1.0e-6_r8) THEN
               min_step_count = min_step_count + 1
               IF (min_step_count > 10) THEN
                  EXIT
               ENDIF
            ELSE
               min_step_count = MAX(0, min_step_count - 2)
            ENDIF
   
            delta_norm_sq = step_scale**2 * delta_norm_sq

            ! 主收敛条件：绝对收敛
            IF (residual_norm_sq < nr_tol_residual**2 .and. delta_norm_sq < nr_tol_delta**2) THEN
               converged = .true.
               EXIT
            ENDIF

            ! 高精度收敛：残差达到极高精度时，完全忽略delta限制
            IF (residual_norm_sq < 1.0e-16_r8) THEN
               converged = .true.
               EXIT
            ENDIF
         
            ! 发散检测 (延迟计算 SQRT)
            IF (newton_iter > 5) THEN
               ! NaN/Inf 检测
               IF (residual_norm /= residual_norm .or. residual_norm > 1.0e10_r8) THEN
                  EXIT
               ENDIF
               
               ! 持续发散检测 (混合版)
               IF (residual_norm > 1.5_r8 * prev_residual_norm .AND. prev_residual_norm > 1.0e-8_r8) THEN
                  IF (residual_norm > 3.0_r8 * prev_residual_norm) THEN
                     EXIT  ! 单次暴涨
                  ELSE
                     residual_div_count = residual_div_count + 1
                     IF (residual_div_count > 3) EXIT  ! 连续温和发散
                  ENDIF
               ELSE
                  residual_div_count = MAX(0, residual_div_count - 1)
               ENDIF
               
               ! 停滞检测 
               IF (ABS(residual_norm - prev_residual_norm) < 1.0e-12_r8 * residual_norm) THEN
                  stagnation_count = stagnation_count + 1
                  IF (stagnation_count > 5 .and. residual_norm < nr_tol_residual * 10.0_r8) THEN
                     converged = .true.
                     EXIT
                  ENDIF
               ELSE
                  stagnation_count = 0
               ENDIF
                  
               prev_residual_norm = residual_norm
            ENDIF
            
         ENDDO  ! End of Newton iterations
         
         IF (converged) THEN
            ANC_final = ANC_iter
            h_final = 10.0_r8 ** StateVector(1)
            h_final = MIN(MAX(h_final, 1.0e-14_r8), 1.0_r8)
            oh_final = MIN(MAX(Kw(j)/h_final, 1.0e-14_r8), 1.0_r8)
            al3_final = MIN(Ksp_gibbsite(j) / (oh_final**3), 1.0_r8)
            ca_final = 10.0_r8 ** StateVector(2)
            ca_final = MAX(ca_final, 1.0e-10_r8)
            DIC_final = 10.0_r8 **StateVector(3)
            DIC_final = MAX(DIC_final, 1.0e-10_r8)
            Bc_final = 10.0_r8 ** StateVector(4)
            Bc_final = MAX(Bc_final, 1.0e-10_r8)

            A_final = h_final**2 + K1(j)*h_final + K1(j)*K2(j)
            h2co3_final = (h_final**2 / A_final) * DIC_final 
            hco3_final = (K1(j)*h_final / A_final) * DIC_final
            co3_final = (K1(j)*K2(j) / A_final) * DIC_final
            
            ca_conc_vr(j,i) = ca_final
            Q_calcite = ca_final * co3_final
            IF (calcite_mass(j,i) > threshold_calcite_mass) THEN
                k_calcite = k0_cal + k_H_cal*h_final + k_OH_cal*oh_final
                rate_cal(j,i) = weight_cal*factor_cal(j)*k_calcite*A_calcite(j) * (1.0_r8 - Q_calcite/Ksp_calcite(j))
                calcite_mass(j,i) = MAX(0.0_r8, calcite_mass(j,i) - rate_cal(j,i) * deltim  * M_CaCO3 * wliq_soisno(j,i))
                calcite_content(j,i) = 100.0_r8 * calcite_mass(j,i) / (BD_all(j,i) * dz_soi(j))
            ELSE
                rate_cal(j,i) = 0.0_r8
            ENDIF

            Bc_conc_vr(j,i) = Bc_final
            Bct_final = ca_final + Bc_final
            base_final = MIN(MAX(SQRT(MAX(Bct_final, 1.0e-20_r8)) / (SQRT(MAX(Bct_final, 1.0e-20_r8)) + K_HBc(j) * h_final), 1.0e-14_r8), 1.0_r8)
            E_Bc(j,i) = base_final * 100_r8
            Bct_conc_vr(j,i) = Bct_final
            
            IF (Bct_final > 1.0e-10_r8 .AND. ca_final > 1.0e-10_r8 .AND. Bc_final > 1.0e-10_r8) THEN
                Bct_ex_final = base_final/2 * CEC_total
                H_ex_final = (1-base_final) * CEC_total
                Bct_ex_flux = weight_cec * factor_cec(j) * (Bct_ex_final - Bct_ex_start) / deltim 
                H_ex_flux = weight_cec * factor_cec(j) * (H_ex_final - H_ex_start) / deltim 
            ELSE
                Bct_ex_flux = 0.0
                H_ex_flux = 0.0
            ENDIF
            rate_Bcex(j,i) = Bct_ex_flux ! K+Ca+Mg
            rate_Hex(j,i) = H_ex_flux
            rate_CO2(j,i) = rate_CO2_explicit
   
            h_conc_vr(j,i) = h_final
            oh_conc_vr(j,i) = oh_final
            al3_conc_vr(j,i) = al3_final
            DIC_conc_vr(j,i) =  DIC_final
            h2co3_conc_vr(j,i) = h2co3_final
            hco3_conc_vr(j,i) = hco3_final
            co3_conc_vr(j,i) = co3_final
            ANC_conc_vr(j,i) = ANC_final
         ELSE
            h_conc_vr(j,i) = h_conc_m1(j,i)
            oh_conc_vr(j,i) = oh_conc_m1(j,i)
            ca_conc_vr(j,i) = ca_conc_m1(j,i)
            DIC_conc_vr(j,i) = DIC_conc_m1(j,i)
            Bc_conc_vr(j,i) = Bc_conc_m1(j,i)
            E_Bc(j,i) = backup_E_Bc
            calcite_mass(j,i) = backup_calcite_mass
            calcite_content(j,i) = backup_calcite_content
            write(*,*) 'Error: Newton-Raphson did not converge at layer', j
         ENDIF
      ENDDO ! nl_soil loop
         
      DO j = 1, nl_soil
         pH_vr(j,i) = -log10(MAX(h_conc_vr(j,i), 1.0e-14_r8))
         ! write(*,*) 'pH', j, pH_vr(j,i)
      ENDDO
      
   END SUBROUTINE soil_acidification_buffering

   ! ----- Subsurface CO2 exchange -----
   SUBROUTINE Gas_liquid_CO2_exchange_M1(i, ps, pe, nl_soil, z_soi, dz_soi, vol_liq, K1, K2, K_H, convert)

   IMPLICIT NONE
   integer,  intent(in)       :: i
   integer,  intent(in)       :: ps                     ! start pft index
   integer,  intent(in)       :: pe                     ! END pft index
   integer,  intent(in)       :: nl_soil
   real(r8), intent(in)       :: z_soi(1:nl_soil)
   real(r8), intent(in)       :: dz_soi(1:nl_soil)
   real(r8), intent(in)       :: vol_liq(1:nl_soil)
   real(r8), intent(in)       :: K1(1:nl_soil), K2(1:nl_soil), K_H(1:nl_soil)
   real(r8), intent(in)       :: convert  ! 默认单位转换常数 (1.957e-3_r8[m2 * atm * g-1 CO2])

   integer :: m, j, ivt, jrt
   real(r8) :: DCO2_s(1:nl_soil), DCO2_s_W(1:nl_soil), DCO2_s_C(1:nl_soil), DCO2_s_gas_C(1:nl_soil), DCO2_s_water_C(1:nl_soil), h2co3_soil(1:nl_soil)
   real(r8) :: h2co3_start, pCO2_atm, pCO2_max, pCO2_min, ss_vliq_min = 5e-2_r8, log10_pCO2
   real(r8) :: roota, factor = 0.5_r8
   real(r8) :: f_a, eps, om_frac, theta, sum_D, sum_dz,D_avg
   real(r8) :: rgas = 8.314467591                 ! universal gas constant (J mol-1 K-1)
      
      pCO2_atm = forc_pco2m(i) / 1.013e5_r8  ! atm

      DO m = ps, pe
         ivt = pftclass(m)
         roota = 0.0_r8
         jrt = 1
         IF (ivt /= 0) THEN
            DO j = 1, nl_soil
               roota = roota + rootfr_p(j, ivt)
               IF (roota >= 0.9_r8) THEN 
                  jrt = j
                  EXIT
               ENDIF
            ENDDO

            !!! calculate soil gas diffusivity (m2/s)
            DO j = 1, nl_soil
               ! method 1: WITCH model
               DCO2_s_W(j) = 0.139e-4_r8 * ((t_soisno(j,i)/273.15)**2) * porsl(j,i) * 0.45_r8

               ! method 2: CLM5 (Wania et al., 2010)
               f_a = 1._r8 - wfc(j,i) / porsl(j,i)
               eps =  porsl(j,i)-wfc(j,i) ! Air-filled fraction of total soil volume

               ! use diffusivity calculation including peat
               IF (organic_max > 0._r8) THEN
                  om_frac = min(OM_density(j,i)/organic_max, 1._r8)
                  ! Use first power, not square as in iniTimeConst
               ELSE
                  om_frac = 1._r8
               ENDIF
               DCO2_s_gas_C(j) = (d_con_g21_CO2 + d_con_g22_CO2*t_soisno(j,i)) * 1.e-4_r8 * &
                  (om_frac * f_a**(10._r8/3._r8) / porsl(j,i)**2 + &
                  (1._r8-om_frac) * eps**2 * f_a**(3._r8 / bsw(j,i)) ) 
               
               DCO2_s_water_C(j) =  ((d_con_w21_CO2 + d_con_w22_CO2*t_soisno(j,i) + d_con_w23_CO2*t_soisno(j,i)**2) * 1.e-9_r8)*porsl(j,i)**2
            ENDDO
      
            ! calculate weighted average diffusion coefficient
            sum_D  = SUM(DCO2_s_gas_C(1:jrt) * dz_soi(1:jrt))
            sum_dz = SUM(dz_soi(1:jrt))
            D_avg  = sum_D / sum_dz

            !!! calculate soil CO2 pressure (ppmv)
            ! WITCH model  默认单位转换常数 (1.957e-3_r8[m2 * atm * g-1 CO2])
            DO j = 1, nl_soil
               IF (soil_er(i) > 1e-12_r8) THEN
                  pCO2_hr(j,i) = convert * (soil_er(i)*(44/12)*(z_soi(j))**2) / (2.0_r8*porsl(j,i)*D_avg)*1e6_r8 !
               ELSE
                  pCO2_hr(j,i) = 0.0_r8
               ENDIF
            ENDDO

            pCO2_max = MIN(pCO2_atm*1e6_r8 + pCO2_hr(jrt,i), 1e5_r8)  
            pCO2_min = pCO2_atm*1e6_r8
            DO j = 1, nl_soil
               pCO2_soil(j,i) = pCO2_atm*1e6_r8
               IF (soil_er(i) > 1e-12_r8 .and. pCO2_max > pCO2_atm*1e6_r8) THEN
                  IF (j <= jrt) THEN
                     pCO2_soil(j,i) = MAX(pCO2_min, pCO2_min + (pCO2_max - pCO2_min) *(z_soi(j) / z_soi(jrt))**1.5_r8) ! 幂函数衰减
                  ELSE
                     pCO2_soil(j,i) = pCO2_max
                  ENDIF
               ENDIF
               ! write(*,*) 'M1 ', j, ': pCO2_soil =', pCO2_soil(j,i), 'ppmv'
            ENDDO
         ENDIF
      ENDDO

   END SUBROUTINE Gas_liquid_CO2_exchange_M1

   SUBROUTINE Gas_liquid_CO2_exchange_M2(i, ps, pe, nl_soil, z_soi, dz_soi, vol_liq, K1, K2, K_H, convert)

   IMPLICIT NONE
   integer,  intent(in)       :: i
   integer,  intent(in)       :: ps                     ! start pft index
   integer,  intent(in)       :: pe                     ! END pft index
   integer,  intent(in)       :: nl_soil
   real(r8), intent(in)       :: z_soi(1:nl_soil)
   real(r8), intent(in)       :: dz_soi(1:nl_soil)
   real(r8), intent(in)       :: vol_liq(1:nl_soil)
   real(r8), intent(in)       :: K1(1:nl_soil), K2(1:nl_soil), K_H(1:nl_soil)
   real(r8), intent(in)       :: convert  ! 默认单位转换常数 (1.957e-3_r8[m2 * atm * g-1 CO2])

   !--- local variables ---
   integer  :: j
   real(r8) :: f_a, om_frac
   real(r8) :: DCO2_s_gas_C(1:nl_soil)   ! effective gas-phase diffusivity [m²/s]
   real(r8) :: DCO2_s_water_C(1:nl_soil) ! effective liquid-phase diffusivity [m²/s]
   real(r8) :: eps(1:nl_soil)          ! air-filled porosity per layer [-]: porsl - wfc
   real(r8) :: pCO2_atm                  ! atmospheric pCO2 [ppmv]
   real(r8) :: C_atm                     ! atmospheric CO2 concentration [mol/m³ air]
   real(r8) :: RT_P                      ! R*T/P_atm [m³/mol]
   real(r8) :: rgas = 8.314467591        ! universal gas constant (J mol-1 K-1)

   !--- tridiagonal system: -aa(j)*C(j-1) + bb(j)*C(j) - cc(j)*C(j+1) = ff(j) ---
   ! aa(j): conductance at top interface of layer j    [m/s]
   ! cc(j): conductance at bottom interface of layer j [m/s]
   ! bb(j) = aa(j) + cc(j)                            [m/s]
   ! ff(j): source term = S_j * dz_j                  [mol CO2/m²/s]
   ! C_sol(j): solution, CO2 concentration             [mol CO2/m³ air]
   real(r8) :: aa(1:nl_soil), bb(1:nl_soil), cc(1:nl_soil), ff(1:nl_soil)
   real(r8) :: C_sol(1:nl_soil)
   real(r8) :: cp(1:nl_soil), dp(1:nl_soil)  ! Thomas algorithm work arrays
   real(r8) :: exchange_calib             ! [s·m²/gC] 校准参数
   real(r8) :: source_raw(1:nl_soil)      ! 原始源项
   real(r8) :: D_times_dz(1:nl_soil)      ! D * dz (扩散能力指标)
   real(r8) :: exchange_limit(1:nl_soil)   ! 交换限制因子

      exchange_calib = DEF_SOILACID_exchange_calib

      pCO2_atm = forc_pco2m(i) / 1.013e5_r8 * 1.0e6_r8  ! [Pa] → [ppmv]

      ! Step A: Effective CO2 diffusivity per layer [m²/s]  (CLM5 method)
      DO j = 1, nl_soil
         f_a     = 1._r8 - wfc(j,i) / porsl(j,i)
         eps(j)     = porsl(j,i) - wfc(j,i)
         IF (organic_max > 0._r8) THEN
            om_frac = MIN(OM_density(j,i)/organic_max, 1._r8)
         ELSE
            om_frac = 1._r8
         ENDIF
         DCO2_s_gas_C(j) = (d_con_g21_CO2 + d_con_g22_CO2*t_soisno(j,i)) * 1.e-4_r8 * &
            ( om_frac * f_a**(10._r8/3._r8) / porsl(j,i)**2 + &
             (1._r8-om_frac) * eps(j)**2 * f_a**(3._r8/bsw(j,i)) )
         DCO2_s_water_C(j) = ((d_con_w21_CO2 + d_con_w22_CO2*t_soisno(j,i) + &
            d_con_w23_CO2*t_soisno(j,i)**2) * 1.e-9_r8) * porsl(j,i)**2
         ! write(*,*) 'M2 ', j, ': DCO2_s_gas_C =', DCO2_s_gas_C(j), ' DCO2_s_water_C =', DCO2_s_water_C(j), DCO2_s_water_C(j) / DCO2_s_gas_C(j)
      ENDDO
         

      ! Step B: Assemble tridiagonal system
      DO j = 1, nl_soil
         ! Source [mol CO2/m²/s]: soil_er_vr [gC/m³/s] / 12.011 [g/mol] × dz [m]
         ! NOTE: no eps division - DCO2_s_gas_C is bulk diffusivity, C_sol will be [mol/m³ air]
         source_raw(j) = MAX(soil_er_vr(j,i), 0.0_r8) / 12.011_r8 * dz_soi(j)

         ! Calculate diffusion capacity indicator [m³/s]
         D_times_dz(j) = DCO2_s_gas_C(j) * dz_soi(j)

         ! Exchange limit: limit source based on diffusion capacity
         ! exchange_calib越小 → 限制越宽松（limit越大）
         exchange_limit(j) = 1.0_r8 / (1.0_r8 + (source_raw(j) * exchange_calib) / D_times_dz(j))

         ! Apply exchange rate limit
         ff(j) = source_raw(j) * exchange_limit(j)

         ! write(*,*) 'M2 ', j, ': er_vr=', soil_er_vr(j,i), ' source_raw=', source_raw(j), ' limit=', exchange_limit(j), ' ff=', ff(j)

         ! Conductance at top interface [m/s]
         IF (j == 1) THEN
            aa(j) = DCO2_s_gas_C(1) / (0.5_r8 * dz_soi(1))
         ELSE
            aa(j) = 1.0_r8 / ( 0.5_r8*dz_soi(j-1)/DCO2_s_gas_C(j-1) &
                              + 0.5_r8*dz_soi(j)  /DCO2_s_gas_C(j) )
         ENDIF

         ! Conductance at bottom interface [m/s]
         IF (j == nl_soil) THEN
            cc(j) = 0.0_r8   ! no-flux lower boundary
         ELSE
            cc(j) = 1.0_r8 / ( 0.5_r8*dz_soi(j)  /DCO2_s_gas_C(j) &
                              + 0.5_r8*dz_soi(j+1)/DCO2_s_gas_C(j+1) )
         ENDIF

         bb(j) = aa(j) + cc(j)
      ENDDO

      ! Incorporate top Dirichlet BC: C(0) = C_atm [mol/m³ air]
      ! C_atm = n/V from ideal gas: pCO2[atm] / (RT/P) = pCO2_frac / (RT/P)
      RT_P  = rgas * MAX(t_soisno(1,i), 250.0_r8) / 1.013e5_r8  ! [m³/mol]
      C_atm = pCO2_atm * 1.0e-6_r8 / RT_P                        ! [ppmv→frac] / [m³/mol] = [mol/m³ air]
      ff(1) = ff(1) + aa(1) * C_atm

      ! Step C: Thomas algorithm (TDMA) — O(N) exact direct solver
      ! Forward sweep
      cp(1) = cc(1) / bb(1)
      dp(1) = ff(1) / bb(1)
      DO j = 2, nl_soil
         cp(j) = cc(j) / ( bb(j) - aa(j)*cp(j-1) )
         dp(j) = ( ff(j) + aa(j)*dp(j-1) ) / ( bb(j) - aa(j)*cp(j-1) )
      ENDDO

      ! Back substitution
      C_sol(nl_soil) = dp(nl_soil)
      DO j = nl_soil - 1, 1, -1
         C_sol(j) = dp(j) + cp(j) * C_sol(j+1)
      ENDDO

      !====================================================================
      ! Step D: Convert C [mol/m³ air] → pCO2 [ppmv]  using pV=nRT
      ! C_sol is in [mol/m³ air] because ff was divided by eps
      ! pCO2 [ppmv] = C [mol/m³] × RT/P [m³/mol] × 1e6
      !====================================================================
      DO j = 1, nl_soil
         RT_P = rgas * MAX(t_soisno(j,i), 250.0_r8) / 1.013e5_r8  ! [m³/mol]
         pCO2_soil(j,i) = C_sol(j) * RT_P * 1.0e6_r8
         pCO2_soil(j,i) = MAX(pCO2_soil(j,i), pCO2_atm)     ! floor at atmospheric
         pCO2_soil(j,i) = MIN(pCO2_soil(j,i), 1.0e5_r8)     ! cap at 100 000 ppmv
         pCO2_hr(j,i)   = pCO2_soil(j,i) - pCO2_atm         ! diagnostic: excess [ppmv]
         ! write(*,*) 'M2 ', j, ': pCO2_soil =', pCO2_soil(j,i), 'ppmv', ' soil_er_vr =', soil_er_vr(j,i), 'gC/m3/s'
      ENDDO

   END SUBROUTINE Gas_liquid_CO2_exchange_M2

   SUBROUTINE calc_equilibrium_constants(T_K, K_H, K1, K2, Ksp_calcite, Kw, &
                                          Kal_1, Kal_2, Kal_3, Ksp_gibbsite)

   IMPLICIT NONE

   real(r8), intent(in)  :: T_K
   real(r8), intent(out) :: K_H, K1, K2, Ksp_calcite, Kw
   real(r8), intent(out) :: Kal_1, Kal_2, Kal_3, Ksp_gibbsite

   real(r8) :: log_KH, log_K1, log_K2, log_Ksp_calcite, log_Kw
   real(r8) :: log_Kal_1, log_Kal_2, log_Kal_3

      ! Plummer, L. N., & Busenberg, E. (1982); Stumm and Morgan (1996)
      log_KH = 108.3865_r8 + 0.01985076_r8 * T_K - 6919.53_r8 / T_K &
      - 40.45154_r8 * log10(T_K) + 669365_r8 / (T_K**2)
      K_H = 10.0_r8 ** log_KH  ! Henry constant (mol/L/atm)  CO2(g) = CO2(aq)

      log_K1 = -356.3094_r8 - 0.06091964_r8 * T_K + 21834.37_r8 / T_K &
      + 126.8339_r8 * log10(T_K) - 1684915_r8 / (T_K**2)
      K1 = 10.0_r8 ** log_K1  ! CO2(aq) + H2O = H+ + HCO3-

      log_K2 = -107.8871_r8 - 0.03252849_r8 * T_K + 5151.79_r8 / T_K &
      + 38.92561_r8 * log10(T_K) - 563713.9_r8 / (T_K**2)
      K2 = 10.0_r8 ** log_K2  ! HCO3- = H+ + CO3 2-

      log_Ksp_calcite = -171.0965_r8 - 0.077993_r8 * T_K + 2839.319_r8 / T_K &
      + 71.595_r8 * log10(T_K)
      Ksp_calcite = 10.0_r8 ** log_Ksp_calcite  ! CaCO3: Ksp_calcite = [Ca2+][CO3-]
      ! Ksp_calcite = 4.8e-9_r8  ! General Chemistry，Pauling，1970

      ! Stumm and Morgan (1996)
      ! log_Kw = -4470.99_r8/T_K + 6.0875_r8 - 0.01706_r8 * T_K
      log_Kw = -283.9710_r8 + 13323.0_r8/T_K - 0.05069842_r8 * T_K &
      + 102.24447_r8 * log10(T_K) - 1119669_r8 / (T_K**2)
      Kw = 10.0_r8 ** log_Kw

      ! Al3+  Stumm and Morgan (1996)
      log_Kal_1 = -38.253_r8 - 656.27_r8/T_K + 14.327_r8 * log10(T_K) ! Al3+ + H2O = AlOH2+ + H+
      Kal_1 = 10.0_r8 ** log_Kal_1
      log_Kal_2 = 88.5_r8 - 9391.6_r8 / T_K - 27.121_r8 * log10(T_K)  ! Al3+ + 2H2O = Al(OH)2+ + 2H+
      Kal_2 = 10.0_r8 ** log_Kal_2
      log_Kal_3 = 226.374_r8 - 18247.8_r8/T_K - 73.597_r8 * log10(T_K) ! Al3+ + 3H2O = Al(OH)3 + 3H+
      Kal_3 = 10.0_r8 ** log_Kal_3

      Ksp_gibbsite = 1.3e-33_r8  ! Al(OH)3 = Al3+ + 3OH-
      ! Ksp_gibbsite = Kw**3 / Kal_3

   END SUBROUTINE calc_equilibrium_constants

   ! 修正的LU分解求解线性方程组 - 适用于酸化-缓冲体系
   SUBROUTINE SolveLinearSystem(n, Jacobian, Residuals, DeltaX, success)
   IMPLICIT NONE

   integer, intent(in) :: n
   real(r8), intent(in) :: Jacobian(n,n)
   real(r8), intent(in) :: Residuals(n)
   real(r8), intent(out) :: DeltaX(n)
   logical, intent(out) :: success  ! 返回求解状态

   real(r8) :: LU(n,n), RHS(n)
   real(r8) :: temp, row_scale(n), col_scale(n)
   real(r8) :: local_max
   real(r8), parameter :: eps = 1.0e-14_r8
   real(r8), parameter :: sing_tol = 1.0e-12_r8
   integer :: m, k, l, pivot_row

      success = .true.
      LU = Jacobian
      RHS = -Residuals

      ! 行列平衡缩放 (改进数值稳定性)
      DO m = 1, n
         row_scale(m) = MAXVAL(ABS(LU(m,:)))
         IF (row_scale(m) < eps) row_scale(m) = 1.0_r8
         col_scale(m) = MAXVAL(ABS(LU(:,m)))
         IF (col_scale(m) < eps) col_scale(m) = 1.0_r8
      ENDDO

      ! 应用平衡缩放
      DO m = 1, n
         LU(m,:) = LU(m,:) / row_scale(m)
         LU(:,m) = LU(:,m) / col_scale(m)
         RHS(m) = RHS(m) / row_scale(m)
      ENDDO

      ! LU分解 (带部分主元) 
      DO m = 1, n-1
         ! 寻找列主元
         pivot_row = m
         DO k = m+1, n
            IF (ABS(LU(k,m)) > ABS(LU(pivot_row,m))) THEN
               pivot_row = k
            ENDIF
         ENDDO
         
         ! 行交换 (同时交换缩放因子)
         IF (pivot_row /= m) THEN
            ! 交换 RHS
            temp = RHS(m)
            RHS(m) = RHS(pivot_row)
            RHS(pivot_row) = temp
            
            ! 交换 LU 行
            DO l = 1, n
               temp = LU(m,l)
               LU(m,l) = LU(pivot_row,l)
               LU(pivot_row,l) = temp
            ENDDO
            
            ! 交换缩放因子
            temp = row_scale(m)
            row_scale(m) = row_scale(pivot_row)
            row_scale(pivot_row) = temp
            
         ENDIF

         ! 奇异性检查
         local_max = MAXVAL(ABS(LU(m:n,m:n)))
         IF (ABS(LU(m,m)) < sing_tol * local_max) THEN
            write(*,*) 'Warning: Near-singular pivot at row ',m, ', value = ', LU(m,m)
            success = .false.
            DeltaX = 0.0_r8
            RETURN
         ENDIF

         ! 高斯消元
         DO l = m+1, n
            IF (ABS(LU(l,m)) > eps) THEN
               temp = LU(l,m) / LU(m,m)
               LU(l,m) = temp
               LU(l,m+1:n) = LU(l,m+1:n) - temp * LU(m,m+1:n)
               RHS(l) = RHS(l) - temp * RHS(m)
            ELSE
               LU(l,m) = 0.0_r8
            ENDIF
         ENDDO
      ENDDO

      ! 检查最后一个对角元
      local_max = ABS(LU(n,n))
      IF (ABS(LU(n,n)) < sing_tol * local_max) THEN
         write(*,*) 'Warning: Singular last pivot, value = ', LU(n,n)
         success = .false.
         DeltaX = 0.0_r8
         RETURN
      ENDIF

      ! 回代求解
      DO m = n, 1, -1
         IF (ABS(LU(m,m)) < eps) THEN
            write(*,*) 'Error: Zero diagonal in back-substitution at row', m
            success = .false.
            DeltaX = 0.0_r8
            RETURN
         ENDIF
         
         ! 向量化点积
         temp = RHS(m)
         DO l = m+1, n
            temp = temp - LU(m,l) * RHS(l)
         ENDDO
         RHS(m) = temp / LU(m,m)
      ENDDO

      ! 恢复原始缩放
      DO m = 1, n
         ! 只除以列缩放,行缩放已经在第一步应用
         DeltaX(m) = RHS(m) / col_scale(m)
      ENDDO

   END SUBROUTINE SolveLinearSystem

   ! 自适应牛顿阻尼 - 带有回溯线搜索（Armijo条件）
   SUBROUTINE AdaptiveDamping(StateVector, DeltaX, current_residuals, current_norm, tolerance, &
                            step_scale, j, i, dt, K1_val, K2_val, Kw_val, Ksp_calcite_val, &
                            Kal_1_val, Kal_2_val, Kal_3_val, Ksp_gibbsite_val, CEC_total, &
                            weight_cal, weight_cec, factor_cal, factor_cec, &
                            k0_cal, k_H_cal, k_OH_cal, A_calcite, K_HBc, factor_co2, &
                            h2co3_soil_val, ANC_start, ca_start, DIC_start, Bc_start, &
                            H_ex_start, Bct_ex_start)
   IMPLICIT NONE
   
   real(r8), intent(in) :: StateVector(4), DeltaX(4), current_residuals(4)
   real(r8), intent(in) :: current_norm, tolerance
   real(r8), intent(inout) :: step_scale
   integer, intent(in) :: j, i
   real(r8), intent(in) :: dt, K1_val, K2_val, Kw_val, Ksp_calcite_val
   real(r8), intent(in) :: Kal_1_val, Kal_2_val, Kal_3_val, Ksp_gibbsite_val, CEC_total
   real(r8), intent(in) :: weight_cal, weight_cec, factor_cal, factor_cec
   real(r8), intent(in) :: k0_cal, k_H_cal, k_OH_cal, A_calcite, K_HBc, factor_co2
   real(r8), intent(in) :: h2co3_soil_val, ANC_start, ca_start, DIC_start, Bc_start
   real(r8), intent(in) :: H_ex_start, Bct_ex_start
   
   real(r8) :: test_state(4), test_residuals(4), test_norm
   real(r8) :: alpha, min_alpha, reduction_factor, sufficient_decrease
   integer  :: max_backtrack, backtrack_iter, k
   logical  :: acceptable_step, bounds_violated
   
      ! 算法参数
      alpha = 1.0_r8
      min_alpha = 1.0e-6_r8
      reduction_factor = 0.5_r8
      sufficient_decrease = 1.0e-4_r8  ! Armijo 条件
      max_backtrack = 15
      acceptable_step = .false.

      ! Armijo线搜索：寻找足够的残差减少
      DO backtrack_iter = 1, max_backtrack
         test_state = StateVector + alpha * DeltaX
         bounds_violated = .false.
         
         ! 检查边界约束：不同变量有不同约束
         IF (test_state(1) < -14.0_r8 .or. test_state(1) > 0.0_r8) THEN  ! logH
            bounds_violated = .true.
         ENDIF
         ! Ca (index 2) 和 Bc (index 4) 的下限为 1.0e-10
         IF (test_state(2) < -10.0_r8 .or. test_state(2) > 2.0_r8) THEN  ! logCa
            bounds_violated = .true.
         ENDIF
         IF (test_state(3) < -10.0_r8 .or. test_state(3) > 2.0_r8) THEN  ! logDIC
            bounds_violated = .true.
         ENDIF
         IF (test_state(4) < -10.0_r8 .or. test_state(4) > 2.0_r8) THEN  ! logBc
            bounds_violated = .true.
         ENDIF
         ! 如果违反边界,缩小步长重试
         IF (bounds_violated) THEN
            alpha = alpha * reduction_factor
            IF (alpha < min_alpha) EXIT
            CYCLE
         ENDIF

         ! 计算新状态下的残差
         CALL EvaluateResiduals(test_state, test_residuals, j, i, dt, &
                               K1_val, K2_val, Kw_val, Ksp_calcite_val, Kal_1_val, Kal_2_val, Kal_3_val, &
                               Ksp_gibbsite_val, CEC_total, weight_cal, weight_cec, &
                               factor_cal, factor_cec, k0_cal, k_H_cal, k_OH_cal, &
                               A_calcite, K_HBc, factor_co2, h2co3_soil_val, &
                               ANC_start, ca_start, DIC_start, Bc_start, H_ex_start, Bct_ex_start)
         
         test_norm = SQRT(MAX(SUM(test_residuals**2), 1.0e-20_r8))

         ! 准则 1: Armijo充分下降条件
         ! ||F(x + α·Δx)||² ≤ ||F(x)||² - c·α·||F(x)||²
         IF (test_norm**2 <= current_norm**2 * (1.0_r8 - sufficient_decrease * alpha)) THEN
            acceptable_step = .true.
            EXIT
         ENDIF
         
         ! 准则 2: 简单下降（当Armijo过于严格时的备选）
         IF (test_norm < current_norm) THEN
            acceptable_step = .true.
            EXIT
         ENDIF
               
         ! 准则3: 近收敛容忍（已接近解）
         IF (test_norm < tolerance * 5.0_r8 .and. &
            test_norm < current_norm * 1.02_r8) THEN
            acceptable_step = .true.
            EXIT
         ENDIF
         
         ! 准则4: 微小增长容忍（数值噪声范围内）- 更严格
         IF (test_norm < current_norm * 1.0001_r8 .and. current_norm < tolerance * 50.0_r8) THEN
            acceptable_step = .true.
            EXIT
         ENDIF
         
         alpha = alpha * reduction_factor
      
         ! 防止步长过小
         IF (alpha < min_alpha) EXIT
      ENDDO
         
      ! 如果找不到合适步长，使用保守的最小步长
      IF (.NOT. acceptable_step) THEN
         alpha = min_alpha
         IF (current_norm > tolerance * 100.0_r8) THEN
         ENDIF
      ENDIF
      step_scale = alpha
      
   END SUBROUTINE AdaptiveDamping

   ! 计算给定状态下的残差
   SUBROUTINE EvaluateResiduals(StateVector, Residuals, j, i, dt, &
                               K1_val, K2_val, Kw_val, Ksp_calcite_val, Kal_1_val, Kal_2_val, Kal_3_val, &
                               Ksp_gibbsite_val, CEC_total, weight_cal, weight_cec, &
                               factor_cal, factor_cec, k0_cal, k_H_cal, k_OH_cal, &
                               A_calcite, K_HBc, factor_co2, h2co3_soil_val, &
                               ANC_start, ca_start, DIC_start, Bc_start, H_ex_start, Bct_ex_start)
   IMPLICIT NONE
   
   real(r8), intent(in) :: StateVector(4)
   real(r8), intent(out) :: Residuals(4)
   integer, intent(in) :: j, i
   real(r8), intent(in) :: dt, K1_val, K2_val, Kw_val, Ksp_calcite_val
   real(r8), intent(in) :: Kal_1_val, Kal_2_val, Kal_3_val, Ksp_gibbsite_val, CEC_total
   real(r8), intent(in) :: weight_cal, weight_cec, factor_cal, factor_cec
   real(r8), intent(in) :: k0_cal, k_H_cal, k_OH_cal, A_calcite, K_HBc, factor_co2
   real(r8), intent(in) :: h2co3_soil_val, ANC_start, ca_start, DIC_start, Bc_start
   real(r8), intent(in) :: H_ex_start, Bct_ex_start
   
   ! 局部变量与主牛顿循环相同
   real(r8) :: logH_iter, logca_iter, logDIC_iter, logBc_iter
   real(r8) :: h_iter, oh_iter, DIC_iter, A_iter
   real(r8) :: h2co3_iter, hco3_iter, co3_iter, ANC_iter
   real(r8) :: al3_iter, al2_iter, al_iter, al_alk_iter
   real(r8) :: ca_iter, Bc_iter, Bct_iter, base_iter
   real(r8) :: rate_CaCO3, Q_calcite, saturation_ratio, k_calcite, CO2_trans
   real(r8) :: Bct_ex_iter, H_ex_iter, H_ex_flux, Bct_ex_flux, ca_ex_flux, Bc_ex_flux
   real(r8) :: pH_iter, al_weight, k_sigmoid = 5.0_r8
   real(r8) :: threshold_calcite_mass = 1.0e-14_r8

      logH_iter = MIN(MAX(StateVector(1), -14.0_r8), 0.0_r8)
      h_iter = MIN(MAX(10.0_r8 ** logH_iter, 1.0e-14_r8), 1.0_r8)
      oh_iter = MIN(MAX(Kw_val / h_iter, 1.0e-14_r8), 1.0_r8)

      logca_iter = MIN(MAX(StateVector(2), -10.0_r8), 2.0_r8)
      ca_iter = MIN(MAX(10.0_r8 ** logca_iter, 1.0e-10_r8), 100.0_r8)
      logDIC_iter = MIN(MAX(StateVector(3), -10.0_r8), 2.0_r8)
      DIC_iter = MIN(MAX(10.0_r8 ** logDIC_iter, 1.0e-10_r8), 100.0_r8)
      logBc_iter = MIN(MAX(StateVector(4), -10.0_r8), 2.0_r8)
      Bc_iter = MIN(MAX(10.0_r8 ** logBc_iter, 1.0e-10_r8), 100.0_r8)

      A_iter = h_iter**2 + K1_val*h_iter + K1_val*K2_val
      h2co3_iter = (h_iter**2/A_iter) * DIC_iter
      hco3_iter = (K1_val*h_iter / A_iter) * DIC_iter
      co3_iter = (K1_val*K2_val / A_iter) * DIC_iter

      al3_iter = MIN(Ksp_gibbsite_val / (oh_iter**3), 1.0_r8)
      al2_iter = al3_iter * Kal_2_val/(h_iter**2)
      al_iter = al3_iter * Kal_1_val/h_iter
      al_alk_iter = 3.0_r8*al3_iter + 2.0_r8*al_iter + al2_iter

      ! IF (-LOG10(h_iter) > 5.5_r8) THEN
      !    al_alk_iter = 0.0_r8
      ! ENDIF

      !=======================================================================
      pH_iter     = -LOG10(h_iter)
      al_weight   = 1.0_r8 / (1.0_r8 + EXP(k_sigmoid * (pH_iter - 5.5_r8)))
      al_alk_iter = al_alk_iter * al_weight
      !=======================================================================
      

      ANC_iter = hco3_iter + 2.0_r8*co3_iter + oh_iter - h_iter - al_alk_iter

      IF (calcite_mass(j,i) > threshold_calcite_mass) THEN
          Q_calcite = ca_iter * co3_iter
          saturation_ratio = Q_calcite / Ksp_calcite_val
          k_calcite = k0_cal + k_H_cal * h_iter + k_OH_cal * oh_iter
          rate_CaCO3 = weight_cal * factor_cal * k_calcite * A_calcite * (1.0_r8 - saturation_ratio)
      ELSE
          rate_CaCO3 = 0.0_r8
      ENDIF

      Bct_iter = ca_iter + Bc_iter
      IF (Bct_iter > 1.0e-10_r8 .AND. ca_iter > 1.0e-10_r8 .AND. Bc_iter > 1.0e-10_r8) THEN
         base_iter = SQRT(MAX(Bct_iter, 1.0e-20_r8)) / (SQRT(MAX(Bct_iter, 1.0e-20_r8)) + K_HBc * h_iter)
         Bct_ex_iter = base_iter/2 * CEC_total
         H_ex_iter = (1-base_iter) * CEC_total
         H_ex_flux = weight_cec * factor_cec * (H_ex_iter - H_ex_start) / dt
         Bct_ex_flux = weight_cec * factor_cec * (Bct_ex_iter - Bct_ex_start) / dt
         ca_ex_flux = Bct_ex_flux * (ca_iter / Bct_iter)
         Bc_ex_flux = Bct_ex_flux * (Bc_iter / Bct_iter)
      ELSE
         H_ex_flux = 0.0_r8
         Bct_ex_flux = 0.0_r8
         ca_ex_flux = 0.0_r8
         Bc_ex_flux = 0.0_r8
      ENDIF

      Residuals(1) = (ANC_iter - ANC_start)/dt - 2.0_r8*rate_CaCO3 - H_ex_flux + h_from_N(j,i)
      Residuals(2) = (ca_iter - ca_start)/dt - rate_CaCO3 + ca_ex_flux          
      Residuals(3) = (DIC_iter - DIC_start)/dt - rate_CaCO3
      Residuals(4) = (Bc_iter - Bc_start)/dt + Bc_ex_flux
      
   END SUBROUTINE EvaluateResiduals

END MODULE MOD_BGC_Soil_BiogeochemNStateUpdate1
#endif
