#include <define.h>

#ifdef CROP
MODULE MOD_CropReadin

!-----------------------------------------------------------------------
   USE MOD_Precision
   IMPLICIT NONE
   SAVE

   ! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: CROP_readin

CONTAINS

   SUBROUTINE CROP_readin ()
!-----------------------------------------------------------------------
! !DESCRIPTION:
!  Read in crop planting date from data, and fertilization from data.
!  Save these data in patch vector.
!
!  Original: Shupeng Zhang, Zhongwang Wei, and Xingjie Lu, 2022
!
!-----------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_SPMD_Task
   USE MOD_LandPatch
   USE MOD_NetCDFSerial
   USE MOD_NetCDFBlock
   USE MOD_SpatialMapping
   USE MOD_Vars_TimeInvariants
   USE MOD_Vars_TimeVariables
   USE MOD_BGC_Vars_TimeVariables

   USE MOD_Vars_Global
   USE MOD_LandPFT
   USE MOD_Vars_PFTimeVariables
   USE MOD_RangeCheck
   USE MOD_Block

   IMPLICIT NONE

   character(len=256) :: file_crop
   type(grid_type)    :: grid_crop
   type(block_data_real8_2d)  :: f_xy_crop
   type(spatial_mapping_type) :: mg2patch_crop
   type(spatial_mapping_type) :: mg2pft_crop
   character(len=256) :: file_irrig
   type(grid_type)    :: grid_irrig
   type(block_data_int32_2d)  :: f_xy_irrig
   type(spatial_mapping_type) :: mg2pft_irrig

   character(len=256) :: file_irrigalloc
   type(grid_type)    :: grid_irrigalloc
   type(block_data_real8_2d)  :: f_xy_irrigalloc
   type(spatial_mapping_type) :: mg2p_irrigalloc

   character(len=256) :: file_fert
   type(grid_type)    :: grid_fert
   type(block_data_real8_2d) :: f_xy_fert
   type(spatial_mapping_type) :: mg2pft_fert

   character(LEN=256) :: file_acid
   type(grid_type)    :: grid_acid
   type(block_data_real8_2d) :: f_xy_acid
   type(spatial_mapping_type) :: mg2pft_acid
   character(len=256) :: file_CEC
   type(grid_type)    :: grid_CEC
   type(block_data_real8_2d)  :: f_xy_CEC
   type(spatial_mapping_type) :: mg2pft_CEC
   character(len=256) :: file_CACO3
   type(grid_type)    :: grid_CACO3
   type(block_data_real8_2d)  :: f_xy_CACO3
   type(spatial_mapping_type) :: mg2pft_CACO3
   character(len=256) :: file_BS
   type(grid_type)    :: grid_BS
   type(block_data_real8_2d)  :: f_xy_BS
   type(spatial_mapping_type) :: mg2pft_BS

   real(r8), allocatable :: pH_vr_tmp(:)
   real(r8), allocatable :: calcite_content_tmp(:)
   real(r8), allocatable :: CEC_tmp(:)
   real(r8), allocatable :: E_Bc_tmp(:)

   real(r8),allocatable :: pdrice2_tmp      (:)
   real(r8),allocatable :: plantdate_tmp    (:)
   real(r8),allocatable :: fertnitro_tmp    (:)
   integer ,allocatable :: irrig_method_tmp (:)

   real(r8),allocatable :: fertilizer_tmp   (:)
   real(r8),allocatable :: manure_tmp       (:)

   ! Local variables
   real(r8), allocatable :: lat(:), lon(:)
   real(r8) :: missing_value
   integer(kind=2) :: missing_value_i2
   integer  :: cft, npatch, ipft, nsl
   character(len=2) :: cx
   integer  :: iblkme, iblk, jblk
   integer  :: maxvalue, minvalue

      ! READ in crops

      file_crop = trim(DEF_dir_runtime) // '/crop/plantdt-colm-64cfts-rice2_fillcoast.nc'

      CALL ncio_read_bcast_serial (file_crop, 'lat', lat)
      CALL ncio_read_bcast_serial (file_crop, 'lon', lon)

      CALL grid_crop%define_by_center (lat, lon)

      IF (p_is_io) THEN
         CALL allocate_block_data  (grid_crop, f_xy_crop)
      ENDIF

      ! missing value
      IF (p_is_master) THEN
         CALL ncio_get_attr (file_crop, 'pdrice2', 'missing_value', missing_value)
      ENDIF
#ifdef USEMPI
      CALL mpi_bcast (missing_value, 1, MPI_REAL8, p_address_master, p_comm_glb, p_err)
#endif
      IF (p_is_io) THEN
         CALL ncio_read_block (file_crop, 'pdrice2', grid_crop, f_xy_crop)
      ENDIF

      CALL mg2patch_crop%build_arealweighted (grid_crop, landpatch)
      CALL mg2patch_crop%set_missing_value   (f_xy_crop, missing_value)

      CALL mg2pft_crop%build_arealweighted   (grid_crop, landpft)
      CALL mg2pft_crop%set_missing_value     (f_xy_crop, missing_value)

      IF (allocated(lon)) deallocate(lon)
      IF (allocated(lat)) deallocate(lat)

      IF (p_is_worker) THEN
         IF (numpatch > 0)  allocate(pdrice2_tmp    (numpatch))
         IF (numpft   > 0)  allocate(plantdate_tmp    (numpft))
         IF (numpft   > 0)  allocate(fertnitro_tmp    (numpft))
         IF (numpft   > 0)  allocate(irrig_method_tmp (numpft))
      ENDIF

      ! (1) Read in plant date for rice2.
      file_crop = trim(DEF_dir_runtime) // '/crop/plantdt-colm-64cfts-rice2_fillcoast.nc'
      IF (p_is_io) THEN
         CALL ncio_read_block (file_crop, 'pdrice2', grid_crop, f_xy_crop)
      ENDIF

      CALL mg2patch_crop%grid2pset (f_xy_crop, pdrice2_tmp)

      IF (p_is_worker) THEN
         DO npatch = 1, numpatch
            IF (pdrice2_tmp(npatch) /= spval) THEN
               pdrice2 (npatch) = int(pdrice2_tmp (npatch))
            ELSE
               pdrice2 (npatch) = 0
            ENDIF
         ENDDO
      ENDIF

#ifdef RangeCheck
      CALL check_vector_data ('plant date value for rice2 ', pdrice2)
#endif

      ! (2) Read in plant date.
      IF (p_is_worker) THEN
         IF (numpft > 0) plantdate_p(:) = -99999999._r8
      ENDIF

      file_crop = trim(DEF_dir_runtime) // '/crop/plantdt-colm-64cfts-rice2_fillcoast.nc'
      DO cft = 15, 78
         write(cx, '(i2.2)') cft
         IF (p_is_io) THEN
            CALL ncio_read_block_time (file_crop, &
               'PLANTDATE_CFT_'//trim(cx), grid_crop, 1, f_xy_crop)
         ENDIF

         CALL mg2pft_crop%grid2pset (f_xy_crop, plantdate_tmp)

         IF (p_is_worker) THEN
            DO ipft = 1, numpft
               IF(landpft%settyp(ipft) .eq. cft)THEN
                  plantdate_p(ipft) = plantdate_tmp(ipft)
                  IF(plantdate_p(ipft) <= 0._r8) THEN
                     plantdate_p(ipft) = -99999999._r8
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
      ENDDO

#ifdef RangeCheck
      CALL check_vector_data ('plantdate_pfts value ', plantdate_p)
#endif

      ! (3) Read in fertlization
      IF (DEF_FERT_SOURCE == 1) THEN
         IF (p_is_worker) THEN
            IF (numpft > 0) fertnitro_p(:) = -99999999._r8
         ENDIF

         file_fert = trim(DEF_dir_runtime) // '/crop/fertnitro_fillcoast.nc'
         DO cft = 15, 78
            write(cx, '(i2.2)') cft
            IF (p_is_io) THEN
               CALL ncio_read_block_time (file_fert, &
                  'CONST_FERTNITRO_CFT_'//trim(cx), grid_crop, 1, f_xy_crop)
            ENDIF

            CALL mg2pft_crop%grid2pset (f_xy_crop, fertnitro_tmp)

            IF (p_is_worker) THEN
               DO ipft = 1, numpft
                  IF(landpft%settyp(ipft) .eq. cft)THEN
                     fertnitro_p(ipft) = fertnitro_tmp(ipft)
                     IF(fertnitro_p(ipft) <= 0._r8) THEN
                        fertnitro_p(ipft) = 0._r8
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDDO

      ELSEIF (DEF_FERT_SOURCE == 2)THEN
         IF (p_is_worker) THEN
            IF (numpft > 0) THEN
               allocate(fertilizer_tmp(numpft))
               allocate(manure_tmp(numpft))
            ENDIF
         ENDIF

         file_fert = trim(DEF_dir_runtime) // '/crop/fertilizer_2015soc.nc'

         CALL ncio_read_bcast_serial (file_fert, 'lat', lat)
         CALL ncio_read_bcast_serial (file_fert, 'lon', lon)
         CALL grid_fert%define_by_center (lat, lon)
         
         IF (p_is_io) THEN
            CALL allocate_block_data (grid_fert, f_xy_fert)
         ENDIF

         CALL mg2pft_fert%build_arealweighted (grid_fert, landpft)

         IF (allocated(lon)) deallocate(lon)
         IF (allocated(lat)) deallocate(lat)

         IF (p_is_worker) THEN
            fertnitro_p(:) = -99999999._r8
            manunitro_p(:) = -99999999._r8
         ENDIF
         
         ! read manure
         IF (p_is_io) THEN
            CALL ncio_read_block (file_fert, 'manure', grid_fert, f_xy_fert)
         ENDIF

         CALL mg2pft_fert%grid2pset (f_xy_fert, manure_tmp)

         IF (p_is_worker) THEN
            DO ipft = 1, numpft
               IF (landpft%settyp(ipft) .ge. 15) THEN
                  manunitro_p(ipft) = manure_tmp(ipft)
                  IF (manunitro_p(ipft) < 0._r8) THEN
                     manunitro_p(ipft) = 0._r8
                  ENDIF
               ENDIF
            ENDDO
         ENDIF

         ! read fertilizer
         DO cft = 1, 64
            IF (p_is_io) THEN
               CALL ncio_read_block_time (file_fert, 'fertilizer', grid_fert, cft, f_xy_fert)
            ENDIF

            CALL mg2pft_fert%grid2pset (f_xy_fert, fertilizer_tmp)

            IF (p_is_worker) THEN
               DO ipft = 1, numpft
                  IF (landpft%settyp(ipft) .eq. cft+14) THEN
                     fertnitro_p(ipft) = fertilizer_tmp(ipft)
                     IF (fertnitro_p(ipft) < 0._r8) THEN
                        fertnitro_p(ipft) = 0._r8
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDDO

         IF (allocated (fertilizer_tmp)) deallocate(fertilizer_tmp)
         IF (allocated (manure_tmp))     deallocate(manure_tmp)
      ENDIF

#ifdef RangeCheck
      CALL check_vector_data ('fert nitro value ', fertnitro_p)
#endif

      ! (4) Read in irrigation method
      !file_irrig = trim(DEF_dir_runtime) // '/crop/surfdata_irrigation_method.nc'
      file_irrig = trim(DEF_dir_runtime) // '/crop/surfdata_irrigation_method_96x144.nc'

      CALL ncio_read_bcast_serial (file_irrig, 'lat', lat)
      CALL ncio_read_bcast_serial (file_irrig, 'lon', lon)

      CALL grid_irrig%define_by_center (lat, lon)

      IF (p_is_io) THEN
         CALL allocate_block_data  (grid_irrig, f_xy_irrig)
      ENDIF

      CALL mg2pft_irrig%build_arealweighted (grid_irrig, landpft)

      IF (allocated(lon)) deallocate(lon)
      IF (allocated(lat)) deallocate(lat)

      IF (p_is_worker) THEN
         IF (numpft > 0) irrig_method_p(:) = -99999999
      ENDIF

      DO cft = 1, N_CFT
         IF (p_is_io) THEN
            CALL ncio_read_block_time (file_irrig, 'irrigation_method', grid_irrig, cft, f_xy_irrig)
         ENDIF

         CALL mg2pft_irrig%grid2pset_dominant (f_xy_irrig, irrig_method_tmp)

         IF (p_is_worker) THEN
            DO ipft = 1, numpft

               IF(landpft%settyp(ipft) .eq. cft + 14)THEN
                  irrig_method_p(ipft) = irrig_method_tmp(ipft)
                  IF(irrig_method_p(ipft) < 0) THEN
                     irrig_method_p(ipft) = -99999999
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
      ENDDO

#ifdef RangeCheck
      CALL check_vector_data ('irrigation method ', irrig_method_p)
#endif

      IF (allocated (pdrice2_tmp  ))    deallocate(pdrice2_tmp  )
      IF (allocated (plantdate_tmp))    deallocate(plantdate_tmp)
      IF (allocated (fertnitro_tmp))    deallocate(fertnitro_tmp)
      IF (allocated (irrig_method_tmp)) deallocate(irrig_method_tmp)

      IF (DEF_IRRIGATION_ALLOCATION == 3) THEN 
         ! (5) Read in irrigation allocated to groundwater
         file_irrigalloc = trim(DEF_dir_runtime) // '/crop/surfdata_irrigation_allocation.nc'

         CALL ncio_read_bcast_serial (file_irrigalloc, 'lat', lat)
         CALL ncio_read_bcast_serial (file_irrigalloc, 'lon', lon)

         CALL grid_irrigalloc%define_by_center (lat, lon)

         IF (p_is_io) THEN
            CALL allocate_block_data (grid_irrigalloc, f_xy_irrigalloc)
         ENDIF

         call mg2p_irrigalloc%build_arealweighted (grid_irrigalloc, landpatch)

         IF (allocated(lon)) deallocate(lon)
         IF (allocated(lat)) deallocate(lat)

         IF (p_is_worker) THEN
            irrig_gw_alloc(:) = 0._r8
         ENDIF

         IF (p_is_io) THEN
            CALL ncio_read_block (file_irrigalloc, 'irrig_gw_alloc', grid_irrigalloc, f_xy_irrigalloc)
         ENDIF

         call mg2p_irrigalloc%grid2pset (f_xy_irrigalloc, irrig_gw_alloc)

#ifdef RangeCheck
         CALL check_vector_data ('irrigation goundwater allocation ', irrig_gw_alloc)
#endif

         ! (6) Read in irrigation allocated to surfacewater
         IF (p_is_worker) THEN
            irrig_sw_alloc(:) = 0._r8
         ENDIF

         IF (p_is_io) THEN
            CALL ncio_read_block (file_irrigalloc, 'irrig_sw_alloc', grid_irrigalloc, f_xy_irrigalloc)
         ENDIF

         call mg2p_irrigalloc%grid2pset (f_xy_irrigalloc, irrig_sw_alloc)

#ifdef RangeCheck
         CALL check_vector_data ('irrigation surfacewater allocation ', irrig_sw_alloc)
#endif
      ENDIF


      ! (5) Read in soil pH
      file_acid = trim(DEF_dir_runtime) // '/soil/PHH2O_16800x43200.nc'

      CALL ncio_read_bcast_serial (file_acid, 'lat', lat)
      CALL ncio_read_bcast_serial (file_acid, 'lon', lon)

      CALL grid_acid%define_by_center (lat, lon)

      IF (p_is_io) THEN
         CALL allocate_block_data  (grid_acid, f_xy_acid)
      ENDIF

      IF (p_is_master) THEN
         CALL ncio_get_attr(file_acid, 'PHH2O', 'missing_value', missing_value)
      ENDIF
#ifdef USEMPI
      CALL mpi_bcast(missing_value, 1, MPI_REAL8, p_address_master, p_comm_glb, p_err)
#endif

      CALL mg2pft_acid%build_arealweighted (grid_acid, landpatch)
   
      IF (allocated(lon)) deallocate(lon)
      IF (allocated(lat)) deallocate(lat)

      IF (p_is_worker) THEN
         IF (numpatch > 0)  allocate(pH_vr_tmp(numpatch))
      ENDIF

      DO nsl = 1, 8
         IF (p_is_io) THEN
            CALL ncio_read_block_time (file_acid, 'PHH2O', grid_acid, nsl, f_xy_acid)
         ENDIF

         CALL mg2pft_acid%set_missing_value(f_xy_acid, missing_value)

         CALL mg2pft_acid%grid2pset (f_xy_acid, pH_vr_tmp)

         IF (p_is_worker) THEN
            DO npatch = 1, numpatch
               IF (pH_vr_tmp(npatch) /= spval) THEN
                  pH_vr(nsl,npatch) = pH_vr_tmp(npatch)/10._r8
                  pH_init_vr(nsl,npatch) = pH_vr_tmp(npatch)/10._r8
               ELSE
                  pH_vr(nsl,npatch) = spval
                  pH_init_vr(nsl,npatch) = spval
               ENDIF
            ENDDO
         ENDIF
      ENDDO

      IF (p_is_worker) THEN
         DO nsl = 9, 2, -1
            pH_vr(nsl,:) = pH_vr(nsl-1,:)
            pH_init_vr(nsl,:) = pH_init_vr(nsl-1,:)
         ENDDO

         DO nsl = nl_soil, 10, -1
            pH_vr(nsl,:) = pH_vr(9,:)
            pH_init_vr(nsl,:) = pH_init_vr(9,:)
         ENDDO
      ENDIF

      IF (allocated (pH_vr_tmp))              deallocate(pH_vr_tmp)


#ifdef RangeCheck
      DO nsl = 1, nl_soil
         CALL check_vector_data ('Soil pH layer '//trim(str(nsl))//' ', pH_vr(nsl,:))
         CALL check_vector_data ('Soil init pH layer '//trim(str(nsl))//' ', pH_init_vr(nsl,:))
      END DO
#endif

      ! (6) Read in soil CACO3 %
      file_CACO3 = trim(DEF_dir_runtime) // '/soil/CACO3_16800x43200.nc'

      CALL ncio_read_bcast_serial (file_CACO3, 'lat', lat)
      CALL ncio_read_bcast_serial (file_CACO3, 'lon', lon)

      CALL grid_CACO3%define_by_center (lat, lon)

      IF (p_is_io) THEN
         CALL allocate_block_data  (grid_CACO3, f_xy_CACO3)
      ENDIF

      IF (p_is_master) THEN
         CALL ncio_get_attr(file_CACO3, 'CACO3', 'missing_value', missing_value)
      ENDIF
#ifdef USEMPI
      CALL mpi_bcast(missing_value, 1, MPI_REAL8, p_address_master, p_comm_glb, p_err)
#endif

      CALL mg2pft_CACO3%build_arealweighted (grid_CACO3, landpatch)
   
      IF (allocated(lon)) deallocate(lon)
      IF (allocated(lat)) deallocate(lat)

      IF (p_is_worker) THEN
         IF (numpatch > 0)  allocate(calcite_content_tmp(numpatch))
      ENDIF

      DO nsl = 1, 8
         IF (p_is_io) THEN
            CALL ncio_read_block_time (file_CACO3, 'CACO3', grid_CACO3, nsl, f_xy_CACO3)
         ENDIF

         CALL mg2pft_CACO3%set_missing_value(f_xy_CACO3, missing_value)

         CALL mg2pft_CACO3%grid2pset (f_xy_CACO3, calcite_content_tmp)

         IF (p_is_worker) THEN
            DO npatch = 1, numpatch
               IF (calcite_content_tmp(npatch) /= spval) THEN
                  calcite_content(nsl,npatch) = calcite_content_tmp(npatch)/100_r8 !%
               ELSE
                  calcite_content(nsl,npatch) = spval
               ENDIF
            ENDDO
         ENDIF
      ENDDO

      IF (p_is_worker) THEN
         DO nsl = 9, 2, -1
            calcite_content(nsl,:) = calcite_content(nsl-1,:)
         ENDDO

         DO nsl = nl_soil, 10, -1
            calcite_content(nsl,:) = calcite_content(9,:)
         ENDDO
      ENDIF

      IF (allocated (calcite_content_tmp))              deallocate(calcite_content_tmp)

#ifdef RangeCheck
      DO nsl = 1, nl_soil
         CALL check_vector_data ('Soil calcite_content layer '//trim(str(nsl))//' ', calcite_content(nsl,:))
      END DO
#endif

      ! (7) Read in soil CEC
      file_CEC = trim(DEF_dir_runtime) // '/soil/CEC_16800x43200.nc'

      CALL ncio_read_bcast_serial (file_CEC, 'lat', lat)
      CALL ncio_read_bcast_serial (file_CEC, 'lon', lon)

      CALL grid_CEC%define_by_center (lat, lon)

      IF (p_is_io) THEN
         CALL allocate_block_data  (grid_CEC, f_xy_CEC)
      ENDIF

      IF (p_is_master) THEN
         CALL ncio_get_attr(file_CEC, 'CEC', 'missing_value', missing_value)
      ENDIF
#ifdef USEMPI
      CALL mpi_bcast(missing_value, 1, MPI_REAL8, p_address_master, p_comm_glb, p_err)
#endif

      CALL mg2pft_CEC%build_arealweighted (grid_CEC, landpatch)
   
      IF (allocated(lon)) deallocate(lon)
      IF (allocated(lat)) deallocate(lat)

      IF (p_is_worker) THEN
         IF (numpatch > 0)  allocate(CEC_tmp(numpatch))
      ENDIF

      DO nsl = 1, 8
         IF (p_is_io) THEN
            CALL ncio_read_block_time (file_CEC, 'CEC', grid_CEC, nsl, f_xy_CEC)
         ENDIF

         CALL mg2pft_CEC%set_missing_value(f_xy_CEC, missing_value)

         CALL mg2pft_CEC%grid2pset (f_xy_CEC, CEC_tmp)

         IF (p_is_worker) THEN
            DO npatch = 1, numpatch
               IF (CEC_tmp(npatch) /= spval) THEN
                  CEC(nsl,npatch) = CEC_tmp(npatch)/100_r8
               ELSE
                  CEC(nsl,npatch) = spval
               ENDIF
            ENDDO
         ENDIF
      ENDDO

      IF (p_is_worker) THEN
         DO nsl = 9, 2, -1
            CEC(nsl,:) = CEC(nsl-1,:)
         ENDDO

         DO nsl = nl_soil, 10, -1
            CEC(nsl,:) = CEC(9,:)
         ENDDO
      ENDIF

      IF (allocated (CEC_tmp))              deallocate(CEC_tmp)

#ifdef RangeCheck
      DO nsl = 1, nl_soil
         CALL check_vector_data ('Soil CEC layer '//trim(str(nsl))//' ', CEC(nsl,:))
      END DO
#endif

      ! (8) Read in soil E_Bc
      file_BS = trim(DEF_dir_runtime) // '/soil/BS_16800x43200.nc'

      CALL ncio_read_bcast_serial (file_BS, 'lat', lat)
      CALL ncio_read_bcast_serial (file_BS, 'lon', lon)

      CALL grid_BS%define_by_center (lat, lon)

      IF (p_is_io) THEN
         CALL allocate_block_data  (grid_BS, f_xy_BS)
      ENDIF

      IF (p_is_master) THEN
         CALL ncio_get_attr(file_BS, 'BS', 'missing_value', missing_value)
      ENDIF
#ifdef USEMPI
      CALL mpi_bcast(missing_value, 1, MPI_REAL8, p_address_master, p_comm_glb, p_err)
#endif

      CALL mg2pft_BS%build_arealweighted (grid_BS, landpatch)
   
      IF (allocated(lon)) deallocate(lon)
      IF (allocated(lat)) deallocate(lat)

      IF (p_is_worker) THEN
         IF (numpatch > 0)  allocate(E_Bc_tmp(numpatch))
      ENDIF

      DO nsl = 1, 8
         IF (p_is_io) THEN
            CALL ncio_read_block_time (file_BS, 'BS', grid_BS, nsl, f_xy_BS)
         ENDIF

         CALL mg2pft_BS%set_missing_value(f_xy_BS, missing_value)

         CALL mg2pft_BS%grid2pset (f_xy_BS, E_Bc_tmp)

         IF (p_is_worker) THEN
            DO npatch = 1, numpatch
               IF (E_Bc_tmp(npatch) /= spval) THEN
                  E_Bc(nsl,npatch) = E_Bc_tmp(npatch)
               ELSE
                  E_Bc(nsl,npatch) = spval
               ENDIF
            ENDDO
         ENDIF
      ENDDO

      IF (p_is_worker) THEN
         DO nsl = 9, 2, -1
            E_Bc(nsl,:) = E_Bc(nsl-1,:)
         ENDDO

         DO nsl = nl_soil, 10, -1
            E_Bc(nsl,:) = E_Bc(9,:)
         ENDDO
      ENDIF

      IF (allocated (E_Bc_tmp))              deallocate(E_Bc_tmp)

#ifdef RangeCheck
      DO nsl = 1, nl_soil
         CALL check_vector_data ('Soil E_Bc layer '//trim(str(nsl))//' ', E_Bc(nsl,:))
      END DO
#endif

   END SUBROUTINE CROP_readin

END MODULE MOD_CropReadin
#endif
