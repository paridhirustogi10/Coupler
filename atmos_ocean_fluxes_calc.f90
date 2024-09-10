module atmos_ocean_fluxes_calc_mod
  use coupler_types_mod, only : coupler_1d_bc_type,&
      & ind_flux,&
      & ind_kw,&
      & ind_u10,&
      & ind_ustar,&
      & ind_hs,&
      & ind_alpha,&
      & ind_pCair,&
      & ind_psurf,&
      & ind_sc_no,&
      & ind_csurf,&
      & ind_flux0,&
      & ind_deltap,&
      & ind_deposition
  use mpp_mod,           only : mpp_error, FATAL
  use constants_mod,     only : wtmair, rdgas, vonkarm
  
  #ifdef INTERNAL_FILE_NML
  use          mpp_mod, only: input_nml_file
  #else
  use          fms_mod, only: open_namelist_file
  #endif
  use        fms_mod, only: close_file, &
                          error_mesg, file_exist, check_nml_error, FATAL, &
                          mpp_pe, mpp_root_pe, &
                          write_version_number, stdlog
  implicit none
  private

  public atmos_ocean_fluxes_calc

  character(len=*), parameter :: mod_name = "cdwfe"

  real, parameter :: epsln=1.0e-30

  character(len=256) :: version = '$Id$'
  character(len=256) :: tagname = '$Name$'
  integer :: kw_dic = 0
  namelist /kw_nml/  kw_dic

contains
  !> \brief Calculate the ocean gas fluxes. Units should be mol/m^2/s, upward flux is positive.
  !
  !! \throw FATAL, "Number of gas fluxes not zero"
  !! \throw FATAL, "Lengths of flux fields do not match"
  !! \throw FATAL, "Unknown implementation ([implementation]) for [name]"
  !! \throw FATAL, "Lengths of flux fields do not match"
  !! \throw FATAL, "Bad parameter ([gas_fluxes%bc(n)%param(1)]) for land_sea_runoff for [name]"
  !! \throw FATAL, "Unknown flux type ([flux_type]) for [name]"
  subroutine atmos_ocean_fluxes_calc(gas_fields_atm, gas_fields_ice,&
      & gas_fluxes, seawater, tsurf, ustar, cd_m)
    type(coupler_1d_bc_type), intent(in)    :: gas_fields_atm !< Structure containing atmospheric surface
                                                              !! variables that are used in the calculation
                                                              !! of the atmosphere-ocean gas fluxes.
    type(coupler_1d_bc_type), intent(in)    :: gas_fields_ice !< Structure containing ice-top and ocean
                                                              !! surface variables that are used in the
                                                              !! calculation of the atmosphere-ocean gas fluxes.
    type(coupler_1d_bc_type), intent(inout) :: gas_fluxes !< Structure containing the gas fluxes between
                                                          !! the atmosphere and the ocean and parameters
                                                          !! related to the calculation of these fluxes.
    real, dimension(:), intent(in)           :: seawater  !< 1 for the open water category, 0 if ice or land.
    real, dimension(:), intent(in)           :: tsurf     !< SST in units of K
    real, dimension(:), intent(in), optional :: ustar, cd_m

    character(len=*), parameter   :: sub_name = 'atmos_ocean_fluxes_calc'
    character(len=*), parameter   :: error_header =&
        & '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

    integer                                 :: n
    integer                                 :: i
    integer                                 :: length
    real, dimension(:), allocatable         :: kw
    real, dimension(:), allocatable         :: cair
    character(len=128)                      :: error_string

    real, parameter :: permeg=1.0e-6
    integer :: io, ierr, unit
    !integer :: kw_dic
    !kwcalc = 2 ! change to be 1 as requested by Brandon, original value is 0 
               ! How to make this function as a run-time parameter?

!------- read namelist for kw_nml ---------
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=kw_nml, iostat=io)
  ierr = check_nml_error(io, 'kw_nml')
#else
      if (file_exist('input.nml')) then
         unit = open_namelist_file ('input.nml')
         ierr=1; do while (ierr /= 0)
            read  (unit, nml=kw_nml, iostat=io, end=10)
            ierr = check_nml_error(io,'kw_nml')
         enddo
  10     call close_file (unit)
      endif
#endif
!------- write version number and namelist ---------
      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number(version, tagname)
           unit = stdlog()
           write (unit, nml=kw_nml)
      endif
    ! Return if no fluxes to be calculated
    if (gas_fluxes%num_bcs .le. 0) return

    if (.not. associated(gas_fluxes%bc)) then
      if (gas_fluxes%num_bcs .ne. 0) then
        call mpp_error(FATAL, trim(error_header) // ' Number of gas fluxes not zero')
      else
        return
      endif
    endif

    do n = 1, gas_fluxes%num_bcs
      ! only do calculations if the flux has not been overridden
      if ( .not. gas_fluxes%bc(n)%field(ind_flux)%override) then
        if (gas_fluxes%bc(n)%flux_type .eq. 'air_sea_gas_flux_generic') then
          length = size(gas_fluxes%bc(n)%field(1)%values(:))

          if (.not. allocated(kw)) then
            allocate( kw(length) )
            allocate ( cair(length) )
          elseif (size(kw(:)) .ne. length) then
            call mpp_error(FATAL, trim(error_header) // ' Lengths of flux fields do not match')
          endif

          if (gas_fluxes%bc(n)%implementation .eq. 'ocmip2') then
            do i = 1, length
              if (seawater(i) == 1.) then
                if (trim(gas_fluxes%bc(n)%name) == 'co2_flux') then
                   if (kw_dic==0) then
                      gas_fluxes%bc(n)%field(ind_kw)%values(i) =&
                           gas_fluxes%bc(n)%param(1) * gas_fields_atm%bc(n)%field(ind_u10)%values(i)**2
                   elseif (kw_dic==1) then
                      gas_fluxes%bc(n)%field(ind_kw)%values(i) =&
                          DM18_mean (gas_fields_atm%bc(n)%field(ind_u10)%values(i))
                   elseif (kw_dic==2) then
                      if (gas_fields_atm%bc(n)%field(ind_hs)%values(i)>0.5.and.&
                          gas_fields_atm%bc(n)%field(ind_hs)%values(i)<30.) then
                    ! this basically says if Hs is not >0.5 m or <39.5 m we
                    ! shift to DM18_mean
                    ! Luc and Brandon found that a factor of 0.775 times
                    ! the coefficients from his paper was needed to match the
                    ! Wanninkhof global mean value (evaluated from JRA55 data
                    ! with the Wanninkhof formula). 0.775 was already applied in
                    ! DM18_mean
                         gas_fluxes%bc(n)%field(ind_kw)%values(i) =0.775*&!(Ustar,HS,SC,Alpha)
                            DM18 (gas_fields_atm%bc(n)%field(ind_ustar)%values(i),&
                                   gas_fields_atm%bc(n)%field(ind_hs)%values(i),&
                                   gas_fields_ice%bc(n)%field(ind_sc_no)%values(i),&
                                   gas_fields_ice%bc(n)%field(ind_alpha)%values(i),&
                                   tsurf(i))
                      else
                         gas_fluxes%bc(n)%field(ind_kw)%values(i) =&
                            DM18_mean (gas_fields_atm%bc(n)%field(ind_u10)%values(i))
                      endif
                   endif !kw_dic
                else
                   gas_fluxes%bc(n)%field(ind_kw)%values(i) =&
                           gas_fluxes%bc(n)%param(1) *gas_fields_atm%bc(n)%field(ind_u10)%values(i)**2
                endif !gas_fluxes%bc(n)%name

                cair(i) = &
                    gas_fields_ice%bc(n)%field(ind_alpha)%values(i) * &
                    gas_fields_atm%bc(n)%field(ind_pCair)%values(i) * &
                    gas_fields_atm%bc(n)%field(ind_psurf)%values(i) * gas_fluxes%bc(n)%param(2)
                gas_fluxes%bc(n)%field(ind_flux)%values(i) =&
                    & gas_fluxes%bc(n)%field(ind_kw)%values(i) *&
                    & sqrt(660. / (gas_fields_ice%bc(n)%field(ind_sc_no)%values(i) + epsln)) *&
                    & (gas_fields_ice%bc(n)%field(ind_csurf)%values(i) - cair(i))
                gas_fluxes%bc(n)%field(ind_flux0)%values(i) =&
                    & gas_fluxes%bc(n)%field(ind_kw)%values(i) *&
                    & sqrt(660. / (gas_fields_ice%bc(n)%field(ind_sc_no)%values(i) + epsln)) *&
                    & gas_fields_ice%bc(n)%field(ind_csurf)%values(i)
                gas_fluxes%bc(n)%field(ind_deltap)%values(i) =&
                    & (gas_fields_ice%bc(n)%field(ind_csurf)%values(i) - cair(i)) / &
                   (gas_fields_ice%bc(n)%field(ind_alpha)%values(i) * permeg + epsln)
              else
                gas_fluxes%bc(n)%field(ind_kw)%values(i) = 0.0
                gas_fluxes%bc(n)%field(ind_flux)%values(i) = 0.0
                gas_fluxes%bc(n)%field(ind_flux0)%values(i) = 0.0
                gas_fluxes%bc(n)%field(ind_deltap)%values(i) = 0.0
                cair(i) = 0.0
              endif
            enddo
          elseif (gas_fluxes%bc(n)%implementation .eq. 'duce') then
            do i = 1, length
              if (seawater(i) == 1.) then
                gas_fluxes%bc(n)%field(ind_kw)%values(i) = &
                    & gas_fields_atm%bc(n)%field(ind_u10)%values(i) /&
                    & (770.+45.*gas_fluxes%bc(n)%param(1)**(1./3.)) *&
                    & 101325./(rdgas*wtmair*1e-3*tsurf(i) *&
                    & gas_fields_ice%bc(n)%field(ind_alpha)%values(i))
                !alpha: mol/m3/atm
                cair(i) = &
                    gas_fields_ice%bc(n)%field(ind_alpha)%values(i) * &
                    gas_fields_atm%bc(n)%field(ind_pCair)%values(i) * &
                    gas_fields_atm%bc(n)%field(ind_psurf)%values(i) * 9.86923e-6
                cair(i) = max(cair(i),0.)
                gas_fluxes%bc(n)%field(ind_flux)%values(i) =&
                    & gas_fluxes%bc(n)%field(ind_kw)%values(i) *&
                    & (max(gas_fields_ice%bc(n)%field(ind_csurf)%values(i),0.) - cair(i))
                gas_fluxes%bc(n)%field(ind_flux0)%values(i) =&
                    & gas_fluxes%bc(n)%field(ind_kw)%values(i) *&
                    & max(gas_fields_ice%bc(n)%field(ind_csurf)%values(i),0.)
                gas_fluxes%bc(n)%field(ind_deltap)%values(i) =&
                    & (max(gas_fields_ice%bc(n)%field(ind_csurf)%values(i),0.) - cair(i)) /&
                    & (gas_fields_ice%bc(n)%field(ind_alpha)%values(i) * permeg + epsln)
              else
                gas_fluxes%bc(n)%field(ind_kw)%values(i) = 0.0
                gas_fluxes%bc(n)%field(ind_flux)%values(i) = 0.0
                gas_fluxes%bc(n)%field(ind_flux0)%values(i) = 0.0
                gas_fluxes%bc(n)%field(ind_deltap)%values(i) = 0.0
                cair(i) = 0.0
              endif
            enddo
          elseif (gas_fluxes%bc(n)%implementation .eq. 'johnson') then
            !f1p: not sure how to pass salinity. For now, just force at 35.
            do i = 1, length
              if (seawater(i) == 1.) then
                !calc_kw(tk,p,u10,h,vb,mw,sc_w,ustar,cd_m)
                gas_fluxes%bc(n)%field(ind_kw)%values(i) =&
                    & calc_kw(tsurf(i),&
                    & gas_fields_atm%bc(n)%field(ind_psurf)%values(i),&
                    & gas_fields_atm%bc(n)%field(ind_u10)%values(i),&
                    & 101325./(rdgas*wtmair*1e-3*tsurf(i)*gas_fields_ice%bc(n)%field(ind_alpha)%values(i)),&
                    & gas_fluxes%bc(n)%param(2),&
                    & gas_fluxes%bc(n)%param(1),&
                    & gas_fields_ice%bc(n)%field(ind_sc_no)%values(i))
                cair(i) =&
                    & gas_fields_ice%bc(n)%field(ind_alpha)%values(i) *&
                    & gas_fields_atm%bc(n)%field(ind_pCair)%values(i) *&
                    & gas_fields_atm%bc(n)%field(ind_psurf)%values(i) * 9.86923e-6
                cair(i) = max(cair(i),0.)
                gas_fluxes%bc(n)%field(ind_flux)%values(i) =&
                    & gas_fluxes%bc(n)%field(ind_kw)%values(i) *&
                    & (max(gas_fields_ice%bc(n)%field(ind_csurf)%values(i),0.) - cair(i))
                gas_fluxes%bc(n)%field(ind_flux0)%values(i) =&
                    & gas_fluxes%bc(n)%field(ind_kw)%values(i) *&
                    & max(gas_fields_ice%bc(n)%field(ind_csurf)%values(i),0.)
                gas_fluxes%bc(n)%field(ind_deltap)%values(i) =&
                    & (max(gas_fields_ice%bc(n)%field(ind_csurf)%values(i),0.) - cair(i)) /&
                    & (gas_fields_ice%bc(n)%field(ind_alpha)%values(i) * permeg + epsln)
              else
                gas_fluxes%bc(n)%field(ind_kw)%values(i) = 0.0
                gas_fluxes%bc(n)%field(ind_flux)%values(i) = 0.0
                gas_fluxes%bc(n)%field(ind_flux0)%values(i) = 0.0
                gas_fluxes%bc(n)%field(ind_deltap)%values(i) = 0.0
                cair(i) = 0.0
              endif
            enddo
          else
            call mpp_error(FATAL, ' Unknown implementation (' //&
                & trim(gas_fluxes%bc(n)%implementation) // ') for ' // trim(gas_fluxes%bc(n)%name))
          endif
        elseif (gas_fluxes%bc(n)%flux_type .eq. 'air_sea_gas_flux') then
          length = size(gas_fluxes%bc(n)%field(1)%values(:))

          if (.not. allocated(kw)) then
            allocate( kw(length) )
            allocate ( cair(length) )
          elseif (size(kw(:)) .ne. length) then
            call mpp_error(FATAL, trim(error_header) // ' Lengths of flux fields do not match')
          endif

          if (gas_fluxes%bc(n)%implementation .eq. 'ocmip2_data') then
            do i = 1, length
              if (seawater(i) == 1.) then
                kw(i) = gas_fluxes%bc(n)%param(1) * gas_fields_atm%bc(n)%field(ind_u10)%values(i)
                cair(i) =&
                    & gas_fields_ice%bc(n)%field(ind_alpha)%values(i) *&
                    & gas_fields_atm%bc(n)%field(ind_pCair)%values(i) *&
                    & gas_fields_atm%bc(n)%field(ind_psurf)%values(i) * gas_fluxes%bc(n)%param(2)
                gas_fluxes%bc(n)%field(ind_flux)%values(i) = kw(i) *&
                    & (gas_fields_ice%bc(n)%field(ind_csurf)%values(i) - cair(i))
              else
                gas_fluxes%bc(n)%field(ind_flux)%values(i) = 0.0
                cair(i) = 0.0
                kw(i) = 0.0
              endif
            enddo
          elseif (gas_fluxes%bc(n)%implementation .eq. 'ocmip2') then
            do i = 1, length
              if (seawater(i) == 1.) then
                kw(i) = gas_fluxes%bc(n)%param(1) * gas_fields_atm%bc(n)%field(ind_u10)%values(i)**2
                cair(i) =&
                    & gas_fields_ice%bc(n)%field(ind_alpha)%values(i) *&
                    & gas_fields_atm%bc(n)%field(ind_pCair)%values(i) *&
                    & gas_fields_atm%bc(n)%field(ind_psurf)%values(i) * gas_fluxes%bc(n)%param(2)
                gas_fluxes%bc(n)%field(ind_flux)%values(i) = kw(i) *&
                    & (gas_fields_ice%bc(n)%field(ind_csurf)%values(i) - cair(i))
              else
                gas_fluxes%bc(n)%field(ind_flux)%values(i) = 0.0
                cair(i) = 0.0
                kw(i) = 0.0
              endif
            enddo
          elseif (gas_fluxes%bc(n)%implementation .eq. 'linear') then
            do i = 1, length
              if (seawater(i) == 1.) then
                kw(i) = gas_fluxes%bc(n)%param(1) *&
                    & max(0.0, gas_fields_atm%bc(n)%field(ind_u10)%values(i) - gas_fluxes%bc(n)%param(2))
                cair(i) =&
                    & gas_fields_ice%bc(n)%field(ind_alpha)%values(i) *&
                    & gas_fields_atm%bc(n)%field(ind_pCair)%values(i) *&
                    & gas_fields_atm%bc(n)%field(ind_psurf)%values(i) * gas_fluxes%bc(n)%param(3)
                gas_fluxes%bc(n)%field(ind_flux)%values(i) = kw(i) *&
                    & (gas_fields_ice%bc(n)%field(ind_csurf)%values(i) - cair(i))
              else
                gas_fluxes%bc(n)%field(ind_flux)%values(i) = 0.0
                cair(i) = 0.0
                kw(i) = 0.0
              endif
            enddo
          else
            call mpp_error(FATAL, ' Unknown implementation (' //&
                & trim(gas_fluxes%bc(n)%implementation) // ') for ' // trim(gas_fluxes%bc(n)%name))
          endif
        elseif (gas_fluxes%bc(n)%flux_type .eq. 'air_sea_deposition') then
          cycle !air_sea_deposition is done in another subroutine
        elseif (gas_fluxes%bc(n)%flux_type .eq. 'land_sea_runoff') then
          if (gas_fluxes%bc(n)%param(1) .le. 0.0) then
            write (error_string, '(1pe10.3)') gas_fluxes%bc(n)%param(1)
            call mpp_error(FATAL, ' Bad parameter (' // trim(error_string) //&
                & ') for land_sea_runoff for ' // trim(gas_fluxes%bc(n)%name))
          endif

          length = size(gas_fluxes%bc(n)%field(1)%values(:))

          if (gas_fluxes%bc(n)%implementation .eq. 'river') then
            do i = 1, length
              if (seawater(i) == 1.) then
                gas_fluxes%bc(n)%field(ind_flux)%values(i) =&
                    & gas_fields_atm%bc(n)%field(ind_deposition)%values(i) /&
                    & gas_fluxes%bc(n)%param(1)
              else
                gas_fluxes%bc(n)%field(ind_flux)%values(i) = 0.0
              endif
            enddo
          else
            call mpp_error(FATAL, ' Unknown implementation (' //&
                & trim(gas_fluxes%bc(n)%implementation) // ') for ' // trim(gas_fluxes%bc(n)%name))
          endif
        else
          call mpp_error(FATAL, ' Unknown flux_type (' // trim(gas_fluxes%bc(n)%flux_type) //&
              & ') for ' // trim(gas_fluxes%bc(n)%name))
        endif
      endif
    enddo

    if (allocated(kw)) then
      deallocate(kw)
      deallocate(cair)
    endif
  end subroutine  atmos_ocean_fluxes_calc

  !> Calculate \f$k_w\f$
  !!
  !! Taken from Johnson, Ocean Science, 2010. (http://doi.org/10.5194/os-6-913-2010)
  !!
  !! Uses equations defined in Liss[1974],
  !! \f[
  !!  F = K_g(c_g - H C_l) = K_l(c_g/H - C_l)
  !! \f]
  !! where \f$c_g\f$ and \f$C_l\f$ are the bulk gas and liquid concentrations, \f$H\f$
  !! is the Henry's law constant (\f$H = c_{sg}/C_{sl}\f$, where \f$c_{sg}\f$ is the
  !! equilibrium concentration in gas phase (\f$g/cm^3\f$ of air) and \f$C_{sl}\f$ is the
  !! equilibrium concentration of unionised dissolved gas in liquid phase (\f$g/cm^3\f$
  !! of water)),
  !! \f[
  !!    1/K_g = 1/k_g + H/k_l
  !! \f]
  !! and
  !! \f[
  !!    1/K_l = 1/k_l + 1/{Hk_g}
  !! \f]
  !! where \f$k_g\f$ and \f$k_l\f$ are the exchange constants for the gas and liquid
  !! phases, respectively.
  real function calc_kw(tk, p, u10, h, vb, mw, sc_w, ustar, cd_m)
    real, intent(in) :: tk !< temperature at surface in kelvin
    real, intent(in) :: p !< pressure at surface in pa
    real, intent(in) :: u10 !< wind speed at 10m above the surface in m/s
    real, intent(in) :: h !< Henry's law constant (\f$H=c_sg/C_sl\f$) (unitless)
    real, intent(in) :: vb !< Molar volume
    real, intent(in) :: mw !< molecular weight (g/mol)
    real, intent(in) :: sc_w
    real, intent(in), optional :: ustar !< Friction velocity (m/s).  If not provided,
                                        !! ustar = \f$u_{10} \sqrt{C_D}\f$.
    real, intent(in), optional :: cd_m !< Drag coefficient (\f$C_D\f$).  Used only if
                                       !! ustar is provided.
                                       !! If ustar is not provided,
                                       !! cd_m = \f$6.1 \times 10^{-4} + 0.63 \times 10^{-4} *u_10\f$

    real :: ra,rl,tc

    tc = tk-273.15
    ra = 1./max(h*calc_ka(tc,p,mw,vb,u10,ustar,cd_m),epsln)
    rl = 1./max(calc_kl(tc,u10,sc_w),epsln)
    calc_kw = 1./max(ra+rl,epsln)
  end function calc_kw

  !> Calculate \f$k_a\f$
  !!
  !! See calc_kw
  real function calc_ka(t, p, mw, vb, u10, ustar, cd_m)
    real, intent(in) :: t !< temperature at surface in C
    real, intent(in) :: p !< pressure at surface in pa
    real, intent(in) :: mw !< molecular weight (g/mol)
    real, intent(in) :: vb !< molar volume
    real, intent(in) :: u10 !< wind speed at 10m above the surface in m/s
    real, intent(in), optional :: ustar !< Friction velocity (m/s).  If not provided,
                                        !! ustar = \f$u_{10} \sqrt{C_D}\f$.
    real, intent(in), optional :: cd_m !< Drag coefficient (\f$C_D\f$).  Used only if
                                       !! ustar is provided.
                                       !! If ustar is not provided,
                                       !! cd_m = \f$6.1 \times 10^{-4} + 0.63 \times 10^{-4} *u_10\f$

    real             :: sc
    real             :: ustar_t, cd_m_t

    if (.not. present(ustar)) then
      !drag coefficient
      cd_m_t = 6.1e-4 +0.63e-4*u10
      !friction velocity
      ustar_t = u10*sqrt(cd_m_t)
    else
      cd_m_t = cd_m
      ustar_t = ustar
    end if
    sc = schmidt_g(t,p,mw,vb)
    calc_ka = 1e-3+ustar_t/(13.3*sqrt(sc)+1/sqrt(cd_m_t)-5.+log(sc)/(2.*vonkarm))
  end function calc_ka

  !> Calculate \f$k_l\f$
  !!
  !! See calc_kw, and Nightingale, Global Biogeochemical Cycles, 2000
  !! (https://doi.org/10.1029/1999GB900091)
  real function calc_kl(t, v, sc)
    real, intent(in) :: t !< temperature at surface in C
    real, intent(in) :: v !< wind speed at surface in m/s
    real, intent(in) :: sc

    calc_kl = (((0.222*v**2)+0.333*v)*(max(sc,epsln)/600.)**(-0.5))/(100.*3600.)
  end function calc_kl

  !> Schmidt number of the gas in air
  real function schmidt_g(t, p, mw, vb)
    real, intent(in) :: t !< temperature at surface in C
    real, intent(in) :: p !< pressure at surface in pa
    real, intent(in) :: mw !< molecular weight (g/mol)
    real, intent(in) :: vb !< molar volume

    real :: d,v

    d = d_air(t,p,mw,vb)
    v = v_air(t)
    schmidt_g = v / d
  end function schmidt_g

  !> From Fuller, Industrial & Engineering Chemistry (https://doi.org/10.1021/ie50677a007)
  real function d_air(t, p, mw, vb)
    real, intent(in) :: t  !< temperature in c
    real, intent(in) :: p  !< pressure in pa
    real, intent(in) :: mw !< molecular weight (g/mol)
    real, intent(in) :: vb !< diffusion coefficient (\f$cm3/mol\f$)

    real, parameter :: ma = 28.97d0 !< molecular weight air in g/mol
    real, parameter :: va = 20.1d0  !< diffusion volume for air (\f$cm^3/mol\f$)

    real            :: pa

    ! convert p to atm
    pa = 9.8692d-6*p
    d_air = 1d-3 *&
        & (t+273.15d0)**(1.75d0)*sqrt(1d0/ma + 1d0/mw)/(pa*(va**(1d0/3d0)+vb**(1d0/3d0))**2d0)
    ! d_air is in cm2/s convert to m2/s
    d_air = d_air * 1d-4
  end function d_air

  !> kinematic viscosity in air
  real function p_air(t)
    real, intent(in) :: t

    real, parameter :: sd_0 = 1.293393662d0,&
        & sd_1 = -5.538444326d-3,&
        & sd_2 = 3.860201577d-5,&
        & sd_3 = -5.2536065d-7
    p_air = sd_0+(sd_1*t)+(sd_2*t**2)+(sd_3*t**3)
  end function p_air

  !> Kinematic viscosity in air (\f$m^2/s\f$
  real function v_air(t)
    real, intent(in) :: t !< temperature in C
    v_air = n_air(t)/p_air(t)
  end function v_air

  !> dynamic viscosity in air
  real function n_air(t)
    real, intent(in) :: t !< temperature in C

    real, parameter :: sv_0 = 1.715747771d-5,&
        & sv_1 = 4.722402075d-8,&
        & sv_2 = -3.663027156d-10,&
        & sv_3 = 1.873236686d-12,&
        & sv_4 = -8.050218737d-14
    ! in n.s/m^2 (pa.s)
    n_air = sv_0+(sv_1*t)+(sv_2*t**2)+(sv_3*t**3)+(sv_4*t**4)
  end function n_air



  !> A basic linear interpolation.  Is there a tool in FMS already?
  real function interp_linear(XI,YI,XO)
    real, intent(in), dimension(:) :: XI, YI
    real, intent(in) :: XO
    integer :: I,N
    N=size(XI)
    if (XO<XI(1)) then
       interp_linear = 0.
       return
    elseif (XO>XI(N)) then
       interp_linear = YI(N)
       return
    else
       do i=2,N
          if (XI(i)>XO) then
             interp_linear =
YI(i-1)+(XO-XI(i-1))*(YI(i)-YI(i-1))/(XI(i)-XI(i-1))
             return
          endif
       enddo
    endif
    interp_linear = -999.
  end function interp_linear

  !> DM18 mean KW in m/s
  real function DM18_mean(U10)
    real, intent(in) :: U10 !< wind speed in m/s
    !/ Defining a look-up table.
    real, parameter, dimension(42) :: &
         U10table = (/&
          0.,            0.64119104,   1.55962283,   2.53868778,   3.52982184,&
          4.51815114,    5.51201397,   6.50414836,   7.48559029,   8.48238844,&
          9.47579109,   10.48318341,  11.47809756,  12.47843195,  13.47547521,&
         14.47356684,   15.47276396,  16.47617237,  17.47051939,  18.46189453,&
         19.46144313,   20.45343019,  21.45670896,  22.45650745,  23.45583207,&
         24.45696779,   25.45506134,  26.4551677,   27.45797992,  28.45880982,&
         29.46211083,   30.45823377,  31.45787607,  32.45438511,  33.45634555,&
         34.45785125,   35.46573798,  36.45633558,  37.47961353,  38.45226419,&
         39.44741611,  100./)
    real, parameter, dimension(42) :: &
         KWtable = (/&
         0.00000000e+00, 8.33780697e-01, 2.01223280e+00, 3.26239166e+00,&
         4.75191128e+00, 6.60394439e+00, 8.82130898e+00, 1.13624181e+01,&
         1.42113667e+01, 1.75314036e+01, 2.13490067e+01, 2.58402723e+01,&
         3.07720960e+01, 3.63315584e+01, 4.25192163e+01, 4.95036107e+01,&
         5.73993068e+01, 6.61515303e+01, 7.60878357e+01, 8.63987500e+01,&
         9.68969127e+01, 1.08379053e+02, 1.20902745e+02, 1.34376992e+02,&
         1.49098682e+02, 1.63985525e+02, 1.78576016e+02, 1.91774283e+02,&
         2.05455138e+02, 2.20407011e+02, 2.34877452e+02, 2.50304792e+02,&
         2.66092352e+02, 2.83209977e+02, 3.00138583e+02, 3.16192852e+02,&
         3.36680945e+02, 3.55027786e+02, 3.76986764e+02, 3.98538587e+02,&
         4.22731874e+02, 1.47210294e+03/)
    real, parameter :: cph_to_ms = 1./(3600.*100.)
    !/ Linearly interpolate over the look-up table
    DM18_mean = interp_linear(U10table,KWtable,U10) * cph_to_ms
  end function DM18_mean

  !> DM18 KW in m/s
  real function DM18(Ustar,HS,SC,Alpha,SST)
    real, intent(in) :: Ustar, &!< wind speed in m s-1, ustar_air
                        Hs, &   !< wave height in m
                        SC, &   !< Dimensionless Schmidt number
                        Alpha, &!< Solubility in mol m-3 atm-1;
                        SST     !< SST in K
! Note that the model output is ustar_water,
! Ustar_air^2*rho_air=Ustar_water^2*rho_water.
! rho_air=1.225 kg/m3; rho_water=1035 kg/m3
    real, parameter :: c_nb = 1.66e-4 !< Non-bubble flux coefficient
(dimensionless)
    real, parameter :: c_b = 1.1e-5 !< Bubble flux coefficient (s2 m-2)
    real, parameter :: SC_CO2=660 !< Schmidt number for CO2 at standard (?)
conditions
    real, parameter :: gamma = 5./3 !<exponent on ustar (DM18)
    real, parameter :: zeta = 2./3  !< exponent on g*Hs (DM18)
    real, parameter :: atm2pa = 1./101325. !< pressure conversion factor
    real, parameter :: R = 8.314 !< Ideal gas constant, units: k*mol/m3/Pa (k is
kelvin)
    real, parameter :: grav=9.81 !<Gravity (m s-2) should be defined elsewhere?
    DM18 = c_nb*Ustar +&
           c_b*sqrt(SC/SC_CO2)*(Ustar**gamma)*((grav*Hs)**zeta)/&
           (Alpha*atm2pa*R*SST)
    !DM18 = 1e-5*Alpha*atm2pa*R*SST
  end function DM18
end module atmos_ocean_fluxes_calc_mod
