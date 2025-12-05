# Changelog

## 0.16 ()

## 0.15 (2024-12-06)

### Changed

- Use January 1st, `annee_ref` as reference date for variable `temps`
  of `restart.nc` (851bb5fc)
- Do not write `iday_end`, iim, jjm, llm, `annee_ref`, daysec, ra,
  preff, rg, rcpd, rkappa, romega into variable `controle` in
  `restart.nc` (731a472a, 5cb2b084, c03b691d, 70ae788b)
- Use global attribute `itau_dyn` rather than variable `controle` in
  `restart.nc` (6ae8ef8b)
- Write 0 at position 24 into variable `controle` in `restart.nc`
  (9ccb8262)

## 0.14 (2024-11-29)

### Added

- Add documentation (ab764e43)

### Removed

- Do not write radsol to restartphy (5d056986)
- Remove useless print (1add0130)

## 0.13 (2024-02-09)

### Changed

- Rename `find_coord` to `nf95_find_coord` (4b779da0)
- Rename NetCDF variables `rain_f` and `snow_f` in `(re)startphy.nc`
  to `rain_fall` and `snow_fall` (dcfe4a2b)
- Set time origin in XIOS from `annee_ref` and `day_ref` (71e0f7d4)
- Set start date in XIOS from `annee_ref` and `day_ini` (06a6364d)
- Stop for bad percentages (8245c594)
- Rename directory `Test_input` to `Tests_input` (0751e983)

### Added

- Check consistency of `iflag_cldcon` and `conv_emanuel` (b732e3e6)
- Send a field to XIOS (c3e18740, f787859b)

### Removed

- Remove possibility to choose calendar (baccbc2f)
- Remove XIOS and MPI from program `test_soil` (645f0de5)

### Fixed

- Correct calendar for XIOS (0889fed3)
- Correct time step for XIOS (8dba6792, c5faa56b)
- Use list-directed formatting in `phylmd/suphec.f90`, removing
  compiler warnings (366093c9)
- Do not assume names for longitude and latitude in `landiceref.nc`
  (9a48a41a)
- Manage temperature in Celsius degrees in `SST.nc`(5fa42d52)
- Correct pctsrf where masque is exactly equal to epsfra (8d824548)

## 0.12 (2022-11-16)

### Changed

- Replace cmake directory by submodule (1900b796, 50001d9b, d36c534f)
- Replace `NR_util` by Jumble (16c9321f)
- Make the final call to caldyn with conser true (therefore adding to
  standard output) (2734d49b)
- Update to Doxygen format `1.9` (d4bc9ca6)
- Replace attribute title by `long_name` in `restartphy.nc` (0042054d)
- Rename NetCDF variable SNOW in `restartphy.nc` to fsnow (367a5778)
- Rename NetCDF variable QS in `restartphy.nc` to fqsurf (8f1cfa4e)
- Rename NetCDF dimension `points_physiques` in `restartphy.nc` to
  klon (19802207)

### Added

- Add Doxygen configuration file (79bec3cf)

### Fixed

- Find NetCDF-Fortran even if not FETCH (fce9560a)
- Define ucont entirely in `dyn3d/covcont.f90` (c84735a4)

### Removed

- Remove sub-dimensioning of radiative transfer to kdlon (5cc969b8)
- Remove computation of aerosol direct effect (so `ok_ade` is removed
  from namelist `physiq_nml`) (0c568711)
- Remove test on convection in `phylmd/phytrac.f90` (9420b7b9)

## 0.11 (2021-12-08)

### Added

- Add option to fetch libraries (98360618, c9d79a3a)

### Fixed

- Add blank (1eb4d7cd)

### Removed

- Leave compiler options to the user (1e7bd9f3)
- Do not find NetCDF (d985ebc9)

### Changed

- Update cmake find modules (19b44a63, 6291fdd4)
- Preprend jumble, `nr_util` and `numer_rec_95` with namespaces
  (c10dfe04, f439722e)

## 0.10 (2021-09-15)

### Removed

- Remove call to `transp_lay` and do not write `ue_lay`, `ve_lay`,
  `uq_lay`, `vq_lay` to `histins.nc` (aa846d9c)

### Changed

- Replace moist by dry static energy transport in variables ue and ve
  of `histins.nc` (77763b16)
- Update compiler options and cmake find modules (c6805a6e)

## 0.9 (2021-05-27)

### Changed

- Collect commands defining `XIOS::XIOS` into a new module
  `FindXIOS.cmake` (327bee46)
- Remove constraint on dimensions of `Relief.nc` (ba252003)
- Transfer variables `cld_lc_lsc`, `cld_lc_con`, `cld_tau_lsc`,
  `cld_tau_con`, `ffallv_lsc`, `ffallv_con`, `coef_eva`, `reevap_ice`,
  `iflag_pdf` of namelist `conf_phys_nml` to new namelist
  `fisrtilp_nml` (152c5d21)

### Fixed

- Do not print `d_q_the` before it is defined in
  `phylmd/Thermcell/calltherm.f90` (22e80ae1)
- Add missing dependencies in CMake files (07e69c50)

### Removed

- Remove option `nsplit_thermals` (1a44b22a)

### Added

- Add output of variables ue, ve, uq, vq, `ue_lay`, `ve_lay`,
  `uq_lay`, `vq_lay` in `histins.nc` (a975ed77)

## 0.8 (2021-04-28)

### Changed

- Move `.cmake` files to new cmake directory (31619bd0)
- Set property `IMPORTED_LINK_INTERFACE_LANGUAGES` instead of naming
  library file `stdc++` in `CMakeLists.txt`

## 0.7 (2021-03-29)

### Changed

- Reference NetCDF95 namespace (f0b4156a)
- Define imported target XIOS in `CMakeLists.txt` (e0fcaae1)
- Call MPI, just initializing and finalizing, in programs `gcm` and
  `test_soil` (462106d9, 8fa9373d)
- Define context for XIOS in programs `gcm` and `test_soil` (af82b745,
  8fa9373d)
- Modify initial rugosity of land and land-ice in `startphy.nc` (ab9bee91)

### Removed

- Do not write dtvr to  an element of `tab_cntrl` in `restart.nc` (9e98c5e3)

## 0.6 (2021-03-19)

### Changed

- Create directory `dynphy_lonlat` and move files into it
- Change calendar name from 360d to `360_day` in history files
- Replace 3.14 by pi in `phylmd/Interface_surf/soil.f90`
- Read `inertie_*` in namelist `soil_nml`
- Accelerate `phylmd/Interface_surf/soil.f90` (a022d0f7)
- Read $\lambda C$ instead of $\sqrt{\lambda C / (1 s)}$ in namelist
  `soil_nml` (6c64d74b)
- Read ngroup in namelist `logic_nml` and check consistency with iim
  and jjm (a283ab79, 2a7d63f1)
- Use module `FindNetCDF_Fortran` in `CMakeLists.txt`

### Fixed

- Correct intent in `phylmd/Interface_surf/climb_hq_up.f90`
- Do not divide by pkf in `phylmd/Interface_surf/climb_hq_down.f90`
- Correct intent in `phylmd/Interface_surf/soil.f90`
- Require XIOS at configuration time (730ebfb4)
- Fix creation of TAGS by CMake (f57bdc62, 6e5e4a1f, d6a57715)

### Removed

- Do not print range of `q_sat` in `dynphy_lonlat/etat0.f90`

### Added

- Add main program `test_soil`
- Add output from `phylmd/Interface_surf/soil.f90` (9eb5b0a4)

## 0.5 (2020-08-10)

### Changed

- Read "sea_ice.nc" and "SST.nc" instead of
  `amipbc_sic_360x180_1979_2008_clim.nc` and
  `amipbc_sst_360x180_1979_2008_clim.nc`
- Accelerate procedure coefkz2
- Read from `soil_nml` instead of `soil.def`
- Rename NetCDF variable TEMPS to time in `limit.nc`

### Fixed

- Circumvent gfortran bug 96436 with `f0.d`
- Define coefm and coefh in `phylmd/Interface_surf/coefkz2.f90`
- Recompute zc and zd from previous time step in
  `phylmd/Interface_surf/soil.f90`
- Initialize therm in `phylmd/Interface_surf/hbtm.f90`
- Initialize fqsurf in `dyn3d/etat0.f90`

### Removed

- Remove printing of fz from `phylmd/Interface_surf/soil.f90`

### Added

- Add longitude and latitude to file `limit.nc`
- Add variable `ffonte` to `histins.nc`

## 0.4 (2020-07-15)

### Changed

- Move directory Namelists and file `tests_LMDZE.json` to directory
  `Test_input`
- Use `find_coord` instead of requiring coordinate names in
  `amipbc_sic_1x1.nc`, `amipbc_sst_1x1.nc` and `'ECDYN.nc`
  
### Fixed

- Fix missing library in `CMakeLists.txt`
- Fix missing initializations in `phylmd/Thermcell/thermcell.f90`

## 0.3 (2020-07-14)

### Changed

- Read spatial dimensions iim, jjm, llm at run time, in namelist
  `dimensions_nml`
- Use pkgconfig to find NetCDF-Fortran

### Added

- Import directory Namelists into the repository
- Import `Test_input` directory

### Removed

- Remove choice mixq in `phylmd/ajsec.f90`
