#' variable define
#' @name all_vari
#' @description In this topic all of the variable in the EDCHM will be defined.
#'
#'
#' @param n_time (int) number of time steps
#' @param n_spat (int) number of spatial units
#' @param Time_dayOfYear_,Time_dayOfYear_d <1, 366> (d) the number of the day in the year between 1 (1 January) and 365 or 366 (31 December)
#' @param ATMOS_precipitation_mm (mm/m2/TS) precipitaion volum
#' @param ATMOS_rain_mm (mm/m2/TS) precipitation in rain form
#' @param ATMOS_SNOW_mm (mm/m2/TS) precipitation in snow form
#' @param ATMOS_temperature_Cel (Cel) the average air temperature in the time phase
#' @param ATMOS_temperatureMax_Cel (Cel) the maximal air temperature in the time phase
#' @param ATMOS_temperatureMin_Cel (Cel) the minimal air temperature in the time phase
#' @param ATMOS_solarRadiat_MJ (MJ/m2/TS) the solar radiation that actually reaches the earths surface
#' @param ATMOS_netRadiat_MJ (MJ/m2/TS) the balance between the energy absorbed, reflected and emitted by the earths surface or the difference between the incoming net shortwave (Rns) and the net outgoing longwave (Rnl) radiation
#' @param ATMOS_relativeHumidity_1 (0, 1) relative humidity
#' @param ATMOS_vaporPress_hPa (hPa) actual vapour pressure, can be calculated by [ATMOS_VaporPress()]
#' @param ATMOS_saturatVaporPress_hPa (hPa) saturation vapour pressure at `ATMOS_temperature_Cel`, can be calculated by [ATMOS_SaturatVaporPress()]
#' @param ATMOS_windSpeed_m_s (m/s) measured wind speed at z m above ground surface
#' @param ATMOS_windMeasureHeight_m (m) height of measurement above ground surface
#' @param ATMOS_windSpeed2m_m_s (m/s) wind speed at 2 m above ground surface
#' @param ATMOS_potentialEvatrans_mm (mm/m2/TS) **potential / reference** evapotranspiration
#' @param LAND_albedo_1 <0, 1> albedo of the region
#' @param LAND_growUpDay_d growup day
#' @param LAND_latitude_Degree,CELL_latitude_deg (degree) average latitude
#' @param LAND_elevation_m (m) average elevation
#' @param LAND_impermeableFrac_1 <0, 1> the maximum impermeable fraction when th soil is fully saturated
#' @param LAND_water_mm (mm/m2) **pounded water** volume in `landLy` and there is no limit, different than `LAND_interceptWater_mm`
#' @param LAND_interceptWater_mm (mm/m2) initial water volume that can be **intercepted**
#' @param LAND_interceptCapacity_mm (mm/m2) average intercept Capacity (maximal storage capacity)
#' @param LAND_actualEvatrans_mm (mm/m2/TS) **actual** evapotranspiration from `landLy`
#' @param LAND_infiltrat_mm (mm/m2/TS) infiltration from `landLy` to `soilLy`
#' @param LAND_leafAreaIndexMin_,LAND_leafAreaIndexMax_ minimal and maximal LAI
#' @param SNOW_ice_mm (mm/m2) water equivalent of **ice** in snowpack
#' @param SOIL_fieldCapacityPerc_1 <0, 1> the relative ratio that the water content can drainage by gravity
#' @param SOIL_water_mm (mm/m2) water volume in `soilLy`
#' @param SOIL_capacity_mm (mm/m2) average soil Capacity (maximal storage capacity)
#' @param SOIL_interflow_mm (mm/m2/TS) subsurface flow directly to the river
#' @param SOIL_actualEvatrans_mm (mm/m2/TS) **actual** evapotranspiration from `soilLy`
#' @param SOIL_potentialPercola_mm <0.01, 7> (mm/m2/TS) **potential** percolation
#' @param SOIL_potentialInteflow_mm <0.01, 7> (mm/m2/TS) **potential** interflow
#' @param SOIL_potentialCapirise_mm <0.01, 7> (mm/m2/TS) **potential** capillary rise
#' @param GROUND_water_mm (mm/m2/TS) water volume in `groundLy`
#' @param GROUND_capacity_mm (mm/m2) water storage capacity in `groundLy`
#' @param GROUND_lateral_mm (mm/m2/TS) lateral flow, exchange with outside region. It can be **NEGATIV**
#' @param GROUND_capillarise_mm (mm/m2/TS) capillary rise from `groundLy` to `soilLy`
#' @param GROUND_potentialLateral_mm <-7, 7> (mm/m2/TS) **potential** lateral flow
#' @param GROUND_potentialBaseflow_mm <0.01, 7> (mm/m2/TS) **potential** baseflow
#' @param CONFLUEN_cellInflow_m3 (m3 / TS) input water volum in every routeline
#' @param CONFLUEN_inputWater_mm,LAND_runoff_mm,GROUND_baseflow_mm (mm/m2) input water volum in every routeline
#' @param CONFLUEN_iuh_1,CONFLUEN_iuhLand_1,CONFLUEN_iuhSoil_1,CONFLUEN_iuhGround_1 (vector of num, sume() = 1) the ratio in every timestep, can be calculated by [confluenIUH_GR4J1()], [confluenIUH_GR4J2()]
#' @param CONFLUEN_responseTime_TS (TS) response or concentration time in every routeline
#' @param water_mm (mm/m2/TS) water volume in `soilLy`, `groundLy` or intercept of `landLy`, same as `SOIL_water_mm`, `GROUND_water_mm` or `LAND_interceptWater_mm`
#' @param capacity_mm (mm/m2) water storage capacity in `soilLy`, `groundLy` or intercept of `landLy`, same as `SOIL_capacity_mm`, `GROUND_capacity_mm` or`LAND_interceptCapacity_mm`
#' @param Lake_cellNumber_int (int) lake cell number
#' @param Lake_water_m3 (m3) lake water volume
#' @param Lake_EVATRANS_mm (mm) evapotranspiration from the lake
#' @param Lake_capacity_m3 (m3) capacity of the lake
#' @param Lake_area_km2 (km2) area of the lake
#' @param Lake_inflow_m3 (m3/TS) inflow of the lake
#' @param RIVER_water_m3 (m3) river water volume
#' @param RIVER_inflow_m3 (m3/TS) river inflow volume
#' @param RIVER_outflow_m3 (m3/TS) river outflow volume
#' @param RIVER_velocity_km (km/TS) river velocity
#' @param RIVER_length_km (km) river length
#' @param Riverlak_cellNumber_int (int) river-lake cell number
#' @param Riverlak_area_km2 (km2) area of the river-lake system
#' @param Riverlak_capacity_m3 (m3) total capacity of the river-lake system
#' @param Riverlak_inflow_m3 (m3) inflow volume to the river-lake system
#' @param Riverlak_water_m3 (m3) water volume in the river-lake system
#' @param RIVER_upstreamInflow_m3 (m3) Upstream flow
#' @param Reservoi_cellNumber_int (int) reservoir cell number
#' @param Reservoi_water_m3 (m3) volume of water in the reservoir
#' @param Reservoi_capacity_m3 (m3) reservoir capacity
#' @param Reservoi_EVATRANS_mm (mm) evapotranspiration from the reservoir
#' @param Reservoi_inflow_m3 (m3/TS) inflow to the reservoir
#' @param Reservoi_demand_m3 (m3/TS) water demand from the reservoir
#' @param Reservoi_potentialVolume_m3 (m3) potential capacity of the reservoir
#' @param Reservoi_area_km2 (km2) surface area of the reservoir
#' @param Reservoi_yearInflow_m3 (m3/TS) annual mean inflow to the reservoir
#' @param Reservoi_yearDemand_m3 (m3/TS) annual mean water demand from the reservoir
#' @param Reservoi_yearRelase_m3 (m3/TS) reference to annual mean water release from the reservoir
#' @param Reservoi_isOperateStart_01 (01) indicates if reservoir operation has started (0 or 1)
#' @param Reservoi_isIrrigate_01 (01) indicates if the reservoir is used for irrigation (0 or 1)
#' @param CELL_landArea_km2 (km2) land area of the basin
#' @param CELL_cellNumberStep_int (int) step number of basin cells
#' @param CELL_inflowCellNumberStep_int (int) step number of inflow cells in the basin
#' @param Upstream_cellNumber_int Upstream cell number
#' @param Upstream_streamflow_m3 Upstream
#' @param if_allVariExport (bool) whether export all the state variables (default in false)
NULL



