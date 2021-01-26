import os
from osgeo import gdal
import numpy as np
import warnings

def main(input_folder, output_folder, Date):

    # Do not show warnings
    warnings.filterwarnings('ignore')

    import pyWAPOR.ETLook as ETLook
    import pyWAPOR.Functions.Processing_Functions as PF
    import pyWAPOR.ETLook.outputs as out

    # Define Date string
    Date_str = "%d%02d%02d" %(Date.year, Date.month, Date.day)

    # Input folder Date
    input_folder_date = os.path.join(input_folder, Date_str)

    ############################ Define inputs ################################

    #input_files
    ALBEDO_filename = os.path.join(input_folder_date, "ALBEDO_%s.tif" %Date_str)
    NDVI_filename = os.path.join(input_folder_date, "NDVI_%s.tif" %Date_str)
    LST_filename = os.path.join(input_folder_date, "LST_%s.tif" %Date_str)
    LAI_filename = os.path.join(input_folder_date, "LAI_%s.tif" %Date_str)
    FVC_filename = os.path.join(input_folder_date, "FVC_%s.tif" %Date_str)
    Time_filename = os.path.join(input_folder_date, "Time_%s.tif" %Date_str)
    Lat_filename = os.path.join(input_folder_date, "Lat_%s.tif" %Date_str)
    Lon_filename = os.path.join(input_folder_date, "Lon_%s.tif" %Date_str)
    DEM_filename = os.path.join(input_folder_date, "DEM.tif")
    Slope_filename = os.path.join(input_folder_date, "Slope.tif")
    Aspect_filename = os.path.join(input_folder_date, "Aspect.tif")
    LandMask_filename = os.path.join(input_folder_date, "LandMask.tif")
    Bulk_filename =os.path.join(input_folder_date, "Bulk_Stomatal_resistance.tif")
    MaxObs_filename = os.path.join(input_folder_date, "Maximum_Obstacle_Height.tif")
    Pair_24_0_filename = os.path.join(input_folder_date, "Pair_24_0_%s.tif" %Date_str)
    Pair_inst_0_filename = os.path.join(input_folder_date, "Pair_inst_0_%s.tif" %Date_str)
    Pair_inst_filename = os.path.join(input_folder_date, "Pair_inst_%s.tif" %Date_str)
    Pre_filename = os.path.join(input_folder_date, "Precipitation_%s.tif" %Date_str)
    Hum_24_filename = os.path.join(input_folder_date, "qv_24_%s.tif" %Date_str)
    Hum_inst_filename = os.path.join(input_folder_date, "qv_inst_%s.tif" %Date_str)
    Tair_24_filename = os.path.join(input_folder_date, "tair_24_%s.tif" %Date_str)
    Tair_max_24_filename = os.path.join(input_folder_date, "tair_max_24_%s.tif" %Date_str)
    Tair_min_24_filename = os.path.join(input_folder_date, "tair_min_24_%s.tif" %Date_str)
    Tair_inst_filename = os.path.join(input_folder_date,"tair_inst_%s.tif" %Date_str)
    Tair_amp_filename = os.path.join(input_folder_date, "Tair_amp_%s.tif" %Date_str)
    Wind_24_filename = os.path.join(input_folder_date, "wind_24_%s.tif" %Date_str)
    Wind_inst_filename = os.path.join(input_folder_date, "wind_inst_%s.tif" %Date_str)
    WatCol_inst_filename = os.path.join(input_folder_date, "wv_inst_%s.tif" %Date_str)
    Trans_24_filename = os.path.join(input_folder_date, "Trans_24_%s.tif" %Date_str)

    ############################ Define outputs ###############################

    # Output folder Date
    output_folder_date = os.path.join(output_folder, Date_str)
    if not os.path.exists(output_folder_date):
        os.makedirs(output_folder_date)

    #output_files
    vc_filename = os.path.join(output_folder_date, "vc_%s.tif" %Date_str)
    lai_filename = os.path.join(output_folder_date, "LAI_%s.tif" %Date_str)
    lai_eff_filename= os.path.join(output_folder_date, "LAI_eff_%s.tif" %Date_str)
    sf_soil_filename = os.path.join(output_folder_date, "sf_soil_%s.tif" %Date_str)
    lat_filename= os.path.join(output_folder_date, "lat_%s.tif" %Date_str)
    slope_filename= os.path.join(output_folder_date, "slope_%s.tif" %Date_str)
    aspect_filename = os.path.join(output_folder_date, "aspect_%s.tif" %Date_str)
    ra_24_toa_filename = os.path.join(output_folder_date, "ra_24_toa_%s.tif" %Date_str)
    ws_filename = os.path.join(output_folder_date, "ws_%s.tif" %Date_str)
    diffusion_index_filename = os.path.join(output_folder_date, "diffusion_index_%s.tif" %Date_str)
    ra_24_filename  = os.path.join(output_folder_date, "ra_24_%s.tif" %Date_str)
    stress_rad_filename = os.path.join(output_folder_date, "stress_rad_%s.tif" %Date_str)
    p_air_24_filename = os.path.join(output_folder_date, "p_air_24_%s.tif" %Date_str)
    vp_24_filename = os.path.join(output_folder_date, "vp_24_%s.tif" %Date_str)
    svp_24_filename = os.path.join(output_folder_date, "svp_24_%s.tif" %Date_str)
    vpd_24_filename = os.path.join(output_folder_date, "vpd_24_%s.tif" %Date_str)
    stress_vpd_filename  = os.path.join(output_folder_date, "stress_vpd_%s.tif" %Date_str)
    stress_temp_filename = os.path.join(output_folder_date, "stress_temp_%s.tif" %Date_str)
    r_canopy_0_filename= os.path.join(output_folder_date, "r_canopy_0_%s.tif" %Date_str)
    t_air_k_24_filename = os.path.join(output_folder_date, "t_air_k_24_%s.tif" %Date_str)
    l_net_filename = os.path.join(output_folder_date, "l_net_%s.tif" %Date_str)
    int_mm_filename = os.path.join(output_folder_date, "int_mm_%s.tif" %Date_str)
    lh_24_filename = os.path.join(output_folder_date, "lh_24_%s.tif" %Date_str)
    int_wm2_filename  = os.path.join(output_folder_date, "int_wm2_%s.tif" %Date_str)
    rn_24_filename = os.path.join(output_folder_date, "rn_24_%s.tif" %Date_str)
    rn_24_canopy_filename= os.path.join(output_folder_date, "rn_24_canopy_%s.tif" %Date_str)
    t_air_k_i_filename = os.path.join(output_folder_date, "t_air_k_i_%s.tif" %Date_str)
    vp_i_filename = os.path.join(output_folder_date, "vp_i_%s.tif" %Date_str)
    ad_moist_i_filename= os.path.join(output_folder_date, "ad_moist_i_%s.tif" %Date_str)
    ad_dry_i_filename = os.path.join(output_folder_date, "ad_dry_i_%s.tif" %Date_str)
    ad_i_filename= os.path.join(output_folder_date, "ad_i_%s.tif" %Date_str)
    u_b_i_bare_filename= os.path.join(output_folder_date, "u_b_i_bare_%s.tif" %Date_str)
    lon_filename= os.path.join(output_folder_date, "lon_%s.tif" %Date_str)
    ha_filename= os.path.join(output_folder_date, "ha_%s.tif" %Date_str)
    ied_filename= os.path.join(output_folder_date, "ied_%s.tif" %Date_str)
    h0_filename = os.path.join(output_folder_date, "h0_%s.tif" %Date_str)
    h0ref_filename = os.path.join(output_folder_date, "h0ref_%s.tif" %Date_str)
    m_filename = os.path.join(output_folder_date, "m_%s.tif" %Date_str)
    rotm_filename = os.path.join(output_folder_date, "rotm_%s.tif" %Date_str)
    Tl2_filename  = os.path.join(output_folder_date, "Tl2_%s.tif" %Date_str)
    B0c_filename = os.path.join(output_folder_date, "B0c_%s.tif" %Date_str)
    Bhc_filename = os.path.join(output_folder_date, "Bhc_%s.tif" %Date_str)
    Dhc_filename = os.path.join(output_folder_date, "Dhc_%s.tif" %Date_str)
    ra_hor_clear_i_filename = os.path.join(output_folder_date, "ra_hor_clear_i_%s.tif" %Date_str)
    emiss_atm_i_filename  = os.path.join(output_folder_date, "emiss_atm_i_%s.tif" %Date_str)
    rn_bare_filename = os.path.join(output_folder_date, "rn_bare_%s.tif" %Date_str)
    rn_full_filename= os.path.join(output_folder_date, "rn_full_%s.tif" %Date_str)
    u_b_i_full_filename = os.path.join(output_folder_date, "u_b_i_full_%s.tif" %Date_str)
    u_star_i_bare_filename = os.path.join(output_folder_date, "u_star_i_bare_%s.tif" %Date_str)
    u_star_i_full_filename = os.path.join(output_folder_date, "u_star_i_full_%s.tif" %Date_str)
    u_i_soil_filename = os.path.join(output_folder_date, "u_i_soil_%s.tif" %Date_str)
    ras_filename  = os.path.join(output_folder_date, "ras_%s.tif" %Date_str)
    raa_filename = os.path.join(output_folder_date, "raa_%s.tif" %Date_str)
    rac_filename= os.path.join(output_folder_date, "rac_%s.tif" %Date_str)
    t_max_bare_filename = os.path.join(output_folder_date, "t_max_bare_%s.tif" %Date_str)
    t_max_full_filename= os.path.join(output_folder_date, "t_max_full_%s.tif" %Date_str)
    w_i_filename = os.path.join(output_folder_date, "w_i_%s.tif" %Date_str)
    t_dew_i_filename = os.path.join(output_folder_date, "t_dew_i_%s.tif" %Date_str)
    t_wet_i_filename = os.path.join(output_folder_date, "t_wet_i_%s.tif" %Date_str)
    t_wet_k_i_filename = os.path.join(output_folder_date, "t_wet_k_i_%s.tif" %Date_str)
    lst_max_filename  = os.path.join(output_folder_date, "lst_max_%s.tif" %Date_str)
    se_root_filename = os.path.join(output_folder_date, "se_root_%s.tif" %Date_str)
    stress_moist_filename= os.path.join(output_folder_date, "stress_moist_%s.tif" %Date_str)
    r_canopy_0_filename= os.path.join(output_folder_date, "r_canopy_0_%s.tif" %Date_str)
    r_canopy_filename= os.path.join(output_folder_date, "r_canopy_%s.tif" %Date_str)
    z_obst_filename = os.path.join(output_folder_date, "z_obst_%s.tif" %Date_str)
    z_oro_filename = os.path.join(output_folder_date, "z_oro_%s.tif" %Date_str)
    z0m_filename = os.path.join(output_folder_date, "z0m_%s.tif" %Date_str)
    ra_canopy_init_filename = os.path.join(output_folder_date, "ra_canopy_init_%s.tif" %Date_str)
    u_b_24_filename = os.path.join(output_folder_date, "u_b_24_%s.tif" %Date_str)
    disp_filename = os.path.join(output_folder_date, "disp_%s.tif" %Date_str)
    u_star_24_init_filename = os.path.join(output_folder_date, "u_star_24_init_%s.tif" %Date_str)
    ad_dry_24_filename = os.path.join(output_folder_date, "ad_dry_24_%s.tif" %Date_str)
    ad_moist_24_filename = os.path.join(output_folder_date, "ad_moist_24_%s.tif" %Date_str)
    ad_24_filename = os.path.join(output_folder_date, "ad_24_%s.tif" %Date_str)
    psy_24_filename = os.path.join(output_folder_date, "psy_24_%s.tif" %Date_str)
    ssvp_24_filename = os.path.join(output_folder_date, "ssvp_24_%s.tif" %Date_str)
    t_24_init_filename = os.path.join(output_folder_date, "t_24_init_%s.tif" %Date_str)
    h_canopy_24_init_filename= os.path.join(output_folder_date, "h_canopy_24_init_%s.tif" %Date_str)
    t_24_filename= os.path.join(output_folder_date, "t_24_%s.tif" %Date_str)
    t_24_mm_filename= os.path.join(output_folder_date, "t_24_mm_%s.tif" %Date_str)
    sf_soil_filename= os.path.join(output_folder_date, "sf_soil_%s.tif" %Date_str)
    rn_24_soil_filename= os.path.join(output_folder_date, "rn_24_soil_%s.tif" %Date_str)
    r_soil_filename= os.path.join(output_folder_date, "r_soil_%s.tif" %Date_str)
    ra_soil_init_filename= os.path.join(output_folder_date, "ra_soil_init_%s.tif" %Date_str)
    u_b_24_filename= os.path.join(output_folder_date, "u_b_24_%s.tif" %Date_str)
    u_star_24_soil_init_filename= os.path.join(output_folder_date, "u_star_24_soil_init_%s.tif" %Date_str)
    g0_bs_filename= os.path.join(output_folder_date, "g0_bs_%s.tif" %Date_str)
    g0_24_filename= os.path.join(output_folder_date, "g0_24_%s.tif" %Date_str)
    e_24_init_filename= os.path.join(output_folder_date, "e_24_init_%s.tif" %Date_str)
    h_soil_24_init_filename= os.path.join(output_folder_date, "h_soil_24_init_%s.tif" %Date_str)
    e_24_filename= os.path.join(output_folder_date, "e_24_%s.tif" %Date_str)
    e_24_mm_filename= os.path.join(output_folder_date, "e_24_mm_%s.tif" %Date_str)
    et_24_mm_filename= os.path.join(output_folder_date, "et_24_mm_%s.tif" %Date_str)
    rn_24_grass_filename= os.path.join(output_folder_date, "rn_24_grass_%s.tif" %Date_str)
    et_ref_24_filename= os.path.join(output_folder_date, "et_ref_24_%s.tif" %Date_str)
    et_ref_24_mm_filename= os.path.join(output_folder_date, "et_ref_24_mm_%s.tif" %Date_str)

    ########################## Open input rasters #############################
    dest_lst = gdal.Open(LST_filename)
    lst = dest_lst.GetRasterBand(1).ReadAsArray()
    lst[lst == -9999] = np.nan

    dest_albedo = gdal.Open(ALBEDO_filename)
    r0 = dest_albedo.GetRasterBand(1).ReadAsArray()
    r0[np.isnan(lst)] = np.nan

    dest_ndvi = gdal.Open(NDVI_filename)
    ndvi = dest_ndvi.GetRasterBand(1).ReadAsArray()
    ndvi[np.isnan(lst)] = np.nan

    dest_lai = gdal.Open(LAI_filename)
    lai = dest_lai.GetRasterBand(1).ReadAsArray()
    lai = np.maximum(0, lai)
    lai[np.isnan(lst)] = np.nan

    dest_fvc = gdal.Open(FVC_filename)
    vc = dest_fvc.GetRasterBand(1).ReadAsArray()
    vc = np.maximum(0, vc)
    vc[np.isnan(lst)] = np.nan

    desttime = gdal.Open(Time_filename)
    dtime = desttime.GetRasterBand(1).ReadAsArray()
    dtime[np.isnan(lst)] = np.nan

    dest_lat = gdal.Open(Lat_filename)
    lat_deg = dest_lat.GetRasterBand(1).ReadAsArray()
    lat_deg[np.isnan(lst)] = np.nan

    dest_lon = gdal.Open(Lon_filename)
    lon_deg = dest_lon.GetRasterBand(1).ReadAsArray()
    lon_deg[np.isnan(lst)] = np.nan

    dest_dem = gdal.Open(DEM_filename)
    z = dest_dem.GetRasterBand(1).ReadAsArray()
    z[np.isnan(lst)] = np.nan

    dest_slope = gdal.Open(Slope_filename)
    slope_deg = dest_slope.GetRasterBand(1).ReadAsArray()
    slope_deg[np.isnan(lst)] = np.nan

    dest_aspect = gdal.Open(Aspect_filename)
    aspect_deg = dest_aspect.GetRasterBand(1).ReadAsArray()
    aspect_deg[np.isnan(lst)] = np.nan

    dest_lm = gdal.Open(LandMask_filename)
    land_mask = dest_lm.GetRasterBand(1).ReadAsArray()
    land_mask[np.isnan(lst)] = np.nan

    #dest_bulk = gdal.Open(Bulk_filename)
    #bulk = dest_bulk.GetRasterBand(1).ReadAsArray()

    dest_maxobs = gdal.Open(MaxObs_filename)
    z_obst_max = dest_maxobs.GetRasterBand(1).ReadAsArray()
    z_obst_max[np.isnan(lst)] = np.nan

    dest_pairsea24 = gdal.Open(Pair_24_0_filename)
    p_air_0_24 = dest_pairsea24.GetRasterBand(1).ReadAsArray()
    p_air_0_24 = ETLook.meteo.air_pressure_kpa2mbar(p_air_0_24)
    p_air_0_24[np.isnan(lst)] = np.nan

    dest_pairseainst = gdal.Open(Pair_inst_0_filename)
    p_air_0_i = dest_pairseainst.GetRasterBand(1).ReadAsArray()
    p_air_0_i = ETLook.meteo.air_pressure_kpa2mbar(p_air_0_i)
    p_air_0_i[np.isnan(lst)] = np.nan

    dest_pairinst = gdal.Open(Pair_inst_filename)
    p_air_i = dest_pairinst.GetRasterBand(1).ReadAsArray()
    p_air_i = ETLook.meteo.air_pressure_kpa2mbar(p_air_i)
    p_air_i[np.isnan(lst)] = np.nan

    dest_precip = gdal.Open(Pre_filename)
    P_24 = dest_precip.GetRasterBand(1).ReadAsArray()
    P_24[np.isnan(lst)] = np.nan

    dest_hum24 = gdal.Open(Hum_24_filename)
    qv_24 = dest_hum24.GetRasterBand(1).ReadAsArray()
    qv_24[np.isnan(lst)] = np.nan

    dest_huminst = gdal.Open(Hum_inst_filename)
    qv_i = dest_huminst.GetRasterBand(1).ReadAsArray()
    qv_i[np.isnan(lst)] = np.nan

    dest_tair24 = gdal.Open(Tair_24_filename)
    t_air_24 = dest_tair24.GetRasterBand(1).ReadAsArray()
    #t_air_24 = ETLook.meteo.disaggregate_air_temperature_daily(t_air_24_coarse, z, z_coarse, lapse)
    t_air_24[np.isnan(lst)] = np.nan

    dest_tair24 = gdal.Open(Tair_max_24_filename)
    t_air_max_24 = dest_tair24.GetRasterBand(1).ReadAsArray()
    t_air_max_24[np.isnan(lst)] = np.nan

    dest_tair24 = gdal.Open(Tair_min_24_filename)
    t_air_min_24 = dest_tair24.GetRasterBand(1).ReadAsArray()
    t_air_min_24[np.isnan(lst)] = np.nan

    dest_tairinst = gdal.Open(Tair_inst_filename)
    t_air_i = dest_tairinst.GetRasterBand(1).ReadAsArray()
    t_air_i[np.isnan(lst)] = np.nan

    dest_tairamp = gdal.Open(Tair_amp_filename)
    t_amp_year = dest_tairamp.GetRasterBand(1).ReadAsArray()
    t_amp_year[np.isnan(lst)] = np.nan

    dest_wind24 = gdal.Open(Wind_24_filename)
    u_24 = dest_wind24.GetRasterBand(1).ReadAsArray()
    u_24[np.isnan(lst)] = np.nan

    dest_windinst = gdal.Open(Wind_inst_filename)
    u_i = dest_windinst.GetRasterBand(1).ReadAsArray()
    u_i[np.isnan(lst)] = np.nan

    dest_watcol = gdal.Open(WatCol_inst_filename)
    wv_i = dest_watcol.GetRasterBand(1).ReadAsArray()
    wv_i[np.isnan(lst)] = np.nan

    dest_trans = gdal.Open(Trans_24_filename)
    trans_24 = dest_trans.GetRasterBand(1).ReadAsArray()
    trans_24[np.isnan(lst)] = np.nan

    # example file
    geo_ex = dest_albedo.GetGeoTransform()
    proj_ex = dest_albedo.GetProjection()

    ########################## Open input constants ###########################

    doy = int(Date.strftime("%j"))
    aod550_i = 0.01 # https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/61/MOD04_L2 heb niet echt een standaard product hiervan gevonden
    se_top = 0.5
    porosity = 0.4

    '''
    http://lawr.ucdavis.edu/classes/SSC100/probsets/pset01.html
    6. Calculate the porosity of a soil sample that has a bulk density of 1.35 g/cm3. Assume the particle density is 2.65 g/cm3.
    Porosity = (1-(r b/r d) x 100 = (1-(1.35/2.65)) x 100 = 49%
    '''

    # Create QC array
    QC = np.ones(lst.shape)
    QC[np.isnan(lst)] = np.nan

    # page 31 flow diagram

    # **effective_leaf_area_index**************************************************
    # constants or predefined:
    nd_min = 0.125
    nd_max = 0.8
    vc_pow = 0.7
    vc_min = 0
    vc_max = 0.9677324224821418
    lai_pow = -0.45

    # **atmospheric canopy resistance***********************************************
    # constants or predefined:
    diffusion_slope = -1.33
    diffusion_intercept = 1.15
    t_opt = 25 # optimal temperature for plant growth
    t_min = 0 # minimal temperature for plant growth
    t_max = 50 # maximal temperature for plant growth
    vpd_slope = -0.3
    rs_min = 70
    rcan_max = 1000000

    # **net radiation canopy******************************************************
    # constants or predefined:
    vp_slope = 0.14
    vp_offset = 0.34
    lw_slope = 1.35
    lw_offset = 0.35
    int_max = 0.2

    # **canopy resistance***********************************************************
    # constants or predefined:
    z_obs = 100
    z_b = 100
    z0m_bare = 0.001
    r0_bare = 0.38
    r0_full = 0.18
    tenacity = 1.5
    disp_bare = 0.0
    disp_full = 0.667
    fraction_h_bare = 0.65
    fraction_h_full = 0.95
    z0m_full = 0.1

    # **initial canopy aerodynamic resistance***********************************************************
    # constants or predefined:
    ndvi_obs_min = 0.25
    ndvi_obs_max = 0.75
    obs_fr = 0.25
    dem_resolution = 250

    # **ETLook.unstable.initial_friction_velocity_daily***********************************************************
    # constants or predefined:
    c1 = 1

    # **ETLook.unstable.transpiration***********************************************************
    # constants or predefined:
    iter_h = 3

    # **ETLook.resistance.soil_resistance***********************************************************
    # constants or predefined:
    r_soil_pow = -2.1
    r_soil_min = 800


    # **ETLook.unstable.initial_sensible_heat_flux_soil_daily***********************************************************
    # constants or predefined:
    #porosity = 0.4 #Note: soil dependent
    #se_top = 1.0 #Note should be input !
    rn_slope = 0.92
    rn_offset = -61.0

    # **ETLook.unstable.evaporation***********************************************************
    # constants or predefined:
    r0_grass = 0.23

    ######################## MODEL ETLOOK #########################################

    # **effective_leaf_area_index**************************************************
    # When using Copernicus data LAI and vegetation cover are estimated using SNAP biophysical
    # processor
    #vc = ETLook.leaf.vegetation_cover(ndvi, nd_min, nd_max, vc_pow)
    #lai = ETLook.leaf.leaf_area_index(vc, vc_min, vc_max, lai_pow)
    lai_eff = ETLook.leaf.effective_leaf_area_index(lai)

    vc[np.isnan(QC)] = np.nan
    lai[np.isnan(QC)] = np.nan
    lai_eff[np.isnan(QC)] = np.nan
    if out.vc == 1:
        PF.Save_as_tiff(vc_filename, vc, geo_ex, proj_ex)
    if out.lai == 1:
        PF.Save_as_tiff(lai_filename, lai, geo_ex, proj_ex)
    if out.lai_eff == 1:
        PF.Save_as_tiff(lai_eff_filename, lai_eff, geo_ex, proj_ex)


    #*******TRANSPIRATION COMPONENT****************************************************************

    # **soil fraction**************************************************************
    sf_soil = ETLook.radiation.soil_fraction(lai)

    sf_soil[np.isnan(QC)] = np.nan
    if out.sf_soil == 1:
        PF.Save_as_tiff(sf_soil_filename, sf_soil, geo_ex, proj_ex)

    # **atmospheric canopy resistance***********************************************
    iesd = ETLook.solar_radiation.inverse_earth_sun_distance(doy)
    sc = ETLook.solar_radiation.seasonal_correction(doy)
    day_angle = ETLook.clear_sky_radiation.day_angle(doy)
    decl = ETLook.solar_radiation.declination(doy)
    lat = ETLook.solar_radiation.latitude_rad(lat_deg)
    slope = ETLook.solar_radiation.slope_rad(slope_deg)
    aspect = ETLook.solar_radiation.aspect_rad(aspect_deg)
    ra_24_toa = ETLook.solar_radiation.daily_solar_radiation_toa(sc, decl, iesd, lat, slope, aspect)
    ws = ETLook.solar_radiation.sunset_hour_angle(lat, decl)
    ra_24_toa_flat = ETLook.solar_radiation.daily_solar_radiation_toa_flat(decl, iesd, lat, ws)
    diffusion_index = ETLook.solar_radiation.diffusion_index(trans_24, diffusion_slope, diffusion_intercept)

    # choose one of the two options below
    #ra_24 = ETLook.solar_radiation.daily_solar_radiation_flat(ra_24_toa_flat, trans_24)
    ra_24 = ETLook.solar_radiation.daily_total_solar_radiation(ra_24_toa, ra_24_toa_flat, diffusion_index, trans_24)
    stress_rad = ETLook.stress.stress_radiation(ra_24)
    p_air_24 = ETLook.meteo.air_pressure_daily(z, p_air_0_24)
    vp_24 = ETLook.meteo.vapour_pressure_from_specific_humidity_daily(qv_24, p_air_24)
    svp_24 = ETLook.meteo.saturated_vapour_pressure_average(
                ETLook.meteo.saturated_vapour_pressure_maximum(t_air_max_24),
                ETLook.meteo.saturated_vapour_pressure_minimum(t_air_min_24))
    vpd_24 = ETLook.meteo.vapour_pressure_deficit_daily(svp_24, vp_24)
    stress_vpd = ETLook.stress.stress_vpd(vpd_24, vpd_slope)
    stress_temp = ETLook.stress.stress_temperature(t_air_24, t_opt, t_min, t_max)
    r_canopy_0 = ETLook.resistance.atmospheric_canopy_resistance(lai_eff, stress_rad, stress_vpd, stress_temp, rs_min, rcan_max)

    # Save as tiff files
    lat[np.isnan(QC)] = np.nan
    slope[np.isnan(QC)] = np.nan
    aspect[np.isnan(QC)] = np.nan
    ra_24_toa[np.isnan(QC)] = np.nan
    ws[np.isnan(QC)] = np.nan
    ra_24_toa_flat[np.isnan(QC)] = np.nan
    diffusion_index[np.isnan(QC)] = np.nan
    ra_24[np.isnan(QC)] = np.nan
    stress_rad[np.isnan(QC)] = np.nan
    p_air_24[np.isnan(QC)] = np.nan
    vp_24[np.isnan(QC)] = np.nan
    svp_24[np.isnan(QC)] = np.nan
    vpd_24[np.isnan(QC)] = np.nan
    stress_vpd[np.isnan(QC)] = np.nan
    stress_temp[np.isnan(QC)] = np.nan
    r_canopy_0[np.isnan(QC)] = np.nan
    if out.lat == 1:
        PF.Save_as_tiff(lat_filename, lat, geo_ex, proj_ex)
    if out.slope == 1:
        PF.Save_as_tiff(slope_filename, slope, geo_ex, proj_ex)
    if out.aspect == 1:
        PF.Save_as_tiff(aspect_filename, aspect, geo_ex, proj_ex)
    if out.ws == 1:
        PF.Save_as_tiff(ws_filename, ws, geo_ex, proj_ex)
    if out.ra_24_toa == 1:
        PF.Save_as_tiff(ra_24_toa_filename, ra_24_toa, geo_ex, proj_ex)
    if out.diffusion_index == 1:
        PF.Save_as_tiff(diffusion_index_filename, diffusion_index, geo_ex, proj_ex)
    if out.ra_24 == 1:
        PF.Save_as_tiff(ra_24_filename, ra_24, geo_ex, proj_ex)
    if out.stress_rad == 1:
        PF.Save_as_tiff(stress_rad_filename, stress_rad, geo_ex, proj_ex)
    if out.p_air_24 == 1:
        PF.Save_as_tiff(p_air_24_filename, p_air_24, geo_ex, proj_ex)
    if out.vp_24 == 1:
        PF.Save_as_tiff(vp_24_filename, vp_24, geo_ex, proj_ex)
    if out.svp_24 == 1:
        PF.Save_as_tiff(svp_24_filename, svp_24, geo_ex, proj_ex)
    if out.vpd_24 == 1:
        PF.Save_as_tiff(vpd_24_filename, vpd_24, geo_ex, proj_ex)
    if out.stress_vpd == 1:
        PF.Save_as_tiff(stress_vpd_filename, stress_vpd, geo_ex, proj_ex)
    if out.stress_temp == 1:
        PF.Save_as_tiff(stress_temp_filename, stress_temp, geo_ex, proj_ex)
    if out.r_canopy_0 == 1:
        PF.Save_as_tiff(r_canopy_0_filename, r_canopy_0, geo_ex, proj_ex)

    # **net radiation canopy******************************************************
    t_air_k_24 = ETLook.meteo.air_temperature_kelvin_daily(t_air_24)
    # select one of the below two
    #l_net = ETLook.radiation.longwave_radiation_fao_etref(t_air_k_24, vp_24, trans_24)
    l_net = ETLook.radiation.longwave_radiation_fao(t_air_k_24, vp_24, trans_24, vp_slope, vp_offset, lw_slope, lw_offset)
    int_mm = ETLook.evapotranspiration.interception_mm(P_24, vc, lai, int_max)
    lh_24 = ETLook.meteo.latent_heat_daily(t_air_24)
    int_wm2 = ETLook.radiation.interception_wm2(int_mm, lh_24)
    rn_24 = ETLook.radiation.net_radiation(r0, ra_24, l_net, int_wm2)
    rn_24_canopy = ETLook.radiation.net_radiation_canopy(rn_24, sf_soil)

    # Save as tiff files
    t_air_k_24[np.isnan(QC)] = np.nan
    l_net[np.isnan(QC)] = np.nan
    int_mm[np.isnan(QC)] = np.nan
    lh_24[np.isnan(QC)] = np.nan
    int_wm2[np.isnan(QC)] = np.nan
    rn_24[np.isnan(QC)] = np.nan
    rn_24_canopy[np.isnan(QC)] = np.nan
    if out.t_air_k_24 == 1:
        PF.Save_as_tiff(t_air_k_24_filename, t_air_k_24, geo_ex, proj_ex)
    if out.l_net == 1:
        PF.Save_as_tiff(l_net_filename, l_net, geo_ex, proj_ex)
    if out.int_mm == 1:
        PF.Save_as_tiff(int_mm_filename, int_mm, geo_ex, proj_ex)
    if out.lh_24 == 1:
        PF.Save_as_tiff(lh_24_filename, lh_24, geo_ex, proj_ex)
    if out.int_wm2 == 1:
        PF.Save_as_tiff(int_wm2_filename, int_wm2, geo_ex, proj_ex)
    if out.rn_24 == 1:
        PF.Save_as_tiff(rn_24_filename, rn_24, geo_ex, proj_ex)
    if out.rn_24_canopy == 1:
        PF.Save_as_tiff(rn_24_canopy_filename, rn_24_canopy, geo_ex, proj_ex)

    # **canopy resistance***********************************************************

    t_air_k_i = ETLook.meteo.air_temperature_kelvin_inst(t_air_i)
    vp_i = ETLook.meteo.vapour_pressure_from_specific_humidity_inst(qv_i, p_air_i)
    ad_moist_i = ETLook.meteo.moist_air_density_inst(vp_i, t_air_k_i)
    ad_dry_i = ETLook.meteo.dry_air_density_inst(p_air_i, vp_i, t_air_k_i)
    ad_i = ETLook.meteo.air_density_inst(ad_dry_i, ad_moist_i)
    u_b_i_bare = ETLook.soil_moisture.wind_speed_blending_height_bare(u_i, z0m_bare, z_obs, z_b)
    lon = ETLook.solar_radiation.longitude_rad(lon_deg)
    ha = ETLook.solar_radiation.hour_angle(sc, dtime, lon)
    I0 = ETLook.clear_sky_radiation.solar_constant()
    ied = ETLook.clear_sky_radiation.inverse_earth_sun_distance(day_angle)
    h0 = ETLook.clear_sky_radiation.solar_elevation_angle(lat, decl, ha)
    h0ref = ETLook.clear_sky_radiation.solar_elevation_angle_refracted(h0)
    m = ETLook.clear_sky_radiation.relative_optical_airmass(p_air_i, p_air_0_i, h0ref)
    rotm = ETLook.clear_sky_radiation.rayleigh_optical_thickness(m)
    Tl2 = ETLook.clear_sky_radiation.linke_turbidity(wv_i, aod550_i, p_air_i, p_air_0_i)
    G0 = ETLook.clear_sky_radiation.extraterrestrial_irradiance_normal(I0, ied)
    B0c = ETLook.clear_sky_radiation.beam_irradiance_normal_clear(G0, Tl2, m, rotm, h0)
    Bhc = ETLook.clear_sky_radiation.beam_irradiance_horizontal_clear(B0c, h0)
    Dhc = ETLook.clear_sky_radiation.diffuse_irradiance_horizontal_clear(G0, Tl2, h0)
    ra_hor_clear_i = ETLook.clear_sky_radiation.ra_clear_horizontal(Bhc, Dhc)
    emiss_atm_i = ETLook.soil_moisture.atmospheric_emissivity_inst(vp_i, t_air_k_i)
    rn_bare = ETLook.soil_moisture.net_radiation_bare(ra_hor_clear_i, emiss_atm_i, t_air_k_i, lst, r0_bare)
    rn_full = ETLook.soil_moisture.net_radiation_full(ra_hor_clear_i, emiss_atm_i, t_air_k_i, lst, r0_full)
    h_bare = ETLook.soil_moisture.sensible_heat_flux_bare(rn_bare, fraction_h_bare)
    h_full = ETLook.soil_moisture.sensible_heat_flux_full(rn_full, fraction_h_full)
    u_b_i_full = ETLook.soil_moisture.wind_speed_blending_height_full_inst(u_i, z0m_full, z_obs, z_b)
    u_star_i_bare = ETLook.soil_moisture.friction_velocity_bare_inst(u_b_i_bare, z0m_bare, disp_bare, z_b)
    u_star_i_full = ETLook.soil_moisture.friction_velocity_full_inst(u_b_i_full, z0m_full, disp_full, z_b)
    L_bare = ETLook.soil_moisture.monin_obukhov_length_bare(h_bare, ad_i, u_star_i_bare, t_air_k_i)
    L_full = ETLook.soil_moisture.monin_obukhov_length_full(h_full, ad_i, u_star_i_full, t_air_k_i)
    u_i_soil = ETLook.soil_moisture.wind_speed_soil_inst(u_i, L_bare, z_obs)
    ras = ETLook.soil_moisture.aerodynamical_resistance_soil(u_i_soil)
    raa = ETLook.soil_moisture.aerodynamical_resistance_bare(u_i, L_bare, z0m_bare, disp_bare, z_obs)
    rac = ETLook.soil_moisture.aerodynamical_resistance_full(u_i, L_full, z0m_full, disp_full, z_obs)
    t_max_bare = ETLook.soil_moisture.maximum_temperature_bare(ra_hor_clear_i, emiss_atm_i, t_air_k_i, ad_i, raa, ras, r0_bare)
    t_max_full = ETLook.soil_moisture.maximum_temperature_full(ra_hor_clear_i, emiss_atm_i, t_air_k_i, ad_i, rac, r0_full)
    w_i = ETLook.soil_moisture.dew_point_temperature_inst(vp_i)
    t_dew_i = ETLook.soil_moisture.dew_point_temperature_inst(vp_i)
    t_wet_i = ETLook.soil_moisture.wet_bulb_temperature_inst(t_air_i, t_dew_i)
    t_wet_k_i = ETLook.meteo.wet_bulb_temperature_kelvin_inst(t_wet_i)
    lst_max = ETLook.soil_moisture.maximum_temperature(t_max_bare, t_max_full, vc)
    lst_min = ETLook.soil_moisture.minimum_temperature(t_wet_k_i, t_air_k_i, vc)
    se_root = ETLook.soil_moisture.soil_moisture_from_maximum_temperature(lst_max, lst, lst_min)
    stress_moist = ETLook.stress.stress_moisture(se_root, tenacity)
    r_canopy_0 = ETLook.resistance.atmospheric_canopy_resistance(lai_eff, stress_rad, stress_vpd, stress_temp, rs_min, rcan_max)
    r_canopy = ETLook.resistance.canopy_resistance(r_canopy_0, stress_moist, rcan_max)

    # Save as tiff files
    t_air_k_i[np.isnan(QC)] = np.nan
    vp_i[np.isnan(QC)] = np.nan
    ad_moist_i[np.isnan(QC)] = np.nan
    ad_dry_i[np.isnan(QC)] = np.nan
    ad_i[np.isnan(QC)] = np.nan
    u_b_i_bare[np.isnan(QC)] = np.nan
    lon[np.isnan(QC)] = np.nan
    ha[np.isnan(QC)] = np.nan
    h0[np.isnan(QC)] = np.nan
    h0ref[np.isnan(QC)] = np.nan
    m[np.isnan(QC)] = np.nan
    rotm[np.isnan(QC)] = np.nan
    Tl2[np.isnan(QC)] = np.nan
    B0c[np.isnan(QC)] = np.nan
    Bhc[np.isnan(QC)] = np.nan
    Dhc[np.isnan(QC)] = np.nan
    ra_hor_clear_i[np.isnan(QC)] = np.nan
    emiss_atm_i[np.isnan(QC)] = np.nan
    rn_bare[np.isnan(QC)] = np.nan
    rn_full[np.isnan(QC)] = np.nan
    u_b_i_full[np.isnan(QC)] = np.nan
    u_star_i_bare[np.isnan(QC)] = np.nan
    u_star_i_full[np.isnan(QC)] = np.nan
    u_i_soil[np.isnan(QC)] = np.nan
    ras[np.isnan(QC)] = np.nan
    raa[np.isnan(QC)] = np.nan
    rac[np.isnan(QC)] = np.nan
    t_max_bare[np.isnan(QC)] = np.nan
    t_max_full[np.isnan(QC)] = np.nan
    w_i[np.isnan(QC)] = np.nan
    t_dew_i[np.isnan(QC)] = np.nan
    t_wet_i[np.isnan(QC)] = np.nan
    t_wet_k_i[np.isnan(QC)] = np.nan
    lst_max[np.isnan(QC)] = np.nan
    se_root[np.isnan(QC)] = np.nan
    stress_moist[np.isnan(QC)] = np.nan
    r_canopy_0[np.isnan(QC)] = np.nan
    r_canopy[np.isnan(QC)] = np.nan
    if out.t_air_k_i == 1:
        PF.Save_as_tiff(t_air_k_i_filename, t_air_k_i, geo_ex, proj_ex)
    if out.vp_i == 1:
        PF.Save_as_tiff(vp_i_filename, vp_i, geo_ex, proj_ex)
    if out.ad_moist_i == 1:
        PF.Save_as_tiff(ad_moist_i_filename, ad_moist_i, geo_ex, proj_ex)
    if out.ad_dry_i == 1:
        PF.Save_as_tiff(ad_dry_i_filename, ad_dry_i, geo_ex, proj_ex)
    if out.ad_i == 1:
        PF.Save_as_tiff(ad_i_filename, ad_i, geo_ex, proj_ex)
    if out.u_b_i_bare == 1:
        PF.Save_as_tiff(u_b_i_bare_filename, u_b_i_bare, geo_ex, proj_ex)
    if out.lon == 1:
        PF.Save_as_tiff(lon_filename, lon, geo_ex, proj_ex)
    if out.ha == 1:
        PF.Save_as_tiff(ha_filename, ha, geo_ex, proj_ex)
    if out.ied == 1:
        PF.Save_as_tiff(ied_filename, ied, geo_ex, proj_ex)
    if out.h0 == 1:
        PF.Save_as_tiff(h0_filename, h0, geo_ex, proj_ex)
    if out.h0ref == 1:
        PF.Save_as_tiff(h0ref_filename, h0ref, geo_ex, proj_ex)
    if out.m == 1:
        PF.Save_as_tiff(m_filename, m, geo_ex, proj_ex)
    if out.rotm == 1:
        PF.Save_as_tiff(rotm_filename, rotm, geo_ex, proj_ex)
    if out.Tl2 == 1:
        PF.Save_as_tiff(Tl2_filename, Tl2, geo_ex, proj_ex)
    if out.B0c == 1:
        PF.Save_as_tiff(B0c_filename, B0c, geo_ex, proj_ex)
    if out.Bhc == 1:
        PF.Save_as_tiff(Bhc_filename, Bhc, geo_ex, proj_ex)
    if out.Dhc == 1:
        PF.Save_as_tiff(Dhc_filename, Dhc, geo_ex, proj_ex)
    if out.ra_hor_clear_i == 1:
        PF.Save_as_tiff(ra_hor_clear_i_filename, ra_hor_clear_i, geo_ex, proj_ex)
    if out.emiss_atm_i == 1:
        PF.Save_as_tiff(emiss_atm_i_filename, emiss_atm_i, geo_ex, proj_ex)
    if out.rn_bare == 1:
        PF.Save_as_tiff(rn_bare_filename, rn_bare, geo_ex, proj_ex)
    if out.rn_full == 1:
        PF.Save_as_tiff(rn_full_filename, rn_full, geo_ex, proj_ex)
    if out.u_b_i_full == 1:
        PF.Save_as_tiff(u_b_i_full_filename, u_b_i_full, geo_ex, proj_ex)
    if out.u_star_i_bare == 1:
        PF.Save_as_tiff(u_star_i_bare_filename, u_star_i_bare, geo_ex, proj_ex)
    if out.u_star_i_full == 1:
        PF.Save_as_tiff(u_star_i_full_filename, u_star_i_full, geo_ex, proj_ex)
    if out.u_i_soil == 1:
        PF.Save_as_tiff(u_i_soil_filename, u_i_soil, geo_ex, proj_ex)
    if out.ras == 1:
        PF.Save_as_tiff(ras_filename, ras, geo_ex, proj_ex)
    if out.raa == 1:
        PF.Save_as_tiff(raa_filename, raa, geo_ex, proj_ex)
    if out.rac == 1:
        PF.Save_as_tiff(rac_filename, rac, geo_ex, proj_ex)
    if out.t_max_bare == 1:
        PF.Save_as_tiff(t_max_bare_filename, t_max_bare, geo_ex, proj_ex)
    if out.t_max_full == 1:
        PF.Save_as_tiff(t_max_full_filename, t_max_full, geo_ex, proj_ex)
    if out.w_i == 1:
        PF.Save_as_tiff(w_i_filename, w_i, geo_ex, proj_ex)
    if out.t_dew_i == 1:
        PF.Save_as_tiff(t_dew_i_filename, t_dew_i, geo_ex, proj_ex)
    if out.t_wet_i == 1:
        PF.Save_as_tiff(t_wet_i_filename, t_wet_i, geo_ex, proj_ex)
    if out.t_wet_k_i == 1:
        PF.Save_as_tiff(t_wet_k_i_filename, t_wet_k_i, geo_ex, proj_ex)
    if out.lst_max == 1:
        PF.Save_as_tiff(lst_max_filename, lst_max, geo_ex, proj_ex)
    if out.se_root == 1:
        PF.Save_as_tiff(se_root_filename, se_root, geo_ex, proj_ex)
    if out.stress_moist == 1:
        PF.Save_as_tiff(stress_moist_filename, stress_moist, geo_ex, proj_ex)
    if out.r_canopy_0 == 1:
        PF.Save_as_tiff(r_canopy_0_filename, r_canopy_0, geo_ex, proj_ex)
    if out.r_canopy == 1:
        PF.Save_as_tiff(r_canopy_filename, r_canopy, geo_ex, proj_ex)

    # **initial canopy aerodynamic resistance***********************************************************

    # When using Copernicus data z_obst_max is already scaled with LAI in herbaceous landcovers
    #z_obst = ETLook.roughness.obstacle_height(ndvi, z_obst_max, ndvi_obs_min, ndvi_obs_max, obs_fr)
    z_obst = z_obst_max
    z_oro = ETLook.roughness.orographic_roughness(slope, dem_resolution) #careful - standard res is set to 250 # !!!
    z0m = ETLook.roughness.roughness_length(lai, z_oro, z_obst, z_obst_max, land_mask)
    ra_canopy_init = ETLook.neutral.initial_canopy_aerodynamic_resistance(u_24, z0m, z_obs)

    # Save as tiff files
    z_obst[np.isnan(QC)] = np.nan
    z_oro[np.isnan(QC)] = np.nan
    z0m[np.isnan(QC)] = np.nan
    ra_canopy_init[np.isnan(QC)] = np.nan
    if out.z_obst == 1:
        PF.Save_as_tiff(z_obst_filename, z_obst, geo_ex, proj_ex)
    if out.z_oro == 1:
        PF.Save_as_tiff(z_oro_filename, z_oro, geo_ex, proj_ex)
    if out.z0m == 1:
        PF.Save_as_tiff(z0m_filename, z0m, geo_ex, proj_ex)
    if out.ra_canopy_init == 1:
        PF.Save_as_tiff(ra_canopy_init_filename, ra_canopy_init, geo_ex, proj_ex)

    # **windspeed blending height daily***********************************************************
    u_b_24 = ETLook.meteo.wind_speed_blending_height_daily(u_24, z_obs, z_b)

    # Save as tiff files
    u_b_24[np.isnan(QC)] = np.nan
    if out.u_b_24 == 1:
        PF.Save_as_tiff(u_b_24_filename, u_b_24, geo_ex, proj_ex)

    # **ETLook.unstable.initial_friction_velocity_daily***********************************************************
    disp = ETLook.roughness.displacement_height(lai, z_obst, land_mask, c1)
    u_star_24_init = ETLook.unstable.initial_friction_velocity_daily(u_b_24, z0m, disp, z_b)

    # Save as tiff files
    disp[np.isnan(QC)] = np.nan
    u_star_24_init[np.isnan(QC)] = np.nan
    if out.disp == 1:
        PF.Save_as_tiff(disp_filename, disp, geo_ex, proj_ex)
    if out.u_star_24_init == 1:
        PF.Save_as_tiff(u_star_24_init_filename, u_star_24_init, geo_ex, proj_ex)

    # **ETLook.neutral.initial_daily_transpiration***********************************************************
    ad_dry_24 = ETLook.meteo.dry_air_density_daily(p_air_24, vp_24, t_air_k_24)
    ad_moist_24 = ETLook.meteo.moist_air_density_daily(vp_24, t_air_k_24)
    ad_24 = ETLook.meteo.air_density_daily(ad_dry_24, ad_moist_24)
    psy_24 = ETLook.meteo.psychrometric_constant_daily(p_air_24, lh_24)
    ssvp_24 = ETLook.meteo.slope_saturated_vapour_pressure_daily(t_air_24)
    t_24_init = ETLook.neutral.initial_daily_transpiration(rn_24_canopy, ssvp_24, ad_24, vpd_24, psy_24, r_canopy, ra_canopy_init)

    # Save as tiff files
    ad_dry_24[np.isnan(QC)] = np.nan
    ad_moist_24[np.isnan(QC)] = np.nan
    ad_24[np.isnan(QC)] = np.nan
    psy_24[np.isnan(QC)] = np.nan
    ssvp_24[np.isnan(QC)] = np.nan
    t_24_init[np.isnan(QC)] = np.nan
    if out.ad_dry_24 == 1:
        PF.Save_as_tiff(ad_dry_24_filename, ad_dry_24, geo_ex, proj_ex)
    if out.ad_moist_24 == 1:
        PF.Save_as_tiff(ad_moist_24_filename, ad_moist_24, geo_ex, proj_ex)
    if out.ad_24 == 1:
        PF.Save_as_tiff(ad_24_filename, ad_24, geo_ex, proj_ex)
    if out.psy_24 == 1:
        PF.Save_as_tiff(psy_24_filename, psy_24, geo_ex, proj_ex)
    if out.ssvp_24 == 1:
        PF.Save_as_tiff(ssvp_24_filename, ssvp_24, geo_ex, proj_ex)
    if out.t_24_init == 1:
        PF.Save_as_tiff(t_24_init_filename, t_24_init, geo_ex, proj_ex)

    # **ETLook.unstable.initial_sensible_heat_flux_canopy_daily***********************************************************
    h_canopy_24_init = ETLook.unstable.initial_sensible_heat_flux_canopy_daily(rn_24_canopy, t_24_init)

    # Save as tiff files
    h_canopy_24_init[np.isnan(QC)] = np.nan
    if out.h_canopy_24_init == 1:
        PF.Save_as_tiff(h_canopy_24_init_filename, h_canopy_24_init, geo_ex, proj_ex)

    # **ETLook.unstable.transpiration***********************************************************

    t_24 = ETLook.unstable.transpiration(rn_24_canopy, ssvp_24, ad_24, vpd_24, psy_24, r_canopy, h_canopy_24_init, t_air_k_24, u_star_24_init, z0m, disp, u_b_24, z_obs, z_b, iter_h)
    t_24_mm = ETLook.unstable.transpiration_mm(t_24, lh_24)

    # Save as tiff files
    t_24[np.isnan(QC)] = np.nan
    t_24_mm[np.isnan(QC)] = np.nan
    if out.t_24 == 1:
        PF.Save_as_tiff(t_24_filename, t_24, geo_ex, proj_ex)
    if out.t_24_mm == 1:
        PF.Save_as_tiff(t_24_mm_filename, t_24_mm, geo_ex, proj_ex)

    #*******EVAPORATION COMPONENT****************************************************************

    # **ETLook.radiation.net_radiation_soil***********************************************************
    sf_soil = ETLook.radiation.soil_fraction(lai)
    rn_24_soil = ETLook.radiation.net_radiation_soil(rn_24, sf_soil)

    # Save as tiff files
    sf_soil[np.isnan(QC)] = np.nan
    rn_24_soil[np.isnan(QC)] = np.nan
    if out.sf_soil == 1:
        PF.Save_as_tiff(sf_soil_filename, sf_soil, geo_ex, proj_ex)
    if out.rn_24_soil == 1:
        PF.Save_as_tiff(rn_24_soil_filename, rn_24_soil, geo_ex, proj_ex)

    # **ETLook.resistance.soil_resistance***********************************************************

    r_soil = ETLook.resistance.soil_resistance(se_top, land_mask, r_soil_pow, r_soil_min)

    # Save as tiff files
    r_soil[np.isnan(QC)] = np.nan
    if out.r_soil == 1:
        PF.Save_as_tiff(r_soil_filename, r_soil, geo_ex, proj_ex)

    # **ETLook.resistance.soil_resistance***********************************************************

    ra_soil_init = ETLook.neutral.initial_soil_aerodynamic_resistance(u_24, z_obs)

    # Save as tiff files
    ra_soil_init[np.isnan(QC)] = np.nan
    if out.ra_soil_init == 1:
        PF.Save_as_tiff(ra_soil_init_filename, ra_soil_init, geo_ex, proj_ex)

    # **ETLook.meteo.wind_speed_blending_height_daily***********************************************************

    u_b_24 = ETLook.meteo.wind_speed_blending_height_daily(u_24, z_obs, z_b)

    # Save as tiff files
    u_b_24[np.isnan(QC)] = np.nan
    if out.u_b_24 == 1:
        PF.Save_as_tiff(u_b_24_filename, u_b_24, geo_ex, proj_ex)

    # **ETLook.unstable.initial_friction_velocity_soil_daily***********************************************************

    u_star_24_soil_init = ETLook.unstable.initial_friction_velocity_soil_daily(u_b_24, disp, z_b)

    # Save as tiff files
    u_star_24_soil_init[np.isnan(QC)] = np.nan
    if out.u_star_24_soil_init == 1:
        PF.Save_as_tiff(u_star_24_soil_init_filename, u_star_24_soil_init, geo_ex, proj_ex)

    # **ETLook.unstable.initial_sensible_heat_flux_soil_daily***********************************************************

    stc = ETLook.radiation.soil_thermal_conductivity(se_top)
    vhc = ETLook.radiation.volumetric_heat_capacity(se_top, porosity)
    dd = ETLook.radiation.damping_depth(stc, vhc)
    g0_bs = ETLook.radiation.bare_soil_heat_flux(doy, dd, stc, t_amp_year, lat)
    g0_24 = ETLook.radiation.soil_heat_flux(g0_bs, sf_soil, land_mask, rn_24_soil, trans_24, ra_24, l_net, rn_slope, rn_offset)
    e_24_init = ETLook.neutral.initial_daily_evaporation(rn_24_soil, g0_24, ssvp_24, ad_24, vpd_24, psy_24, r_soil, ra_soil_init)
    h_soil_24_init = ETLook.unstable.initial_sensible_heat_flux_soil_daily(rn_24_soil, e_24_init, g0_24)

    # Save as tiff files
    g0_bs[np.isnan(QC)] = np.nan
    g0_24[np.isnan(QC)] = np.nan
    e_24_init[np.isnan(QC)] = np.nan
    h_soil_24_init[np.isnan(QC)] = np.nan
    if out.g0_bs == 1:
        PF.Save_as_tiff(g0_bs_filename, g0_bs, geo_ex, proj_ex)
    if out.g0_24 == 1:
        PF.Save_as_tiff(g0_24_filename, g0_24, geo_ex, proj_ex)
    if out.e_24_init == 1:
        PF.Save_as_tiff(e_24_init_filename, e_24_init, geo_ex, proj_ex)
    if out.h_soil_24_init == 1:
        PF.Save_as_tiff(h_soil_24_init_filename, h_soil_24_init, geo_ex, proj_ex)

    # **ETLook.unstable.evaporation***********************************************************

    e_24 = ETLook.unstable.evaporation(rn_24_soil, g0_24, ssvp_24, ad_24, vpd_24, psy_24, r_soil, h_soil_24_init, t_air_k_24, u_star_24_soil_init, disp, u_b_24, z_b, z_obs, iter_h)
    e_24_mm = ETLook.unstable.evaporation_mm(e_24, lh_24)
    et_24_mm = ETLook.evapotranspiration.et_actual_mm(e_24_mm, t_24_mm)

    # Save as tiff files
    e_24[np.isnan(QC)] = np.nan
    e_24_mm[np.isnan(QC)] = np.nan
    et_24_mm[np.isnan(QC)] = np.nan
    if out.e_24 == 1:
        PF.Save_as_tiff(e_24_filename, e_24, geo_ex, proj_ex)
    if out.e_24_mm == 1:
        PF.Save_as_tiff(e_24_mm_filename, e_24_mm, geo_ex, proj_ex)
    if out.et_24_mm == 1:
        PF.Save_as_tiff(et_24_mm_filename, et_24_mm, geo_ex, proj_ex)

    # **ETLook.unstable.evaporation***********************************************************
    rn_24_grass = ETLook.radiation.net_radiation_grass(ra_24, l_net, r0_grass)
    et_ref_24 = ETLook.evapotranspiration.et_reference(rn_24_grass, ad_24, psy_24, vpd_24, ssvp_24, u_24)
    et_ref_24_mm = ETLook.evapotranspiration.et_reference_mm(et_ref_24, lh_24)

    # Save as tiff files
    rn_24_grass[np.isnan(QC)] = np.nan
    et_ref_24[np.isnan(QC)] = np.nan
    et_ref_24_mm[np.isnan(QC)] = np.nan
    if out.rn_24_grass == 1:
        PF.Save_as_tiff(rn_24_grass_filename, rn_24_grass, geo_ex, proj_ex)
    if out.et_ref_24 == 1:
        PF.Save_as_tiff(et_ref_24_filename, et_ref_24, geo_ex, proj_ex)
    if out.et_ref_24_mm == 1:
        PF.Save_as_tiff(et_ref_24_mm_filename, et_ref_24_mm, geo_ex, proj_ex)
    return()
