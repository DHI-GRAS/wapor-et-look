"""
===========================================
ETLook functions (:mod:`ETLook`)
===========================================

.. currentmodule:: ETLook

Within the :mod:`ETLook` module all physical and empirical
functions related to the calculation of the soil moisture, interceptin,
evaporation and transpiration are provided. These functions listed here
can be used to build function chains.

Instantaneous Radiation (:mod:`ETLook.clear_sky_radiation`)
==============================================================
.. automodule:: ETLook.clear_sky_radiation

.. autosummary::
   :toctree: generated/

   beam_irradiance_horizontal_clear
   beam_irradiance_normal_clear
   day_angle
   declination
   diffuse_irradiance_horizontal_clear
   extraterrestrial_irradiance_normal
   hour_angle
   inverse_earth_sun_distance
   linke_turbidity
   ra_clear_horizontal
   rayleigh_optical_thickness
   relative_optical_airmass
   solar_constant
   solar_elevation_angle
   solar_elevation_angle_refracted

Evapotranspiration (:mod:`ETLook.evapotranspiration`)
===========================================================================

.. automodule:: ETLook.evapotranspiration

.. autosummary::
    :toctree: generated/

    et_actual_mm
    et_reference
    et_reference_mm
    interception_mm

Vegetation Cover (:mod:`ETLook.leaf`)
===========================================================

.. automodule:: ETLook.leaf

.. autosummary::
   :toctree: generated/

   vegetation_cover
   leaf_area_index
   effective_leaf_area_index

Meteorology (:mod:`ETLook.meteo`)
==========================================================

.. automodule:: ETLook.meteo

.. autosummary::
   :toctree: generated/

   air_density
   air_density_daily
   air_density_inst
   air_pressure
   air_pressure_daily
   air_temperature_kelvin
   air_temperature_kelvin_daily
   air_temperature_kelvin_inst
   disaggregate_air_temperature
   disaggregate_air_temperature_daily
   disaggregate_air_temperature_inst
   disaggregate_dew_point_temperature_inst
   dry_air_density
   dry_air_density_daily
   dry_air_density_inst
   latent_heat
   latent_heat_daily
   moist_air_density
   moist_air_density_daily
   moist_air_density_inst
   psychrometric_constant
   psychrometric_constant_daily
   saturated_vapour_pressure
   saturated_vapour_pressure_minimum
   saturated_vapour_pressure_maximum
   saturated_vapour_pressure_average
   saturated_vapour_pressure_daily
   slope_saturated_vapour_pressure
   slope_saturated_vapour_pressure_daily
   vapour_pressure_deficit
   vapour_pressure_deficit_daily
   vapour_pressure_from_specific_humidity
   vapour_pressure_from_specific_humidity_daily
   vapour_pressure_from_specific_humidity_inst
   wet_bulb_temperature_kelvin_inst
   wind_speed_blending_height
   wind_speed_blending_height_daily



Net Available Energy (:mod:`ETLook.radiation`)
=================================================================

.. automodule:: ETLook.radiation

.. autosummary::
    :toctree: generated/

    bare_soil_heat_flux
    damping_depth
    interception_wm2
    longwave_radiation_fao
    longwave_radiation_fao_etref
    net_radiation
    net_radiation_canopy
    net_radiation_grass
    net_radiation_soil
    soil_fraction
    soil_heat_flux
    soil_thermal_conductivity
    volumetric_heat_capacity

Roughness (:mod:`ETLook.roughness`)
=============================================

.. automodule:: ETLook.roughness

.. autosummary::
   :toctree: generated/

   roughness_length
   obstacle_height
   displacement_height

Solar radiation (:mod:`ETLook.solar_radiation`)
=================================================================

.. automodule:: ETLook.solar_radiation

.. autosummary::
    :toctree: generated/

    aspect_rad
    cosine_solar_zenith_angle
    daily_solar_radiation_flat
    daily_solar_radiation_toa
    daily_solar_radiation_toa_flat
    daily_total_solar_radiation
    declination
    diffusion_index
    hour_angle
    inst_solar_radiation_toa
    inverse_earth_sun_distance
    latitude_rad
    seasonal_correction
    slope_rad
    sunset_hour_angle

Plant stress (:mod:`ETLook.stress`)
===================================================

.. automodule:: ETLook.stress

.. autosummary::
    :toctree: generated/

    stress_moisture
    stress_radiation
    stress_temperature
    stress_vpd

Canopy and Soil Resistance (:mod:`ETLook.resistance`)
===========================================================

.. automodule:: ETLook.resistance

.. autosummary::
    :toctree: generated/

    atmospheric_canopy_resistance
    canopy_resistance
    soil_resistance

Neutral Atmosphere (:mod:`ETLook.neutral`)
=====================================================

.. automodule:: ETLook.neutral

.. autosummary::
    :toctree: generated/

    initial_canopy_aerodynamic_resistance
    initial_daily_evaporation
    initial_daily_evaporation_mm
    initial_daily_transpiration
    initial_daily_transpiration_mm
    initial_soil_aerodynamic_resistance

Unstable Atmosphere (:mod:`ETLook.unstable`)
=======================================================

.. automodule:: ETLook.unstable

.. autosummary::
    :toctree: generated/

    evaporation
    evaporation_mm
    friction_velocity
    initial_friction_velocity_daily
    initial_sensible_heat_flux_canopy_daily
    initial_sensible_heat_flux_soil_daily
    monin_obukhov_length
    ra_canopy
    ra_soil
    stability_correction_heat_obs
    stability_factor
    stability_parameter
    stability_parameter_obs
    transpiration
    transpiration_mm


Soil Moisture (:mod:`ETLook.soil_moisture`)
====================================================================

.. automodule:: ETLook.soil_moisture

.. autosummary::
    :toctree: generated/

    psi_m
    psi_h
    aerodynamical_resistance_bare
    aerodynamical_resistance_full
    aerodynamical_resistance_soil
    atmospheric_emissivity_inst
    friction_velocity_bare_inst
    friction_velocity_full_inst
    initial_friction_velocity_inst
    minimum_temperature
    maximum_temperature
    maximum_temperature_bare
    maximum_temperature_full
    monin_obukhov_length_bare
    monin_obukhov_length_full
    net_radiation_bare
    net_radiation_full
    sensible_heat_flux_bare
    sensible_heat_flux_full
    soil_moisture_from_maximum_temperature
    wind_speed_blending_height_bare
    wind_speed_blending_height_full_inst
    wind_speed_soil_inst


"""
from pyWAPOR.ETLook import ETLook_code, solar_radiation, clear_sky_radiation, meteo, radiation, evapotranspiration, soil_moisture, leaf, stress, resistance, roughness, neutral, unstable, outputs

__all__ = ['ETLook_code', 'solar_radiation', 'clear_sky_radiation', 'meteo, radiation', 'evapotranspiration', 'soil_moisture', 'leaf', 'stress', 'resistance', 'roughness', 'neutral', 'unstable', 'outputs']

__version__ = '0.1'
