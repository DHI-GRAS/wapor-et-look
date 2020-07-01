# temperature conversions
zero_celcius = 273.15     # 0 degrees C in K

# reference values
t_ref = 293.15            # reference temperature 20 degrees celcius
p_ref = 1013.25           # reference pressure in mbar
z_ref = 0                 # sea level m

lapse = -0.0065           # lapse rate K m-1
g = 9.807                 # gravity m s-2
gc_spec = 287.0           # gas constant J kg-1 K-1
gc_dry = 2.87             # dry air gas constant mbar K-1 m3 kg-1
gc_moist = 4.61           # moist air gas constant mbar K-1 m3 kg-1
r_mw = 0.622              # ratio water particles/ air particles
sh = 1004.0               # specific heat J kg-1 K-1
lh_0 = 2501000.0          # latent heat of evaporation at 0 C [J/kg]
lh_rate = -2361           # rate of latent heat vs temperature [J/kg/C]
power = (g/(-lapse*gc_spec))
k = 0.41                  # karman constant (-)
sol = 1367                # maximum solar radiation at top of atmosphere W m-2
sb = 5.67e-8              # stefan boltzmann constant
day_sec = 86400.0         # seconds in a day
year_sec = day_sec * 365  # seconds in a year

absorbed_radiation = 0.48 # biomass factor
conversion = 0.864        # conversion biomass calculation from g s-1 m-2
                          # to kg ha-1 d-1

z0_soil = 0.001           # soil roughness m
