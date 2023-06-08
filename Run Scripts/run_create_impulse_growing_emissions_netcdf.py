import DataLoader
import numpy as np
import Create_NetCDF
import datetime
import pandas as pd

import Emissions
import Plotter
import GrowEmissions
from skyfield.api import load
import math
import netCDF4 as nc

# Netcdf files can be made with impulse emissions
# Furthermore, netcdf files can be made using historical reentry location data or from sampling a reentry location distribution

def main():
  #Load necessary data for space debris aluminum emissions
  rb_ablation_profile = DataLoader.load_rb_ablation_profile("../esa_rocket_body_ablation_profile.csv")
  satellite_ablation_profile=DataLoader.load_satellite_ablation_profile("../esa_satellite_ablation_profile.csv")
  df_groundtrack_heatmap = DataLoader.load_groundtrack_heatmap("../GroundtrackFiles/normalized_reentry_location_data.csv")
  file_path ="../NetcdfFiles/"
  grid_data_file_path = "../../WACCM-Cheyenne/Grid Template Files/coords_0.95x1.25_L70_c150828.nc"
  emission_data = pd.read_csv('../Reentry_Forecasting/reentry_forecast_optimistic.csv')

  #Setting up the grids for the netcdf file (lat, lon, alt)
  altitude_step = 2 #Note that altitude step must be 1 km for the volume grid to match Mike Mills script, but can be any value. Verify the ablation profile sampling is correct with new value. Given in km
  altitude_range = np.arange(30, 94, altitude_step)

  rebase_year = 2000
  D_gn = 1.1e-7  #value from physical properities file for mam4 mode1 as defined in rad_climate in user_nl_cam (dgnum) units: m
  sigma_g = 1.6 #value from physical properities file for mam4 mode1 as defined in rad_climate in user_nl_cam (sigmag) units: none
  diameter = D_gn * math.exp(1.5 * (math.log(sigma_g))**2) #m
  chunking_size = 50
  num_files = 10
  #Impulse Emissions
  file_path_impulse = file_path + "ImpulseEmission/"
  #Create_NetCDF.create_aluminum_emission_forecast(emission_data, satellite_ablation_profile, rb_ablation_profile, altitude_range, altitude_step, file_path, grid_data_file_path, diameter, rebase_year)
  #Create_NetCDF.update_num_emissions(file_path + "Alumina_num_sampled_forecast_1_.nc", 7.805670344214621e-7, 8.873880134692096e-8, 2600, 'Alumina_num', file_path + "updated.nc")

  Emissions.sample_location(df_groundtrack_heatmap)
if __name__ == "__main__":
  main()
