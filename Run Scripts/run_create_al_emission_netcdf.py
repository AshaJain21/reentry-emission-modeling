import DataLoader
import numpy as np
import Create_NetCDF
import datetime
import pandas as pd
import Plotter
from skyfield.api import load


# Netcdf files can be made with impulse emissions or with grouped emissions
# Furthermore, netcdf files can be made using historical reentry location data or from sampling a reentry location distribution

def impulse_emission(species_fractions, df_reentry_tle_mass_data, satellite_ablation_profile,
                                rb_ablation_profile, altitude_range, altitude_step, lons_range, lats_range,
                                file_path, time_slice, ts, tle_expiration, grid_data_file_path, rebase_year = None):
  emissions = Create_NetCDF.create_aluminum_emission_hist(species_fractions, df_reentry_tle_mass_data, satellite_ablation_profile,
                                rb_ablation_profile, altitude_range, altitude_step, lons_range, lats_range,
                                file_path, time_slice, ts, tle_expiration, grid_data_file_path, emission_time_groups=None, rebase_year = rebase_year)
  return emissions

def grouped_emissions(species_fractions, df_reentry_tle_mass_data,satellite_ablation_profile,rb_ablation_profile, altitude_range, altitude_step,lons_range, lats_range,file_path, time_slice, ts, tle_expiration, grid_data_file_path, emission_time_groups):
  emissions = Create_NetCDF.create_aluminum_emission_hist(species_fractions, df_reentry_tle_mass_data,satellite_ablation_profile,rb_ablation_profile, altitude_range, altitude_step,lons_range, lats_range,file_path, time_slice, ts, tle_expiration, grid_data_file_path, emission_time_groups)
  return emissions

def main():
  #Load necessary data for space debris aluminum emissions
  rb_ablation_profile = DataLoader.load_rb_ablation_profile("/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/WACCM-Emissions/esa_rocket_body_ablation_profile.csv")
  satellite_ablation_profile=DataLoader.load_satellite_ablation_profile("/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/WACCM-Emissions/esa_satellite_ablation_profile.csv")

  df_reentry_tle_mass_data, df_reentry_mass, df_reentries = DataLoader.read_in_clean_datasets("/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/WACCM-Emissions/")
  tle_expiration = 5  # 5 days after tle epoch, tle is considered bad
  ts = load.timescale()

  #Al2O3 Species Fraction
  species_fractions = (1)

  #Setting up the grids for the netcdf file (lat, lon, alt)
  altitude_step = 1 #Note that altitude step must be 1 km for the volume grid to match Mike Mills script, but can be any value. Given in km
  altitude_range = np.arange(30, 92, altitude_step)

  file_path ="/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/WACCM-Emissions/NetcdfFiles/"
  grid_data_file_path = "/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/WACCM-Cheyenne/Grid Template Files/coords_1.9x2.5_L88_c150828_netcdf4.nc"
  #Create array with times of interest for emissions with specific emission groupings (ie per month, per year, per day)
    #emission_time_groups = pd.date_range(start='2020-01-01', end='2020-01-01', periods = 12) #linear spacing of time for N "periods"
    #emission_time_groups = [pd.to_datetime(datetime.datetime(2020,1,1)) + pd.DateOffset(months = i) for i in range(13)]
  #emission_time_groups = [pd.to_datetime(datetime.datetime(2020,1,1)) + pd.DateOffset(years = i) for i in range(2)]
  time_slice = [pd.to_datetime(datetime.datetime(2014,1,1)), pd.to_datetime(datetime.datetime(2015,1,1))]

  #Grouped Emissions
  file_path_grouped = file_path + "GroupedEmission/"
  # emissions_object_list = grouped_emissions(species_fractions, df_reentry_tle_mass_data, satellite_ablation_profile, rb_ablation_profile,altitude_range, altitude_step, lons_range, lats_range, file_path_grouped, time_slice, ts, grid_data_file_path, tle_expiration,emission_time_groups)
  # Plotter.plot_emission_fast(emissions_object_list, emission_time_groups)

  #Impulse Emissions
  file_path_impulse = file_path + "ImpulseEmission/"
  rebase_year = None
  impulse_emission(species_fractions, df_reentry_tle_mass_data, satellite_ablation_profile,rb_ablation_profile, altitude_range, altitude_step, lons_range, lats_range,file_path_impulse, time_slice, ts, tle_expiration, grid_data_file_path, rebase_year = rebase_year)

if __name__ == "__main__":
    main()
