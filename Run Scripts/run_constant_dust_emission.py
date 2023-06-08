import datetime
import pandas as pd
import numpy as np
import netCDF4 as nc
import Create_NetCDF
global avagadros
avagadros = 6.0221408e+23

def main():
  dust_constant_mass_value = 1.33e6/12 #yearly emitted mass
  dust_constant_number_value = dust_constant_mass_value
  part_cm3_s = dust_constant_number_value / 1000 /avagadros
  print("Part_cm3_s: " + str(part_cm3_s))
  print("Number emission: " + str(dust_constant_number_value))
  altitude_step = 2
  altitude_range = np.arange(30, 92, altitude_step)

  lats_range = list(np.arange(-90,90, 5)) #matches step size of plane et al 2021 aluminum paper
  lons_range = list(np.arange(0,360, 2.5)) #matches step size of plane et al 2021 aluminum paper


  file_path ="/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/WACCM-Emissions/NetcdfFiles/"
  #duration = [pd.to_datetime(datetime.datetime(1800, 1, 15, 00, 00, 1)), pd.to_datetime(datetime.datetime(1980, 1, 1, 1, 59, 59)), pd.to_datetime(datetime.datetime(1980, 1, 1, 2, 00, 00)),pd.to_datetime(datetime.datetime(1980, 1, 1, 4, 59, 59)), pd.to_datetime(datetime.datetime(1980, 1, 1, 5, 00, 00)), pd.to_datetime(datetime.datetime(2100, 12, 15, 23, 59, 59))]

  duration = [pd.to_datetime(datetime.datetime(1800, 1, 1, 00, 00, 1)), pd.to_datetime(datetime.datetime(2000, 1, 1, 1, 59, 59)), pd.to_datetime(datetime.datetime(2000, 1, 1, 2, 00, 00)),pd.to_datetime(datetime.datetime(2000, 1, 1, 2, 59, 59)), pd.to_datetime(datetime.datetime(2000, 1, 1, 3, 00, 00)), pd.to_datetime(datetime.datetime(2100, 12, 15, 23, 59, 59))]
  #duration = [pd.to_datetime(datetime.datetime(1800, 1, 1, 00, 00, 1)), pd.to_datetime(datetime.datetime(2000, 1, 1, 2, 00, 00)), pd.to_datetime(datetime.datetime(2000, 1, 1, 2, 30, 00)),pd.to_datetime(datetime.datetime(2000, 1, 1, 3, 00, 00)),  pd.to_datetime(datetime.datetime(2100, 12, 15, 23, 59, 59))]

  Create_NetCDF.create_constant_dust_emission(dust_constant_mass_value,dust_constant_number_value, altitude_range, altitude_step, lons_range, lats_range, duration,
                                file_path)

  #Create_NetCDF.create_output_emission_csv(dust_constant_mass_value,dust_constant_number_value, altitude_range, lons_range, lats_range, duration)

  # grid_data_file_path = "../../WACCM-Cheyenne/Grid Template Files/coords_0.95x1.25_L70_c150828.nc"
  # altitude_depth_km = 2
  # lower_boundary_altitudes = [48,50]
  # dataset = nc.Dataset(grid_data_file_path)
  # altitude_depth_cm = altitude_depth_km * 1e5
  # gw = np.asarray(dataset['gw'])
  # grid_lats = np.asarray(dataset['lat'])
  # grid_lons = np.asarray(dataset['lon'])
  # column_volume = np.zeros([len(lower_boundary_altitudes), len(grid_lats), len(grid_lons)])
  #
  # # These equations come from Mike Mills script
  # for a in range(len(lower_boundary_altitudes)):
  #   earth_radius = 6.37122e8  # cm at the surface of Earth
  #   # altitudes are in km, so converting to cm with 1e5
  #   # SA_earth = 4.0 * np.pi * ( (earth_radius + a*1e5)** 2)
  #   SA_earth = 4.0 * np.pi * ((earth_radius) ** 2)
  #   constant = SA_earth / len(grid_lons) / float(sum(gw))
  #   for l in range(len(grid_lats)):
  #     column_volume[a, l, :] = np.ones([1, len(grid_lons)]) * gw[l] * constant * altitude_depth_cm
  #
  # print()
if __name__ == "__main__":
  main()




