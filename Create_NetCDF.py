import numpy as np
import netCDF4 as nc
import pandas as pd
import Emissions
import os
import csv
import DataLoader
import math
import xarray as xr


global avagadros
avagadros = 6.0221408e+23


def find_index(value, array):
  array = np.asarray(array)
  idx = (np.abs(array - value)).argmin()
  return idx

def find_time_index(timestamp, date_array, datesec_array):
  datesec_array = np.asarray(datesec_array)
  date_array = np.asarray(date_array)
  ts_date = timestamp.strftime("%Y%m%d")
  ts_datesec = get_number_of_seconds_since_midnight(timestamp)
  #find where datsec is exact match
  dates_with_matching_datesecs = date_array[np.where(datesec_array == ts_datesec)]
  #Check to see if any of those indices match the date
  indices_of_matching_date = np.where(dates_with_matching_datesecs == ts_date)
  return indices_of_matching_date.values[0]


def check_num_changes(changed_file, old_file):
  changed_dataset = nc.Dataset(changed_file, 'r')
  old_dataset = nc.Dataset(old_file, 'r')
  diff_times = [0,1,2,3,4,5,6, 7,8,9,10]
  for t in diff_times:
    diff = changed_dataset['Alumina_num'][t][:][:][:] - old_dataset['Alumina_num'][t][:][:][:]
    print(sum(sum(sum(changed_dataset['Alumina_num'][t][:][:][:]))))
    print(sum(sum(sum(old_dataset['Alumina_num'][t][:][:][:]))))
    print( sum(sum(sum(diff))))

def update_num_emissions(emission_netcdf_filepath, old_aerosol_diameter, new_aerosol_diameter, density, variable_name, output_file):
  old_mass_particle = (math.pi/6) * density * old_aerosol_diameter**3
  old_mass_particle_to_new_mass_particle_ratio =  old_mass_particle /  ( (math.pi/6) * density * new_aerosol_diameter**3)
  print(old_mass_particle_to_new_mass_particle_ratio)
  chunk_size = 10
  netcdf_data = nc.Dataset(emission_netcdf_filepath, 'r')
  output_data = nc.Dataset(output_file, 'w')
  num_iterations = math.ceil(len(netcdf_data[variable_name])/chunk_size)

  # copy dimensions
  for name in netcdf_data.dimensions:
    if name == 'time':
      output_data.createDimension(name, None)
    else:
      output_data.createDimension(name, (netcdf_data.dimensions[name].size))

  for name in netcdf_data.variables:
    if name == 'date' or name == 'datesec':
      print(name)
      output_data.createVariable(name, 'i4', ('time',))
    else:
      output_data.createVariable(name, netcdf_data.variables[name].dtype.str[1:], netcdf_data.variables[name].dimensions)
      output_data.variables[name].long_name = netcdf_data.variables[name].long_name
      if not name == variable_name:
        output_data.variables[name][:] = netcdf_data.variables[name][:]
        output_data.variables[name].unit = netcdf_data.variables[name].unit
      else:
        output_data.variables[name].units = netcdf_data.variables[name].units
  output_data.close()

  for i in range(num_iterations):
    output_data = nc.Dataset(output_file, 'r+')
    print(sum(sum(sum(netcdf_data.variables[variable_name][2]))))
    if chunk_size * (i+1) < netcdf_data.variables[variable_name].shape[0]:
      new_data = netcdf_data.variables[variable_name][chunk_size*i:chunk_size*(i+1)][:][:][:] *  old_mass_particle_to_new_mass_particle_ratio
      print(sum(sum(sum(new_data[2]))))
      output_data.variables['date'][chunk_size*i:chunk_size*(i+1)]  = netcdf_data['date'][chunk_size*i:chunk_size*(i+1)]
      output_data.variables['datesec'][chunk_size*i:chunk_size*(i+1)] = netcdf_data['datesec'][chunk_size*i:chunk_size*(i+1)]
      output_data.variables[variable_name][chunk_size*i:chunk_size*(i+1)] = new_data
      print(sum(sum(sum(output_data.variables[variable_name][2]))))
    else:
      new_data  = np.multiply(netcdf_data.variables[variable_name][chunk_size * i:][:][:][:], old_mass_particle_to_new_mass_particle_ratio)
      output_data.variables['date'][chunk_size * i:] = netcdf_data['date'][chunk_size * i:chunk_size * (i + 1)]
      output_data.variables['datesec'][chunk_size * i:] = netcdf_data['datesec'][chunk_size * i:chunk_size * (i + 1)]
      output_data.variables[variable_name][chunk_size * i:] = new_data

    output_data.close()
  netcdf_data.close()



def setup_netcdf(file_name, grid_lats, grid_lons, altitude_range, alt_int_vector, species_name, unit, molecular_weight):
  if os.path.exists(file_name):
    os.remove(file_name)
  ds = nc.Dataset(file_name, 'w', format='NETCDF4_CLASSIC')
  # Creating the file dimensions
  time = ds.createDimension('time', None)
  latitude = ds.createDimension('lat', len(grid_lats))
  longitude = ds.createDimension('lon', len(grid_lons))
  altitude = ds.createDimension('altitude', len(altitude_range))
  altitude_int = ds.createDimension('altitude_int', len(alt_int_vector))

  # Creating the file variables
  species = ds.createVariable(species_name, 'f4', ('time', 'altitude', 'lat', 'lon'))
  species.units = unit
  species.long_name = species_name
  species.molecular_weight = molecular_weight
  species.standard_name = 'Injection rate'
  date = ds.createVariable('date', 'i4', ('time',))
  date.units = 'YYYMMDD'
  date.long_name = 'date'
  datesec = ds.createVariable('datesec', 'i4', ('time',))
  datesec.units = 'seconds since midnight'
  datesec.long_name = "datesec"
  lats = ds.createVariable('lat', 'f4', ('lat',))
  lats.unit = 'degrees_north'
  lats.long_name = 'lat'
  lons = ds.createVariable('lon', 'f4', ('lon',))
  lons.unit = 'degrees_east'
  lons.long_name = 'lons'
  altitudes = ds.createVariable('altitude', 'f4', ('altitude',))
  altitudes.unit = 'km'
  altitudes.long_name = 'altitude'
  altitude_int = ds.createVariable('altitude_int', 'f4', ('altitude_int',))
  altitude_int.unit = 'km'
  altitude_int.long_name = 'altitude_interface'

  # Filling the file variables
  altitudes[:] = altitude_range
  altitude_int[:] = alt_int_vector
  lons[:] = grid_lons
  lats[:] = grid_lats

  # Adding very early date
  date[0] = 18000115
  datesec[0] = 1
  species[0, :, :, :] = np.zeros([1, len(altitudes[:]), len(lats[:]), len(lons[:])])

  ds.close()


def create_aluminum_emission_forecast(emission_data, satellite_ablation_profile, rb_ablation_profile, altitude_range,altitude_step, output_file_path, grid_filename, aerosol_emission_diameter, rebase_year=None):
  names = ['Alumina_mass', 'Alumina_num']
  units = ['molecules/cm3/s', '(particles/cm2/s)(molecules/mole)(g/kg)']
  molecular_weights = [135.064, 135.064]  # grams per mole for dst_a1
  densities = [0.0026, 0.0026]  # kg/cm^3 for dst_a1 (from kg/m^3 in physical properties file for dust)
  alt_int_vector = np.arange(min(altitude_range) - altitude_step / 2, max(altitude_range) + altitude_step,altitude_step)
  column_area, grid_lats, grid_lons = compute_grid_volumes(grid_filename, altitude_step, alt_int_vector)
  emission_duration = 1800
  df_groundtrack_heatmap = DataLoader.load_groundtrack_heatmap("../GroundtrackFiles/normalized_reentry_location_data.csv")

  # Preparing emission data
  emissions = Emissions.sample_forecasted_emissions(emission_data, satellite_ablation_profile, rb_ablation_profile, df_groundtrack_heatmap, altitude_range, rebase_year)

  # Saving emission characteristics
  emission_lats = [e.latitude for e in emissions]
  emission_lons = [e.longitude for e in emissions]
  emission_times = [e.rounded_time for e in emissions]
  emission_mass = [sum(e.emitted_mass) for e in emissions]
  emission_type = [e.type for e in emissions]
  emission_data = pd.DataFrame({'time': emission_times, 'latitude': emission_lats, 'longitude': emission_lons, 'emitted_mass': emission_mass,'emission_type': emission_type})
  emission_data.to_csv("./emission_data.csv")

  # Pre process to group emissions by emission date
  emissions_per_date = {}
  for e in emissions:
    if e.rounded_time not in emissions_per_date:
      emissions_per_date[e.rounded_time] = []
    emissions_per_date[e.rounded_time].append(e)

  sorted_emission_date_keys = sorted(emissions_per_date.keys())

  print("Finished preparing emission data")

  for species_index in range(len(names)):
    species_name = names[species_index]
    file_name = output_file_path + species_name + "_sampled_forecast.nc"
    setup_netcdf(file_name, grid_lats, grid_lons, altitude_range, alt_int_vector, species_name, units[species_index], molecular_weights[species_index])
    written_date_to_index = {}
    for key_index in range(len(sorted_emission_date_keys)):
        ds = nc.Dataset(file_name, 'r+')
        current_date_size = len(ds.variables['date'])
        current_species_size = len(ds.variables[species_name])

        emissions_at_date = emissions_per_date[sorted_emission_date_keys[key_index]]
        first_zero_pad = emissions_at_date[0].rounded_time - pd.Timedelta(seconds=1)
        start_emission_time = emissions_at_date[0].rounded_time
        end_emission_time = emissions_at_date[0].rounded_time + pd.Timedelta(seconds=(emission_duration - 1))
        last_zero_pad = emissions_at_date[0].rounded_time + pd.Timedelta(seconds=(emission_duration))
        e_times = [first_zero_pad, start_emission_time, end_emission_time, last_zero_pad]
        emission_time_indices = []
        overlap_count = 0

        for e_t in e_times:
          if e_t in written_date_to_index.keys():
            emission_time_indices.append(written_date_to_index[e_t])
            overlap_count = overlap_count +1
          else:
            emission_time_indices.append(current_date_size)
            written_date_to_index[e_t] = current_date_size
            current_date_size = current_date_size+1

        ds.variables['date'][min(emission_time_indices) : max(emission_time_indices) + 1] = [e.strftime("%Y%m%d") for e in e_times]
        ds.variables['datesec'][min(emission_time_indices) : max(emission_time_indices)+ 1] = [get_number_of_seconds_since_midnight(e) for e in e_times]

        emission_map = np.zeros((4, len(ds.variables['altitude'][:]), len(ds.variables['lat']), len(ds.variables['lon'][:])))
        for emission in emissions_at_date:
          e_latitude = emission.latitude
          e_longitude = emission.longitude
          e_lat_index = find_index(e_latitude, ds.variables['lat'][:])
          e_lon_index = find_index(e_longitude, ds.variables['lon'][:])

          for e_alt_index in range(len(emission.altitudes)):
            alt_index = find_index(emission.altitudes[e_alt_index], ds.variables['altitude'][:])
            grid_volume = lookup_grid_volume(e_latitude, e_longitude, column_area, grid_lats, grid_lons, emission.altitudes[e_alt_index], alt_int_vector)
            emission_flux = (emission.emitted_mass[e_alt_index]) * 1000 * avagadros / molecular_weights[species_index] / emission_duration / grid_volume #molecules_per_second_percm3

            if 'num' in species_name:
              mass_particle = (math.pi / 6) * densities[species_index] * aerosol_emission_diameter ** 3
              emission_flux = emission_flux * molecular_weights[species_index] / mass_particle #particles_per_second_per_cm3 * (molecules/mole) * (g/kg)

            emission_map[1][alt_index][e_lat_index][e_lon_index] = emission_flux
            emission_map[2][alt_index][e_lat_index][e_lon_index] = emission_flux

        if overlap_count > 0:
          ds.variables[species_name][min(emission_time_indices) : min(emission_time_indices) + overlap_count]  = np.add(emission_map[0 : overlap_count][:][:][:], ds.variables[species_name][min(emission_time_indices) : min(emission_time_indices) + overlap_count][:][:][:])
          ds.variables[species_name][min(emission_time_indices) + overlap_count : max(emission_time_indices) + 1]  = emission_map[overlap_count : ][:][:][:]
        else:
          ds.variables[species_name][min(emission_time_indices) : max(emission_time_indices) + 1]  = emission_map

        ds.close()
    print("Closing " + species_name)

    #Adding far out date to end of emission file
    ds = nc.Dataset(file_name, 'r+')
    current_date_size = len(ds.variables['date'])
    current_species_size = len(ds.variables[species_name])
    ds.variables['date'][current_date_size] = 21001215
    ds.variables['datesec'][current_date_size] = 86399
    ds.variables[species_name][current_species_size, :, :, :] = np.zeros([1, len(ds.variables['altitudes'][:]), len(ds.variables['lat'][:]), len(ds.variables['lon'][:])])
    ds.close()

    print("Finished")


def create_aluminum_emission_hist(species_fractions, df_reentry_tle_mass, satellite_ablation_profile,
                                  rb_ablation_profile, altitude_range, altitude_step, file_path, time_slice, ts,
                                  tle_expiration, grid_filename, emission_time_groups=None, rebase_year=None):
  names = ['Al', 'Alumina_mass', 'Alumina_num']
  units = ['atoms/s/cm3', 'molecules/s/cm3', '(particles/cm2/s)(molecules/mole)(g/kg)']
  molecular_weights = [26.9815, 135.064, 135.064]  # grams per mole
  densities = [0, 0.0026, 0.0026]  # kg/cm^3
  diameters = [0, 0.11, 0.11]  # cm^3
  emission_duration = 1800  # in seconds (default time step of WACCM)
  if emission_time_groups is not None:
    time_slice = [min(emission_time_groups), max(emission_time_groups)]
    emission_duration = [(emission_time_groups[i + 1] - emission_time_groups[i]).totalseconds for i in
                         range(len(emission_time_groups) - 1)]

  emissions = Emissions.determine_emissions(df_reentry_tle_mass, altitude_range, satellite_ablation_profile,
                                            rb_ablation_profile, ts, tle_expiration, emission_time_groups, time_slice)

  alt_int_vector = np.arange(min(altitude_range) - altitude_step / 2, max(altitude_range) + altitude_step,
                             altitude_step)
  column_area, grid_lats, grid_lons = compute_grid_volumes(grid_filename, altitude_step, alt_int_vector)
  for i in range(len(names)):
    print("Starting " + names[i])
    species_name = names[i]
    file_name = file_path + species_name + ".nc"
    if os.path.exists(file_name):
      os.remove(file_name)
    ds = nc.Dataset(file_name, 'w', format='NETCDF4')
    # Writing Global Attributes
    # ds.input_method = "SERIAL"

    # Creating the file dimensions
    time = ds.createDimension('time', None)
    latitude = ds.createDimension('lat', len(grid_lats))
    longitude = ds.createDimension('lon', len(grid_lons))
    altitude = ds.createDimension('altitude', len(altitude_range))
    altitude_int = ds.createDimension('altitude_int', len(alt_int_vector))

    # Creating the file variables
    species = ds.createVariable(species_name, 'f8', ('time', 'altitude', 'lat', 'lon'))
    species.units = units[i]
    species.long_name = species_name
    species.molecular_weight = molecular_weights[i]
    species.standard_name = 'Injection rate'
    date = ds.createVariable('date', 'i4', ('time',))
    date.units = 'YYYMMDD'
    date.long_name = 'date'
    datesec = ds.createVariable('datesec', 'i4', ('time',))
    datesec.units = 'seconds since midnight'
    datesec.long_name = "datesec"
    lats = ds.createVariable('lat', 'f8', ('lat',))
    lats.unit = 'degrees_north'
    lats.long_name = 'lat'
    lons = ds.createVariable('lon', 'f8', ('lon',))
    lons.unit = 'degrees_east'
    lons.long_name = 'lons'
    altitudes = ds.createVariable('altitude', 'f8', ('altitude',))
    altitudes.unit = 'km'
    altitudes.long_name = 'altitude'
    altitude_int = ds.createVariable('altitude_int', 'f8', ('altitude_int',))
    altitude_int.unit = 'km'
    altitude_int.long_name = 'altitude_interface'

    # Filling the file variables
    altitudes[:] = altitude_range
    altitude_int[:] = alt_int_vector

    lons[:] = grid_lons
    lats[:] = grid_lats

    if emission_time_groups is not None:
      # Grouped Emissions
      date[:] = [int(emission_time_groups[i].strftime("%Y%m%d")) for i in range(len(emission_time_groups) - 1)]
      if 'num' in species_name:
        # Number Emissions
        emission_map = construct_grouped_num_emission_map(emission_time_groups, emissions, species_fractions[i - 1],
                                                          altitudes, lats, lons, molecular_weights[i - 1], densities[i],
                                                          diameters[i], emission_duration, column_area, grid_lats,
                                                          grid_lons, alt_int_vector)
      else:
        # Mass Emissions
        emission_map = construct_grouped_emission_map(emission_time_groups, emissions, species_fractions[i], altitudes,
                                                      lats, lons, molecular_weights[i], emission_duration, column_area,
                                                      grid_lats, grid_lons, alt_int_vector)
    else:
      # Impulse Emissions
      emission_times = [e.rounded_time for e in emissions]
      if 'num' in species_name:
        # Number Emissions
        emission_map, zero_padded_dates, zero_padded_datesecs = construct_impulse_num_emission_map(emissions,
                                                                                                   emission_times,
                                                                                                   species_fractions[
                                                                                                     i - 1], altitudes,
                                                                                                   lats, lons,
                                                                                                   molecular_weights[
                                                                                                     i - 1],
                                                                                                   densities[i],
                                                                                                   diameters[i],
                                                                                                   emission_duration,
                                                                                                   column_area,
                                                                                                   grid_lats, grid_lons,
                                                                                                   alt_int_vector,
                                                                                                   rebase_year)
        date[:] = zero_padded_dates
        datesec[:] = zero_padded_datesecs
      else:
        # Mass Emissions
        emission_map, zero_padded_dates, zero_padded_datesecs = construct_impulse_emission_map(emissions,
                                                                                               emission_times,
                                                                                               species_fractions[i],
                                                                                               altitudes, lats, lons,
                                                                                               molecular_weights[i],
                                                                                               emission_duration,
                                                                                               column_area, grid_lats,
                                                                                               grid_lons,
                                                                                               alt_int_vector,
                                                                                               rebase_year)
        date[:] = zero_padded_dates
        datesec[:] = zero_padded_datesecs
    species[:, :, :, :] = emission_map
    ds.close()
    print("Closing " + species_name)
  zero_padded_dates
  return emissions


def construct_grouped_num_emission_map(emission_time_groups, emissions, species_fraction, altitudes, latitudes,
                                       longitudes, molecular_weight, density, diameter, emission_duration, column_area,
                                       grid_lats, grid_lons, alt_int_vector):
  emission_map = np.zeros((len(emission_time_groups), len(altitudes[:]), len(latitudes[:]), len(longitudes[:])))
  for time_index in range(len(emission_time_groups) - 1):
    for emission in emissions[time_index]:
      emission_lat = emission.latitude
      lat_index = find_index(emission_lat, latitudes[:])
      emission_long = emission.longitude + 180
      lon_index = find_index(emission_long, longitudes[:])
      for j in range(len(emission.altitudes)):
        alt_index = find_index(emission.altitudes[j], altitudes[:])
        grid_volume = lookup_grid_volume(emission_lat, emission_long, column_area, grid_lats, grid_lons,
                                         emission.altitudes[j], alt_int_vector)
        particles_per_second_percm3 = (emission.emitted_mass[
                                         j] * species_fraction) * 1000 * avagadros / molecular_weight / \
                                      emission_duration[time_index] / grid_volume  # molecules/cm^3/s
        mass_particle = density * np.pi * (diameter ** 3) / 6
        emission_number = particles_per_second_percm3 * molecular_weight / mass_particle
        emission_map[time_index][alt_index][lat_index][lon_index] = emission_number
  return emission_map


def construct_grouped_emission_map(emission_time_groups, emissions, species_fraction, altitudes, latitudes, longitudes,
                                   molecular_weight, emission_duration, column_area, grid_lats, grid_lons,
                                   alt_int_vector):
  emission_map = np.zeros((len(emission_time_groups), len(altitudes[:]), len(latitudes[:]), len(longitudes[:])))

  for time_index in range(len(emission_time_groups) - 1):
    for emission in emissions[time_index]:
      emission_lat = emission.latitude
      lat_index = find_index(emission_lat, latitudes[:])
      emission_long = emission.longitude + 180
      lon_index = find_index(emission_long, longitudes[:])
      for j in range(len(emission.altitudes)):
        alt_index = find_index(emission.altitudes[j], altitudes[:])
        grid_volume = lookup_grid_volume(emission_lat, emission_long, column_area, grid_lats, grid_lons,
                                         emission.altitudes[j], alt_int_vector)
        emission_map[time_index][alt_index][lat_index][lon_index] = (emission.emitted_mass[
                                                                       j] * species_fraction) * 1000 * avagadros / molecular_weight / \
                                                                    emission_duration[time_index] / grid_volume
  return emission_map


def construct_impulse_num_emission_map(emissions, species_fraction, altitudes, latitudes, longitudes,molecular_weight, density, diameter, emission_duration, column_area, grid_lats,grid_lons, alt_int_vector, written_dates, rebase_year):
  # Find overlapped emissions and new date emissions


  emission_map = np.zeros((len(emissions) * 4, len(altitudes[:]), len(latitudes[:]), len(longitudes[:])))
  zero_padded_dates = [0] * (len(emissions) * 4)
  zero_padded_datesecs = [0] * (len(emissions) * 4)
  emission_index = 0
  overlapping_emissions = []
  new_dates = []
  for time_index in range(0, len(emissions) * 4, 4):
    emission_lat = emissions[emission_index].latitude
    lat_index = find_index(emission_lat, latitudes[:])
    emission_long = emissions[emission_index].longitude + 180
    lon_index = find_index(emission_long, longitudes[:])
    emission_date_info = emissions[emission_index].rounded_time
    if rebase_year is not None and rebase_year > 1677:
      emission_date_info = emission_date_info.replace(year=rebase_year)
      emission_date = emission_date_info.strftime("%Y%m%d")
    elif rebase_year is not None and rebase_year < 1677:
      emission_date = emission_date_info.strftime("%Y%m%d")
      emission_date = str(rebase_year).rjust(4, '0') + emission_date[4:]

    # Writing the date
    zero_padded_dates[time_index] = int((emission_date_info - pd.Timedelta(seconds=1)).strftime("%Y%m%d"))
    zero_padded_dates[time_index + 1] = int(emission_date)
    zero_padded_dates[time_index + 2] = int((emission_date_info + pd.Timedelta(seconds=emission_duration - 1)).strftime("%Y%m%d"))
    zero_padded_dates[time_index + 3] = int((emission_date_info + pd.Timedelta(seconds=emission_duration)).strftime("%Y%m%d"))

    # Writing the datesec
    zero_padded_datesecs[time_index] = get_number_of_seconds_since_midnight(
      emission_date_info - pd.Timedelta(seconds=1))
    zero_padded_datesecs[time_index + 1] = get_number_of_seconds_since_midnight(emission_date_info)
    zero_padded_datesecs[time_index + 2] = get_number_of_seconds_since_midnight(
      emission_date_info + pd.Timedelta(seconds=emission_duration - 1))
    zero_padded_datesecs[time_index + 3] = get_number_of_seconds_since_midnight(
      emission_date_info + pd.Timedelta(seconds=emission_duration))

    for e_alt_index in range(len(emissions[emission_index].altitudes)):
      alt_index = find_index(emissions[emission_index].altitudes[e_alt_index], altitudes[:])
      grid_volume = lookup_grid_volume(emission_lat, emission_long, column_area, grid_lats, grid_lons,emissions[emission_index].altitudes[e_alt_index], alt_int_vector)
      molecules_per_second_percm3 = (emissions[emission_index].emitted_mass[e_alt_index] * species_fraction) * 1000 * avagadros / molecular_weight / emission_duration / grid_volume
      mass_particle = (math.pi / 6) * density * diameter ** 3
      particles_per_second_per_cm3 = molecules_per_second_percm3 * molecular_weight / mass_particle
      emission_map[time_index + 1][alt_index][lat_index][lon_index] = particles_per_second_per_cm3
      emission_map[time_index + 2][alt_index][lat_index][lon_index] = particles_per_second_per_cm3
    emission_index = emission_index + 1

  return emission_map, zero_padded_dates, zero_padded_datesecs, overlapping_emissions, new_dates


def construct_impulse_emission_map(emissions, emission_times, species_fraction, altitudes, latitudes, longitudes,
                                   molecular_weight, emission_duration, column_area, grid_lats, grid_lons,
                                   alt_int_vector, written_dates, rebase_year):
  # Find overlapped emissions and new date emissions
  new_dates = []
  overlapped_emissions = []
  unique_emissions = []
  for e in emissions:
    if e.rounded_time in written_dates or e.rounded_time in new_dates:
      overlapped_emissions.append(e)
    else:
      new_dates.append(e.rounded_time)
      unique_emissions.append(e)

  emission_map = np.zeros((len(unique_emissions) * 4, len(altitudes[:]), len(latitudes[:]), len(longitudes[:])))
  zero_padded_dates = [0] * (len(unique_emissions) * 4)
  zero_padded_datesecs = [0] * (len(unique_emissions) * 4)
  emission_index = 0
  overlapping_emissions = []
  new_dates = []
  for time_index in range(0, len(unique_emissions) * 4, 4):
    emission_lat = unique_emissions[emission_index].latitude
    lat_index = find_index(emission_lat, latitudes[:])
    emission_long = unique_emissions[emission_index].longitude + 180
    lon_index = find_index(emission_long, longitudes[:])
    emission_date_info = unique_emissions[emission_index].rounded_time
    if rebase_year is not None and rebase_year > 1677:
      emission_date_info = emission_date_info.replace(year=rebase_year)
      emission_date = emission_date_info.strftime("%Y%m%d")
    elif rebase_year is not None and rebase_year < 1677:
      emission_date = emission_date_info.strftime("%Y%m%d")
      emission_date = str(rebase_year).rjust(4, '0') + emission_date[4:]

    # Writing the date
    zero_padded_dates[time_index] = int((emission_date_info - pd.Timedelta(seconds=1)).strftime("%Y%m%d"))
    zero_padded_dates[time_index + 1] = int(emission_date)
    zero_padded_dates[time_index + 2] = int(
      (emission_date_info + pd.Timedelta(seconds=emission_duration - 1)).strftime("%Y%m%d"))
    zero_padded_dates[time_index + 3] = int(
      (emission_date_info + pd.Timedelta(seconds=emission_duration)).strftime("%Y%m%d"))
    # Writing the datesec
    zero_padded_datesecs[time_index] = get_number_of_seconds_since_midnight(
      emission_date_info - pd.Timedelta(seconds=1))
    zero_padded_datesecs[time_index + 1] = get_number_of_seconds_since_midnight(emission_date_info)
    zero_padded_datesecs[time_index + 2] = get_number_of_seconds_since_midnight(
      emission_date_info + pd.Timedelta(seconds=emission_duration - 1))
    zero_padded_datesecs[time_index + 3] = get_number_of_seconds_since_midnight(
      emission_date_info + pd.Timedelta(seconds=emission_duration))

    for e_alt_index in range(len(emissions[emission_index].altitudes)):
      alt_index = find_index(emissions[emission_index].altitudes[e_alt_index], altitudes[:])
      grid_volume = lookup_grid_volume(emission_lat, emission_long, column_area, grid_lats, grid_lons,
                                       emissions[emission_index].altitudes[e_alt_index], alt_int_vector)
      molecules_per_second_percm3 = (emissions[emission_index].emitted_mass[
                                       e_alt_index] * species_fraction) * 1000 * avagadros / molecular_weight / emission_duration / grid_volume
      emission_map[time_index + 1][alt_index][lat_index][lon_index] = molecules_per_second_percm3
      emission_map[time_index + 2][alt_index][lat_index][lon_index] = molecules_per_second_percm3
    emission_index = emission_index + 1

  # Handling Overlapped Emissions
  for o_emissions in overlapped_emissions:
    emitted_flux = []
    emission_lat = o_emissions.latitude
    emission_long = o_emissions.longitude + 180
    for e_alt_index in range(len(o_emissions.altitudes)):
      grid_volume = lookup_grid_volume(emission_lat, emission_long, column_area, grid_lats, grid_lons,
                                       emissions[emission_index].altitudes[e_alt_index], alt_int_vector)
      molecules_per_second_percm3 = (o_emissions.emitted_mass[
                                       e_alt_index] * species_fraction) * 1000 * avagadros / molecular_weight / emission_duration / grid_volume
      emitted_flux.append(molecules_per_second_percm3)
    o_emissions.set_emission_flux(emitted_flux)

  return emission_map, zero_padded_dates, zero_padded_datesecs, overlapping_emissions, new_dates


def get_number_of_days_since_1750(emission_date):
  reference_date = (pd.to_datetime('1750-01-01 00:00:00')).to_pydatetime()
  return float((emission_date.to_pydatetime() - reference_date).days) + (
            emission_date.to_pydatetime() - reference_date).seconds / 86400.0


def get_number_of_seconds_since_midnight(emission_date):
  return emission_date.second + 60 * emission_date.minute + 3600 * emission_date.hour


def create_constant_dust_emission(dust_constant_mass_value, dust_constant_number_value, altitude_range, altitude_step,
                                  lons_range, lats_range, duration, file_path):
  contant_values = [dust_constant_mass_value, dust_constant_number_value]
  filenames = ["t_constant_dust_mass_emission_", "t_constant_dust_number_emission_"]

  # Read in latitudes and longtitude based on grid
  dataset = nc.Dataset(
    "/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/WACCM-Cheyenne/Grid Template Files/coords_0.95x1.25_L70_c150828.nc")
  grid_lats = np.asarray(dataset['lat'])
  grid_lons = np.asarray(dataset['lon'])

  # Read in pressure levels
  pressure_levs_dataset = nc.Dataset(
    "../../WACCM-Cheyenne/WACCM Output/testing_impulse_dust/six_hour_testing/MW/cesm220_FWsc2000climo_f09_f09_mg17.cam.h3.0001-01-01-00000.nc")
  pressure_levs = np.asarray(pressure_levs_dataset['lev'])
  pressure_ilevs = np.asarray(pressure_levs_dataset['ilev'])

  for i in range(len(contant_values)):
    file_name = file_path + filenames[i] + str(duration[0].strftime("%Y")) + "-" + str(
      duration[1].strftime("%Y")) + ".nc"
    if os.path.exists(file_name):
      os.remove(file_name)
    os.chmod(file_path, 0o777)
    ds = nc.Dataset(file_name, 'w', format='NETCDF4')

    alt_int_vector = np.arange(min(altitude_range) - altitude_step / 2, max(altitude_range) + altitude_step,
                               altitude_step)

    # Creating the file dimensions
    time = ds.createDimension('time', None)
    latitude = ds.createDimension('lat', len(grid_lats))
    longitude = ds.createDimension('lon', len(grid_lons))
    # lev = ds.createDimension('lev', len(pressure_levs))
    # ilev = ds.createDimension('ilev', len(pressure_ilevs))
    alt = ds.createDimension('altitude', len(altitude_range))
    alt_int = ds.createDimension('altitude_int', len(alt_int_vector))

    # Creating the file variables
    species = ds.createVariable("dust", 'f8', ('time', 'altitude', 'lat', 'lon'))
    species.units = "(particles/cm3/s)(molecules/mole)(g/kg)"
    species.long_name = "AlSiO5 - dust - injection rate"
    species.molecular_weight = '135.064'

    date = ds.createVariable('date', 'i4', ('time',))
    date.units = 'YYYMMDD'
    date.long_name = 'date'

    datesec = ds.createVariable('datesec', 'i4', ('time',))
    datesec.units = 'seconds since midnight'
    datesec.long_name = "datesec"

    lats = ds.createVariable('lat', 'f8', ('lat',))
    lats.unit = 'degrees_north'
    lats.long_name = 'lat'
    lons = ds.createVariable('lon', 'f8', ('lon',))
    lons.unit = 'degrees_east'
    lons.long_name = 'lons'
    # levs = ds.createVariable('lev', 'f8', ('lev',))
    # ilevs = ds.createVariable('ilev', 'f8', ('ilev',))
    alts = ds.createVariable('altitude', 'f8', ('altitude',))
    alts_int = ds.createVariable('altitude_int', 'f8', ('altitude_int',))

    # Filling the file variables
    # levs[:] = pressure_levs
    # ilevs[:] = pressure_ilevs
    alts[:] = altitude_range
    alts_int[:] = alt_int_vector

    lons[:] = grid_lons
    lats[:] = grid_lats

    # Assuming constant emissions for the entire year
    print(duration)
    date[:] = [int(duration[i].strftime("%Y%m%d")) for i in range(len(duration))]
    datesec[:] = [int(get_number_of_seconds_since_midnight(duration[i])) for i in range(len(duration))]

    emission_map = np.zeros((len(date[:]), len(alts[:]), len(lats[:]), len(lons[:])))

    for time_index in range(2, len(date[:]) - 2):
      print(str(date[time_index]) + ": " + str(datesec[time_index]))
      # Assumed emission latitude to equator (0-1), alt = 50km
      emission_lat = find_index(grid_lats[2], lats[:])
      emission_alt = find_index(49, alts[:])
      emission_lon = find_index(grid_lons[3], lons[:])
      emission_map[time_index][emission_alt][emission_lat][emission_lon] = contant_values[i]

    # If emissions are shorter than the length of the run, need to add zero as last value
    species[:, :, :, :] = emission_map
    ds.close()


def compute_grid_volumes(grid_filename, altitude_depth_km, lower_boundary_altitudes):
  dataset = nc.Dataset(grid_filename)
  altitude_depth_cm = altitude_depth_km * 1e5
  gw = np.asarray(dataset['gw'])
  grid_lats = np.asarray(dataset['lat'])
  grid_lons = np.asarray(dataset['lon'])
  column_volume = np.zeros([len(lower_boundary_altitudes), len(grid_lats), len(grid_lons)])

  # These equations come from Mike Mills script
  for a in range(len(lower_boundary_altitudes)):
    earth_radius = 6.37122e8  # cm at the surface of Earth
    # altitudes are in km, so converting to cm with 1e5
    # SA_earth = 4.0 * np.pi * ( (earth_radius + a*1e5)** 2)
    SA_earth = 4.0 * np.pi * ((earth_radius) ** 2)
    constant = SA_earth / len(grid_lons) / float(sum(gw))
    for l in range(len(grid_lats)):
      column_volume[a, l, :] = np.ones([1, len(grid_lons)]) * gw[l] * constant * altitude_depth_cm
  save_to_file(column_volume, "./column_volumes.csv")
  return (column_volume, grid_lats, grid_lons)


def save_to_file(data, output_filename):
  data_shape = data.shape
  flatten_data = np.asarray(data)
  flatten_data = flatten_data.flatten()
  df = pd.DataFrame({'levs':data_shape[0], 'lat':data_shape[1], 'lon':data_shape[2], 'data':flatten_data})
  df.to_csv(output_filename)


def lookup_grid_volume(lat, lon, column_area, grid_lats, grid_lons, altitude, altitude_range):
  grid_lat_array = np.asarray(grid_lats)
  grid_lons_array = np.asarray(grid_lons)
  grid_alt_array = np.asarray(altitude_range)
  lat_idx = (np.abs(grid_lat_array - lat)).argmin()
  lon_idx = (np.abs(grid_lons_array - lon)).argmin()
  alt_idx = (np.abs(grid_alt_array - altitude)).argmin()
  return column_area[alt_idx, lat_idx, lon_idx]


def create_output_emission_csv(dust_constant_mass_value, dust_constant_number_value, altitude_range, lons_range,
                               lats_range, duration):
  contant_values = [dust_constant_mass_value, dust_constant_number_value]
  filepath = "/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/WACCM-Emissions/NCL Emission Generator/"
  filenames = ["altitudes.csv", "date.csv", "datesec.csv", "emission.csv"]
  pd.DataFrame(altitude_range).to_csv(filepath + filenames[0])
  dates = [int(duration[i].strftime("%Y%m%d")) for i in range(len(duration))]
  pd.DataFrame(dates).to_csv(filepath + filenames[1])
  datesec = [int(get_number_of_seconds_since_midnight(duration[i])) for i in range(len(duration))]
  pd.DataFrame(datesec).to_csv(filepath + filenames[2])

  f = open(filepath + filenames[3], 'w')
  writer = csv.writer(f)

  # If emissions are shorter than the length of the run, then the last date needs to be a zero ie len(times[:]) -1 in the loop
  # Zero padding before emissions
  for time_index in range(0, 2):
    print(str(dates[time_index]) + ": " + str(datesec[time_index]))
    # Assumed emission latitude to equator (0-1), alt = 50km
    e_lats = [-85, -90, -80]
    for e_lat in e_lats:
      emission_lat = find_index(e_lat, lats_range)
      emission_alt = find_index(50, altitude_range)
      e_lons = [0, 2.5, 5, 7.5]
      for j in e_lons:
        emission_lon = find_index(j, lons_range)
        data = [str(dates[time_index]), str(emission_alt), str(emission_lat), str(emission_lon), str(0)]
        writer.writerow(data)

  for time_index in range(2, len(dates) - 2):
    print(str(dates[time_index]) + ": " + str(datesec[time_index]))
    # emission_map = np.zeros((len(date[:]), len(altitudes[:]), len(lats[:]), len(lons[:])))
    # Assumed emission latitude to equator (0-1), alt = 50km
    e_lats = [-85, -90, -80]
    for e_lat in e_lats:
      emission_lat = find_index(e_lat, lats_range)
      emission_alt = find_index(50, altitude_range)
      e_lons = [0, 2.5, 5, 7.5]
      for j in e_lons:
        emission_lon = find_index(j, lons_range)
        data = [str(dates[time_index]), str(emission_alt), str(emission_lat), str(emission_lon),
                str(dust_constant_mass_value)]
        writer.writerow(data)

  for time_index in range(len(dates) - 2, len(dates)):
    print(str(dates[time_index]) + ": " + str(datesec[time_index]))
    # Assumed emission latitude to equator (0-1), alt = 50km
    e_lats = [-85, -90, -80]
    for e_lat in e_lats:
      emission_lat = find_index(e_lat, lats_range)
      emission_alt = find_index(50, altitude_range)
      e_lons = [0, 2.5, 5, 7.5]
      for j in e_lons:
        emission_lon = find_index(j, lons_range)
        data = [str(dates[time_index]), str(emission_alt), str(emission_lat), str(emission_lon), str(0)]
        writer.writerow(data)
  f.close()




