import numpy as np
import pandas as pd
from skyfield.api import wgs84
from skyfield.toposlib import ITRSPosition
from skyfield import positionlib as poslib
from skyfield import constants as constants
from Emission import Emission
import Propogator
import Plotter
import random

def find_index(value, array):
  array = np.asarray(array)
  idx = (np.abs(array - value)).argmin()
  return idx

def sample_date(emission_year):

  start = pd.Timestamp(year =emission_year, day =1, month = 1, hour = 0, second = 0, minute = 0)
  end = pd.Timestamp(year = emission_year, day = 31, month = 12, hour = 23, minute = 55, second = 0)
  delta = end - start
  diff_seconds = delta.days * 24*60*60 + delta.seconds
  random_second = random.randrange(start = 0, stop = diff_seconds)
  sampled_date = start + pd.Timedelta(seconds = random_second)

  return sampled_date
  #return pd.Timestamp(year = 2000, day = 1, month = 1, hour = 2, second = 0, minute = 0)

def create_cdf(df_groundtrack_heatmap):
  pdf = df_groundtrack_heatmap['normalized_bins']
  cdf = np.cumsum(pdf)
  return cdf

def sample_location(df_groundtrack_heatmap):
  cdf_distribution = create_cdf(df_groundtrack_heatmap)
  random_number = np.random.uniform(0,1, 1)
  index_of_closest_cdf_value = find_index(random_number, cdf_distribution)
  lat = df_groundtrack_heatmap['Bin Lat'].values[index_of_closest_cdf_value]
  lon =df_groundtrack_heatmap['Bin Long'].values[index_of_closest_cdf_value]
  lon = lon + abs(min(df_groundtrack_heatmap['Bin Long']))
  Plotter.plot_cum_dist(cdf_distribution)
  return (lat, lon)

def determine_Al_mass(reentry_object):
  # From SMAD, structure and mechanisms comprise, on average, 24% of dry mass. Pg 948 SMAD
  # Al Injected Mass per Year Calculation Assumes:
  # 1. that all of the structure and mechanism satellite mass is Aluminium, rho = 0.24 of dry mass
  # 3. rocket body frame is aluminium and structural mass minus engine mass is the aluminum mass
  # 4. rocket body engine is alpha = 470kg / (4000-470) kg = 0.13 engine to structure ratio (from Merlin 1D Dry Mass and Falcon 9 Upper Stage)

  alpha = 0.7
  rho = 0.21
  if reentry_object['Type'] == 'R':
    al = reentry_object['DryMass'] *  alpha
    search_for_centaur = reentry_object['Discos Object Name'].__contains__("CENTAUR")
    if search_for_centaur == True:
      search_for_centaur = not (reentry_object['Discos Object Name'].__contains__("PAYLOAD FAIRING"))
    search_for_delta = reentry_object['Discos Object Name'].__contains__("DELTA II")
    search_for_antares = reentry_object['Discos Object Name'].__contains__("ANTARES")
    if search_for_delta == True or search_for_centaur == True or search_for_antares == True:
      al = reentry_object['DryMass'] * 0.05
      #print("Ignoring: " + reentry_object['Discos Object Name'] + " : " + str(reentry_object['DryMass']))
  else:
    al = reentry_object['DryMass'] * rho
  return al


def update_geocentric_object(geocentric_object):
  position = wgs84.subpoint_of(geocentric_object)
  altitude_km = wgs84.height_of(geocentric_object).km
  new_geocentric_object = geocentric_object
  if altitude_km > 120:
    altitude_km = 120
    earth_position = wgs84.latlon(position.latitude.degrees, position.longitude.degrees, altitude_km * 1000)
    updated_gcrs_position = ITRSPosition(earth_position.itrs_xyz).at(geocentric_object.t)
    new_geocentric_object = poslib.build_position(updated_gcrs_position.position.km / constants.AU_KM,geocentric_object.velocity.km_per_s / constants.AU_KM * constants.DAY_S,t=geocentric_object.t, center=399)
  return new_geocentric_object

def calculate_reentry_trajectory(reentry_object, trajectory_calculator, ts):
  # Euler method to solve equations of motion assuming:
  # 1. Spherical, non-rotating Earth
  # 2. Constant Ballasitc Coefficient
  # 3. No initial accerlation of reentry object
  # 4. Cd is constant and equals 1.2
  # 5. Assuming all objects have a circular cross sectional area. Assumed a default diameter for rocket bodies and for satellites

  #Output File to save results
  output_filename = "./Trajectory Logging/" + str(reentry_object['norad']) + "_trajectory_output.txt"
  output_file = open(output_filename, "w")

  #Generating Parameters
  Cd = 1.2
  if  reentry_object['Diameter'] > 0:
    S = reentry_object['Diameter']**2 * (np.pi/4)
  else:
    if(reentry_object['Type']=='S'):
      uniform_density = 1.33/10e3 #assumes homogenous density in 1U cubesat weighing 1.33kg, cite nasa cubesat overview https://www.nasa.gov/mission_pages/cubesats/overview
    else:
      uniform_density = 2350 / (np.pi / 4 * 2.7**2 * 6.7) #assumes homogenous density in SL-4 Rocket Body (a kind representing the majority by mass distribution)
    volume = reentry_object['DryMass'] * uniform_density
    S = volume ** (2 / 3) * (np.pi / 4)

  B = reentry_object['DryMass'] / (Cd * S)
  dt = 0.1  # dt must be in terms of seconds
  maxIters = 100000 # max iters must be longer than 5 mins
  end_altitude = 20  # km
  mu = 3.986004418e14 #m3s-2

  #Force initial position to be at 120 km altitude
  geocentric_obj = update_geocentric_object(reentry_object['Geocentric Position Object'])
  params = {'R': 6378e3, 'L/D': 0, 'g0': 9.81, 'B': B, 'GeocentricObject':geocentric_obj, 'dt':dt, 'max_iters':maxIters, 'end_altitude':end_altitude, 'S':S, 'Cd':Cd, 'mu':mu, 'outputfile':output_file, 'ts':ts}

  #Generating initial states
  initial_position =  np.array(geocentric_obj.position.km *1000) #meters
  initial_velocity = np.array(geocentric_obj.velocity.km_per_s * 1000) # m/s
  # initial_height = wgs84.height_of(geocentric_obj).km * 1000 # m
  # initial_flight_path_angle = 0 #assumed fpa since uncontrolled reentry
  initial_density = 0
  # initial_gravity = 9.81 * (params['R'] / (params['R']+initial_height)) ** 2
  initial_states = [initial_position, initial_velocity, initial_density ]
  trajectory  = trajectory_calculator.forward_euler(initial_states, params)

  trajectory_data = pd.DataFrame(trajectory, columns = ['position','velocity', 'density','altitude', 'latitude','longitude'])
  output_file.close()
  return trajectory_data

def determine_mass_fraction_loss_for_altitude_range(altitude_range, ablation_profile):
  mass_fraction_loss_over_altitude_range = np.zeros(len(altitude_range))
  ablation_profile = ablation_profile.sort_values(by='Altitude')

  for i in range(0, len(altitude_range)):
    altitude_index = np.argmax(altitude_range[i] < ablation_profile['Altitude'])
    lower_altitude = ablation_profile['Altitude'][altitude_index - 1]
    upper_altitude = ablation_profile['Altitude'][altitude_index]
    if abs(lower_altitude - altitude_range[i]) < abs(upper_altitude - altitude_range[i]):
      mass_fraction_loss = ablation_profile['Al Mass Loss Fraction'][altitude_index - 1]
    else:
      mass_fraction_loss = ablation_profile['Al Mass Loss Fraction'][altitude_index]
    mass_fraction_loss_over_altitude_range[i] = mass_fraction_loss

  Plotter.plot_ablation_over_altitude(ablation_profile['Altitude'], ablation_profile['Al Mass Loss Fraction'], altitude_range, mass_fraction_loss_over_altitude_range)

  return mass_fraction_loss_over_altitude_range

def compute_altitude_injection(row, satellite_mass_fraction_loss, rocketbody_mass_fraction_loss):
  al_composition = determine_Al_mass(row)
  if row['Type'] in 'S':
    al_contribution =  np.multiply(al_composition, satellite_mass_fraction_loss)
  else:
    al_contribution = np.multiply(al_composition, rocketbody_mass_fraction_loss)
  return al_contribution

def calculate_total_al_injection_per_object(df_reentry_mass, satellite_ablation_profile, rb_ablation_profile):
  sat_total_al_mass_fraction_loss = np.trapz(satellite_ablation_profile['Al Mass Loss Fraction'],
                                             x=satellite_ablation_profile['Altitude'])
  rb_total_al_mass_fraction_loss = np.trapz(rb_ablation_profile['Al Mass Loss Fraction'],
                                            x=rb_ablation_profile['Altitude'])
  al_mass = []
  for index, row in df_reentry_mass.iterrows():
    if row['Type'] == 'R':
      al = determine_Al_mass(row) * 0.1 #rb_total_al_mass_fraction_loss
    else:
      al = determine_Al_mass(row) * sat_total_al_mass_fraction_loss
    al_mass.append(al)

  df_reentry_mass['Al Mass Contribution'] = al_mass
  return df_reentry_mass


def track_altitude_injection(altitude_range, df_reentry_mass, satellite_ablation_profile, rb_ablation_profile, year=None):
    df_track = df_reentry_mass.copy()
    satellite_al_contribution = np.zeros(len(altitude_range), dtype=np.dtype(float))
    rocket_body_al_contribution = np.zeros(len(altitude_range), dtype=np.dtype(float))

    if year is not None:
      df_track = df_reentry_mass.loc[df_reentry_mass['Year'] == year]

    satellite_mass_fraction_loss = determine_mass_fraction_loss_for_altitude_range(altitude_range, satellite_ablation_profile)
    rocketbody_mass_fraction_loss = determine_mass_fraction_loss_for_altitude_range(altitude_range, rb_ablation_profile)

    for index, row in df_track.iterrows():
      if row['Type'] in 'S':
        satellite_al_contribution = np.add(satellite_al_contribution,
                                           np.multiply(determine_Al_mass(row), satellite_mass_fraction_loss))
      else:
        rocket_body_al_contribution = np.add(rocket_body_al_contribution,
                                             np.multiply(determine_Al_mass(row), rocketbody_mass_fraction_loss))
    return (satellite_al_contribution, rocket_body_al_contribution)

def find_time_group_index(emission_time_groups, emitted_time):
  for time_index, e in enumerate(emission_time_groups):
    if emitted_time < e:
      return time_index-1

def round_reentry_time_to_nearest_30min(reentry_time):
  minutes = reentry_time.minute
  if minutes <14:
    rounded_time = reentry_time - pd.Timedelta(minutes = minutes)
  elif minutes >= 15 and minutes < 45:
    rounded_time = reentry_time + pd.Timedelta(minutes = 30 - minutes)
  else:
    rounded_time = reentry_time + pd.Timedelta(minutes= 60 - minutes)
  rounded_time = rounded_time - pd.Timedelta(seconds = rounded_time.second)
  return rounded_time

#Creating data for netCDF file - EMISSIONS
def determine_emissions(df_reentry_tle_mass, altitude_range, satellite_ablation_profile, rb_ablation_profile, ts, tle_expiration, emission_time_groups = None, time_slice = None):
  satellite_mass_fraction_loss = determine_mass_fraction_loss_for_altitude_range(altitude_range, satellite_ablation_profile)
  rocketbody_mass_fraction_loss = determine_mass_fraction_loss_for_altitude_range(altitude_range, rb_ablation_profile)

  if time_slice is not None:
    df_emissions = df_reentry_tle_mass.loc[pd.to_datetime(df_reentry_tle_mass['DECAY_EPOCH']) >= time_slice[0]].loc[pd.to_datetime(df_reentry_tle_mass['DECAY_EPOCH']) < time_slice[1]]
    df_emissions.sort_values(by = ['DECAY_EPOCH'], inplace=True)
  else:
    df_emissions = df_reentry_tle_mass

  if emission_time_groups is not None:
    #Grouped Emissions
    emissions = np.empty(shape=(len(emission_time_groups) - 1, 0), dtype='object').tolist()
    for index, row in df_emissions.iterrows():
      emitted = Emission()
      (lat, lon, reentry_time, error_code, geocentric) = Propogator.satellite_propogation(row, shouldRandomSample=False, ts=ts, tle_expiration=tle_expiration)
      emitted.latitude = float(lat.degrees)
      emitted.longitude = float(lon.degrees)
      emitted.rounded_time = round_reentry_time_to_nearest_30min(reentry_time)
      emitted.time = reentry_time
      emitted.altitudes = altitude_range
      emitted.norad = row['norad']
      if row['Type'] == 'S':
        emitted.emitted_mass = np.multiply(satellite_mass_fraction_loss, determine_Al_mass(row))
      else:
        emitted.emitted_mass = np.multiply(rocketbody_mass_fraction_loss, determine_Al_mass(row))
      time_group_index = find_time_group_index(emission_time_groups, reentry_time)
      emissions[time_group_index].append(emitted)
  else:
      #Impulse Emissions
      emissions = []
      for index, row in df_emissions.iterrows():
        emitted = Emission()
        (lat, lon, reentry_time, error_code, geocentric) = Propogator.satellite_propogation(row,shouldRandomSample=False,ts=ts,tle_expiration=tle_expiration)
        emitted.latitude = float(lat.degrees)
        emitted.longitude = float(lon.degrees)
        emitted.rounded_time = round_reentry_time_to_nearest_30min(reentry_time)
        emitted.time = reentry_time
        emitted.altitudes = altitude_range
        emitted.norad = row['norad']
        if row['Type'] == 'S':
          emitted.emitted_mass = np.multiply(satellite_mass_fraction_loss, determine_Al_mass(row))
        else:
          emitted.emitted_mass = np.multiply(rocketbody_mass_fraction_loss, determine_Al_mass(row))
        emissions.append(emitted)

  return emissions

def sample_forecasted_emissions(emission_data, satellite_ablation_profile, rb_ablation_profile, df_groundtrack_heatmap, altitude_range, rebase_year):
  satellite_mass_fraction_loss = determine_mass_fraction_loss_for_altitude_range(altitude_range, satellite_ablation_profile)
  rocketbody_mass_fraction_loss = determine_mass_fraction_loss_for_altitude_range(altitude_range, rb_ablation_profile)
  altitude_step = abs(altitude_range[1]-altitude_range[0])
  alpha = 0.7
  rho = 0.21
  total_reentry_events = sum(emission_data['Total Reentry Events'])
  satelltie_reentry_events = sum(emission_data['Total Reentry Events'].loc[emission_data['Type'].str.contains('S')])
  rocket_body_reentry_events = sum(emission_data['Total Reentry Events'].loc[emission_data['Type'].str.contains('R')])
  print("Total Reentry Events: " + str(total_reentry_events))
  print("Satellite Events: " + str(satelltie_reentry_events))
  print("Rocket Body Events: " + str(rocket_body_reentry_events))
  print("Satellite Mass: " + str(emission_data['Avg Mass (kg)'].loc[emission_data['Type'].str.contains('S')].values[0]))
  print("Rocket Body Mass: " + str(emission_data['Avg Mass (kg)'].loc[emission_data['Type'].str.contains('R')].values[0]))

  emissions = []
  for rocket_body_index in range(rocket_body_reentry_events):
    reentry_lat, reentry_lon = sample_location(df_groundtrack_heatmap)
    reentry_time = round_reentry_time_to_nearest_30min(sample_date(rebase_year))
    reentry_emitted_aluminum_mass = altitude_step * np.multiply(rocketbody_mass_fraction_loss, emission_data['Avg Mass (kg)'].loc[emission_data['Type'].str.contains('R')].values[0] * alpha)
    emitted = Emission()
    emitted.setup(latitude=reentry_lat, longitude=reentry_lon, time=reentry_time, altitudes=altitude_range,emitted_mass=reentry_emitted_aluminum_mass, type= 'R')
    emissions.append(emitted)
  for sat_index in range(satelltie_reentry_events):
    reentry_lat, reentry_lon = sample_location(df_groundtrack_heatmap)
    reentry_time = round_reentry_time_to_nearest_30min(sample_date(rebase_year))
    reentry_emitted_aluminum_mass =  altitude_step * np.multiply(satellite_mass_fraction_loss, emission_data['Avg Mass (kg)'].loc[emission_data['Type'].str.contains('S')].values[0] * rho)
    emitted = Emission()
    emitted.setup(latitude=reentry_lat, longitude=reentry_lon, time=reentry_time, altitudes=altitude_range,emitted_mass=reentry_emitted_aluminum_mass, type ='S')
    emissions.append(emitted)
  emissions.sort(key=lambda x: x.time, reverse=False)
  return emissions

