class Emission:
  def __init__(self):
    self.time = ""
    self.rounded_time = ""
    self.latitude = -1
    self.longitude = -1
    self.altitudes = []
    self.emitted_mass = []
    self.norad = 0
    self.type = ''
    self.emission_flux = []

  def setup(self, time, latitude, longitude, altitudes, emitted_mass, type):
    self.time = time
    self.rounded_time = time
    self.latitude = latitude
    self.longitude = longitude
    self.altitudes = altitudes
    self.emitted_mass = emitted_mass
    self.norad = 0
    self.type = type

  def set_emission_flux(self, flux):
    self.emission_flux = flux