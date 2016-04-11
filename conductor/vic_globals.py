""" This module represents the VIC "global parameters file", essentially
  all of the settings that affect VIC at run time
"""

__all__ = ['Global']

from os.path import isdir, isfile, basename, dirname
from collections import OrderedDict
from datetime import date

class Scalar(object):
  def __init__(self, type_, value=None):
    self.type_ = type_
    self.value = value
  def __set__(self, instance, value):
    try:
      self.value = self.type_(value)
    except:
      raise ValueError("Cannot convert '{}' to type {}".format(value, self.type_))
  def __get__(self, instance, cls):
    return self.value
  def __str__(self, instance, cls, name=''):
    if self.value is None:
      return ''
    return '{} {}\n'.format(name.upper(), self.value)

class Boolean(Scalar):
  def __init__(self, value=None):
    super().__init__(bool, value)
  def __set__(self, instance, value):
    # catch strings which represent 'Falsy' values
    if value in ('FALSE', 'False', 'false', '0'):
      value = False
    super().__set__(instance, value)
  def __str__(self, instance, cls, name=''):
    value = 'TRUE' if self.value else 'FALSE'
    return '{} {}\n'.format(name.upper(), value)

class Filename(Scalar):
  def __init__(self, value=None):
    self.type_ = str
    self.value = value
  def __set__(self, instance, value):
    if not isdir(dirname(value)):
      raise ValueError("Cannot set parameter to a file in a non-existant \
        directory: {}".format(dirname(value)))
    super().__set__(instance, value)

class Mapping(object):
  def __init__(self):
    self.dict_ = {}
  def __set__(self, instance, value):
    try:
      key, value = value.split(None, 1)
    except:
      raise ValueError('Assigment to {} requires a key and value separated \
        by whitespace'.format(self.__class__.__name__))
    self.dict_[key] = value
  def __get__(self, instance, cls):
    return self.dict_
  def __str__(self, instance, cls, name=''):
    return '\n'.join('{} {} {}'.format(name.upper(), k, v) for k, v in \
      self.dict_.items()) + '\n'

class List(object):
  def __init__(self):
    self.value = []
  def __set__(self, instance, value):
    self.value.append(value)
  def __get__(self, instance, cls):
    return self.value
  def __str__(self, instance, cls, name=''):
    return '\n'.join('{} {}'.format(name.upper(), item) for item in self.value) + '\n'

class OutfileList(object):
  def __init__(self):
    self.value = OrderedDict()
  def __set__(self, instance, value):
    ''' Assume assignment will happen in order and just maintain
      state along the way
    '''
    try:
      filename, num_vars = value.split()
      self.value[filename] = []
    except ValueError: # append a variable
      last_key = list(self.value.keys())[-1]
      self.value[last_key].append(value)
  def __get__(self, instance, cls):
    return self.value
  def __str__(self, instance, cls, name=None):
    rv = ''
    for filename, varlist in self.value.items():
      rv += 'OUTFILE {} {}\n'.format(filename, len(varlist))
      for var in varlist:
        rv += 'OUTVAR {}\n'.format(var)
    return rv

class AttributeOrderDict(dict):
  """Dict-like object used for recording attribute definition order.
  """

  def __init__(self, no_special_methods=True, no_callables=True):
    self.member_order = []
    self.no_special_methods = no_special_methods
    self.no_callables = no_callables
    super().__init__()

  def __setitem__(self, key, value):
    skip = False
    if self.no_callables:
      if callable(value):
        skip = True
    # Skip special methods if not wanted.
    if self.no_special_methods:
      if key.startswith('__') and key.endswith('__'):
        skip = True
    # Skip properties
    if isinstance(value, property):
      skip = True
    if not skip:
      self.member_order.append(key)
    super().__setitem__(key, value)


class OrderedMeta(type):
  """Meta class that helps to record attribute definition order.
  """

  @classmethod
  def __prepare__(mcs, name, bases, **kwargs):
    return AttributeOrderDict(**kwargs)

  def __new__(mcs, name, bases, cdict, **kwargs):
    cls = type.__new__(mcs, name, bases, cdict)
    cls.member_order = cdict.member_order
    cls._closed = True
    return cls

  # Needed to use up kwargs.
  def __init__(cls, name, bases, cdict, **kwargs):
    super().__init__(name, bases, cdict)

  def __setattr__(cls, name, value):
    # Later attribute additions go through here.
    if getattr(cls, '_closed', False):
      raise AttributeError(
        'Cannot set attribute after class definition.')
    super().__setattr__(name, value)
  
class Global(metaclass = OrderedMeta):
  """An object representation of the global parameter file for the
     Variable Inflow and Capacity (VIC) hydrologic model. This object
     contains attributes for everything that affects an invocation of VIC.
  """
  time_step = Scalar(int)
  snow_step = Scalar(int)
  startyear = Scalar(int)
  startmonth = Scalar(int)
  startday = Scalar(int)
  starthour = Scalar(int)
  endyear = Scalar(int)
  endmonth = Scalar(int)
  endday = Scalar(int)
  full_energy = Boolean()
  frozen_soil = Boolean()
  no_flux = Boolean()
  dist_prcp = Boolean()
  corrprec = Boolean()
  min_wind_speed = Scalar(float)
  prec_expt = Scalar(float)
  #open_ground_id = Scalar(int)
  glacier_id = Scalar(int)
  glacier_accum_start_year = Scalar(int)
  glacier_accum_start_month = Scalar(int)
  glacier_accum_start_day = Scalar(int)
  glacier_accum_interval = Scalar(int)
  output_force = Boolean()
  init_state = Filename()
  statename = Filename() # filename prefix
  stateyear = Scalar(int)
  statemonth = Scalar(int)
  stateday = Scalar(int)
  state_format = Scalar(str) # netcdf|ascii
  grid_decimal = Scalar(int)
  wind_h = Scalar(int)
  measure_h = Scalar(int)
  alma_input = Boolean()
  forcing1 = Filename()
  force_format = Scalar(str) # netcdf|ascii
  force_endian = Scalar(str) # little|big
  n_types = Scalar(int)
  force_type = Mapping()
  force_dt = Scalar(int)
  forceyear = Scalar(int)
  forcemonth = Scalar(int)
  forceday = Scalar(int)
  forcehour = Scalar(int)
  nlayer = Scalar(int)
  nodes = Scalar(int)
  soil = Filename()
  baseflow = Scalar(str)
  arc_soil = Boolean()
  vegparam = Filename()
  vegparam_lai = Boolean()
  lai_src = Scalar(str)
  veglib = Filename()
  root_zones = Scalar(int)
  snow_band = Scalar(str) # int, str
  result_dir = Filename() # dirname
  out_step = Scalar(int)
  skipyear = Scalar(int)
  compress = Boolean()
  output_format = Scalar(str) # netcdf|ascii
  alma_output = Boolean()
  prt_header = Boolean()
  prt_snow_band = Boolean()
  netcdf_attribute = Mapping()
  n_outfiles = Scalar(int)
  outfiles = OutfileList()

  # FIXME: Rewrite all of these *date setter/getters to all use the same logic
  @property
  def startdate(self):
    return date(self.startyear, self.startmonth, self.startday)

  @startdate.setter
  def startdate(self, value):
    self.startyear, self.startmonth, self.startday = value.year, value.month, value.day

  @property
  def enddate(self):
    return date(self.endyear, self.endmonth, self.endday)

  @enddate.setter
  def enddate(self, value):
    self.endyear, self.endmonth, self.endday = value.year, value.month, value.day

  @property
  def statedate(self):
    return date(self.stateyear, self.statemonth, self.stateday)

  @statedate.setter
  def statedate(self, value):
    self.stateyear, self.statemonth, self.stateday = value.year, value.month, value.day

  @property
  def glacier_accum_startdate(self):
    return date(self.glacier_accum_start_year,
          self.glacier_accum_start_month,
          self.glacier_accum_start_day)

  @glacier_accum_startdate.setter
  def glacier_accum_startdate(self, value):
    self.glacier_accum_startyear, self.glacier_accum_startmonth,\
      self.glacier_accum_startday = value.year, value.month, value.day

  def __init__(self, input_stream):
    for idx, line in enumerate(input_stream):
      if line.isspace() or line.startswith('#'):
        continue

      name, value = line.split(None, 1)
      name = name.lower()
      value = value.strip()

      # This set of keys are "special". Use a different key that knows
      # how to handle them. This a hack. Don't ever do something like
      # this again.
      if name in ('outfile', 'outvar'):
        name = 'outfiles'

      setattr(self, name, value)

  def _get_descriptor(self, name):
    ''' descriptors have getters that return their contents rather than the \
      object itself. Use this if you *really* need the object itself (only \
      if you know what you are doing!)
    '''
    return self.__class__.__dict__[name]

  def _str_member(self, member):
    ''' Returns the string representation of an member as it would be written \
      in the VIC global file
    '''
    cls = self.__class__
    return cls.__dict__[member].__str__(self, cls, member)

  def __str__(self):
    return ''.join([ self._str_member(member) for member in self.member_order ])

  def write(self, filename):
    with open(filename, 'w') as f:
      f.write(str(self))
