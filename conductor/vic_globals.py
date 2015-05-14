''' This module represents the VIC "global parameters file", essentially
    all of the settings that affect VIC at run time
'''

__all__ = ['get_global_parms', 'update_global_parms', 'write_global_parms_file', 'Global']

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
            raise ValueError("Cannot set parameter to a file in a non-existant directory: {}".format(dirname(value)))
        super().__set__(instance, value)

class Mapping(object):
    def __init__(self):
        self.dict_ = {}
    def __set__(self, instance, value):
        try:
            key, value = value.split(None, 1)
        except:
            raise ValueError('Assigment to {} requires a key and value separated by whitespace'.format(self.__class__.__name__))
        self.dict_[key] = value
    def __get__(self, instance, cls):
        return self.dict_
    def __str__(self, instance, cls, name=''):
        return '\n'.join('{} {} {}'.format(name.upper(), k, v) for k, v in self.dict_.items())

class List(object):
    def __init__(self):
        self.value = []
    def __set__(self, instance, value):
        self.value.append(value)
    def __get__(self, instance, cls):
        return self.value
    def __str__(self, instance, cls, name=''):
        return '\n'.join('{} {}'.format(name.upper(), item) for item in self.value)

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
    
class Global(metaclass=OrderedMeta):
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
        self.glacier_accum_startyear, self.glacier_accum_startmonth, self.glacier_accum_startday = value.year, value.month, value.day

    def __init__(self, filename):
        with open(filename, 'r') as f:
            for line in f:
                if line.isspace() or line.startswith('#'):
                    continue

                name, value = line.split(None, 1)
                name = name.lower()
                value = value.strip()

                if name in ('outfile', 'outvar'):
                    name = 'outfiles'

                setattr(self, name, value)

    def _get_descriptor(self, name):
        ''' descriptors have getters that return their contents rather
            than the object itself. Use this if you *really* need the
            object itself (only if you know what you are doing!)
        '''
        return self.__class__.__dict__[name]

    def _str_member(self, member):
        ''' Get the string representation of an object 
        '''
        cls = self.__class__
        return cls.__dict__[member].__str__(self, cls, member)

    def __str__(self):
        return ''.join([ self._str_member(member) for member in self.member_order ])

# To have nested ordered defaultdicts
class OrderedDefaultdict(OrderedDict):
    # from: http://stackoverflow.com/questions/4126348/how-do-i-rewrite-this-function-to-implement-ordereddict/4127426#4127426
    def __init__(self, *args, **kwargs):
        if not args:
            self.default_factory = None
        else:
            if not (args[0] is None or isinstance(args[0], collections.Callable)):
                raise TypeError('first argument must be callable or None')
            self.default_factory = args[0]
            args = args[1:]
        super(OrderedDefaultdict, self).__init__(*args, **kwargs)
    def __missing__ (self, key):
        if self.default_factory is None:
            raise KeyError(key)
        self[key] = default = self.default_factory()
        return default
    def __reduce__(self):  # optional, for pickle support
        args = (self.default_factory,) if self.default_factory else ()
        return self.__class__, args, None, None, iter(self.items())

def get_global_parms(global_parm_file):
    """ Parses the initial VIC global parameters file created by the
        user with the settings for the entire VIC-RGM run
    """
    g = OrderedDefaultdict()
    n_outfile_lines = 0
    with open(global_parm_file, 'r') as f:
        for line in f:
            if line.isspace() or line.startswith('#'):
                continue

            name, value = line.split(None, 1)
            value = value.strip()

            # special case because there are multiple occurrences,
            # not consecutive
            if name == 'OUTFILE':
                n_outfile_lines += 1
                name = 'OUTFILE_' + str(n_outfile_lines)
            # special case because multiple OUTVAR lines follow each OUTFILE line
            elif name == 'OUTVAR':
                name = 'OUTVAR_' + str(n_outfile_lines)

            if name in g:
                # print('parm {} exists already'.format(name))
                if type(g[name]) != list:
                    g[name] = [ g[name] ]
                g[name].append(value)
            else:
                g[name] = value

            # We need to create a placeholder in this position for INIT_STATE if
            # it doesn't exist in the initial global parameters file, to be used
            # for all iterations after the first year
            #
            # OUTPUT_FORCE should always immediately precede INIT_STATE in the
            # global file
            if name == 'OUTPUT_FORCE': 
                g['INIT_STATE'] = []

    return g

def update_global_parms(global_parm_dict, temp_vpf, temp_snb, num_snow_bands,
                        start_date, end_date, init_state_file, state_date):
    """ Updates the global_parms dict at the beginning of an annual
        iteration
    """
    g = global_parm_dict
    # Set start and end dates for the upcoming VIC run (should be one year long, except for the spin-up run)
    g['STARTYEAR'] = [[str(start_date.year)]]
    g['STARTMONTH'] = [[str(start_date.month)]]
    g['STARTDAY'] = [[str(start_date.day)]]
    g['ENDYEAR'] = [[str(end_date.year)]]
    g['ENDMONTH'] = [[str(end_date.month)]]
    g['ENDDAY'] = [[str(end_date.day)]]
    # Set new output state file parameters for upcoming VIC run
    g['STATEYEAR'] = [[str(state_date.year)]]
    g['STATEMONTH'] = [[str(state_date.month)]]
    g['STATEDAY'] = [[str(state_date.day)]]
    g['VEGPARAM'] = [[temp_vpf]]
    g['SNOW_BAND'] = [[num_snow_bands, temp_snb]]
    # All VIC iterations except for the initial spin-up period have to load a saved state from the previous
    if init_state_file:
        g['INIT_STATE'] = [[init_state_file]]
        
def write_global_parms_file(global_parm_dict, temp_gpf):
    """ Reads existing global_parms OrderedDict and writes out a new temporary
        VIC Global Parameter File temp_gpf for feeding into VIC
    """
    g = global_parm_dict
    with open(temp_gpf, 'w') as f:
        for name, value in g.items():
            # Don't output this key if it's not set
            if name == 'INIT_STATE' and not g['INIT_STATE']:
                continue
            # We'll deal with these separately at the end
            if name.startswith('OUTFILE') or name.startswith('OUTVAR'):
                continue

            # The main event; either write out a list as multiple lines
            # or a scalar as one line
            if type(value) == list:
                for item in value:
                    f.write("{} {}\n".format(name, item))
            else:
                f.write("{} {}\n".format(name, value))

        # Special case for writing the hiearchical variables per output file
        for i in range(1, int(g['N_OUTFILES']) + 1):
            key = 'OUTFILE_{}'.format(i)
            f.write("OUTFILE {}\n".format(g[key]))
            for outvar in g['OUTVAR_{}'.format(i)]:
                f.write("OUTVAR {}\n".format(outvar))
