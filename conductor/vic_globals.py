''' This module represents the VIC "global parameters file", essentially
    all of the settings that affect VIC at run time
'''

__all__ = ['get_global_parms', 'update_global_parms', 'write_global_parms_file']

from collections import OrderedDict

# James's assessment of proper typing in the global file
# {'TIME_STEP': int,
#  'SNOW_STEP': int,
#  'STARTYEAR': int,
#  'STARTMONTH': int,
#  'STARTDAY': int,
#  'STARTHOUR': int,
#  'ENDYEAR': int,
#  'ENDMONTH': int,
#  'ENDDAY': int,
#  'FULL_ENERGY': bool,
#  'FROZEN_SOIL': bool,
#  'NO_FLUX': bool,
#  'DIST_PRCP': bool,
#  'CORRPREC': bool,
#  'MIN_WIND_SPEED': float,
#  'PREC_EXPT': float,
#  'GLACIER_ID': int,
#  'GLACIER_ACCUM_START_YEAR': int,
#  'GLACIER_ACCUM_START_MONTH': int,
#  'GLACIER_ACCUM_START_DAY': int,
#  'GLACIER_ACCUM_INTERVAL': int,
#  'OUTPUT_FORCE': bool,
#  'INIT_STATE': [], #??
#  'STATENAME': str # filename
#  'STATEYEAR': int,
#  'STATEMONTH': int,
#  'STATEDAY': int,
#  'STATE_FORMAT': str # NETCDF|??|??,
#  'GRID_DECIMAL': int #?,
#  'WIND_H': int # ?,
#  'MEASURE_H': int # ?,
#  'ALMA_INPUT': bool,
#  'FORCING1': str, # filename
#  'FORCE_FORMAT': str # NETCDF|??|??,
#  'FORCE_ENDIAN': str # LITTLE|BIG,
#  'N_TYPES': int,
#  'FORCE_TYPE': dict,
#  'FORCE_DT': int,
#  'FORCEYEAR': int,
#  'FORCEMONTH': int,
#  'FORCEDAY': int,
#  'FORCEHOUR': int,
#  'NLAYER': int,
#  'NODES': int,
#  'SOIL': str, # filename
#  'BASEFLOW': str,
#  'ARC_SOIL': bool,
#  'VEGPARAM': str, # filename
#  'VEGPARAM_LAI': bool,
#  'LAI_SRC': str,
#  'VEGLIB': str, # filename
#  'ROOT_ZONES': int,
#  'SNOW_BAND': (int, str ) #filename
#  'RESULT_DIR': str, # dirname
#  'OUT_STEP': int,
#  'SKIPYEAR': int,
#  'COMPRESS': bool,
#  'OUTPUT_FORMAT': str, # NETCDF|??
#  'ALMA_OUTPUT': bool,
#  'PRT_HEADER': bool,
#  'PRT_SNOW_BAND': bool,
#  'NETCDF_ATTRIBUTE': dict,
#  'N_OUTFILES': int,
#  'OUTFILE_?': # list of outfiles
#  'OUTFILE_N': (str (prefix), int)
#  'OUTVAR_N': list, # of strings
# }

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

            # special case because there are multiple occurrences, not consecutiev
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

            # We need to create a placeholder in this position for INIT_STATE if it doesn't exist in the initial
            # global parameters file, to be used for all iterations after the first year
            if name == 'OUTPUT_FORCE': # OUTPUT_FORCE should always immediately precede INIT_STATE in the global file
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
    """ Reads existing global_parms OrderedDict and writes out a new temporary VIC Global Parameter File
        temp_gpf for feeding into VIC
    """
    g = global_parm_dict
    with open(temp_gpf, 'w') as f:
        writer = csv.writer(f, delimiter=' ')
        for name, value in g.items():
            #print('write_g_file: name: {} value: {}'.format(name, value))
            num_parm_lines = len(g[name])
            if name == 'INIT_STATE' and not g['INIT_STATE']:
                pass
            elif name.startswith('OUTFILE_'):
                line = []
                line.append('OUTFILE')
                for value in value[0]:
                    line.append(value)
                writer.writerow(line)
            elif name.startswith('OUTVAR_'):
                for line_num in range(num_parm_lines):
                    line = []
                    line.append('OUTVAR')
                    for value in value[line_num]:
                        line.append(value)
                        writer.writerow(line)
            elif num_parm_lines == 1:
                line = []
                line.append(name)
                for value in value[0]:
                    line.append(value)
                writer.writerow(line)
            elif num_parm_lines > 1:
                for line_num in range(num_parm_lines):
                    line = []
                    line.append(name)
                    for value in value[line_num]:
                        line.append(value)
                    writer.writerow(line)
