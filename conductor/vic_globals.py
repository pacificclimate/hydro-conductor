''' This module represents the VIC "global parameters file", essentially
    all of the settings that affect VIC at run time
'''

__all__ = ['get_global_parms', 'update_global_parms', 'write_global_parms_file']

from collections import OrderedDict

def get_global_parms(global_parm_file):
    """ Parses the initial VIC global parameters file created by the user with the settings for the entire VIC-RGM run """
    global_parms = OrderedDefaultdict()
    n_outfile_lines = 0
    with open(global_parm_file, 'r') as f:
        for line in f:
            #print('line: {}'.format(line))
            if not (line.isspace() or line.startswith('#')):
                split_line = line.split()
                #print('columns: {}'.format(split_line))
                parm_name = split_line[0]
                if parm_name == 'OUTFILE': # special case because there are multiple occurrences, not consecutive
                        n_outfile_lines += 1
                        parm_name = 'OUTFILE_' + str(n_outfile_lines)
                elif parm_name == 'OUTVAR': # special case because multiple OUTVAR lines follow each OUTFILE line
                        parm_name = 'OUTVAR_' + str(n_outfile_lines)
                try:
                    if global_parms[parm_name]: # if we've already read one or more entries of this parm_name
#                        print('parm {} exists already'.format(parm_name))
                        global_parms[parm_name].append(split_line[1:])
                except:
                    global_parms[parm_name] = []
                    global_parms[parm_name].append(split_line[1:])
                    #print('global_parms[{}]: {}'.format(parm_name,global_parms[parm_name]))
                # We need to create a placeholder in this position for INIT_STATE if it doesn't exist in the initial
                # global parameters file, to be used for all iterations after the first year
                if parm_name == 'OUTPUT_FORCE': # OUTPUT_FORCE should always immediately precede INIT_STATE in the global file
                    global_parms['INIT_STATE'] = []
    return global_parms

def update_global_parms(global_parms, temp_vpf, temp_snb, num_snow_bands, start_date, end_date, init_state_file, state_date):
    """ Updates the global_parms dict at the beginning of an annual iteration """
    # Set start and end dates for the upcoming VIC run (should be one year long, except for the spin-up run)
    global_parms['STARTYEAR'] = [[str(start_date.year)]]
    global_parms['STARTMONTH'] = [[str(start_date.month)]]
    global_parms['STARTDAY'] = [[str(start_date.day)]]
    global_parms['ENDYEAR'] = [[str(end_date.year)]]
    global_parms['ENDMONTH'] = [[str(end_date.month)]]
    global_parms['ENDDAY'] = [[str(end_date.day)]]
    # Set new output state file parameters for upcoming VIC run
    global_parms['STATEYEAR'] = [[str(state_date.year)]]
    global_parms['STATEMONTH'] = [[str(state_date.month)]]
    global_parms['STATEDAY'] = [[str(state_date.day)]]
    global_parms['VEGPARAM'] = [[temp_vpf]]
    global_parms['SNOW_BAND'] = [[num_snow_bands, temp_snb]]
    # All VIC iterations except for the initial spin-up period have to load a saved state from the previous
    if init_state_file:
        global_parms['INIT_STATE'] = [[init_state_file]]
        
def write_global_parms_file(global_parms, temp_gpf):
    """ Reads existing global_parms OrderedDict and writes out a new temporary VIC Global Parameter File 
        temp_gpf for feeding into VIC 
    """
    with open(temp_gpf, 'w') as f:
        writer = csv.writer(f, delimiter=' ')
        for parm_name, parm_value in global_parms.items():
            #print('write_global_parms_file: parm_name: {} parm_value: {}'.format(parm_name, parm_value))
            num_parm_lines = len(global_parms[parm_name])
            if parm_name == 'INIT_STATE' and not global_parms['INIT_STATE']:
                pass
            elif parm_name[0:8] == 'OUTFILE_':
                line = []
                line.append('OUTFILE')
                for value in parm_value[0]:
                    line.append(value)
                writer.writerow(line)
            elif parm_name[0:7] == 'OUTVAR_':
                for line_num in range(num_parm_lines):
                    line = []
                    line.append('OUTVAR')
                    for value in parm_value[line_num]:
                        line.append(value)
                        writer.writerow(line)
            elif num_parm_lines == 1:
                line = []
                line.append(parm_name)
                for value in parm_value[0]:
                    line.append(value)
                writer.writerow(line)
            elif num_parm_lines > 1:
                for line_num in range(num_parm_lines):
                    line = []
                    line.append(parm_name)
                    for value in parm_value[line_num]:
                        line.append(value)
                    writer.writerow(line)
