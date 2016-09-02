
vic_full_path = '/home/mfischer/code/vic/vicNl'
rgm_full_path = '/home/mfischer/code/rgm/rgm'
output_path = '/home/mfischer/vic_dev/out/testing/'
temp_files_path = output_path + 'hydrocon_temp/'

# These are set in the VIC header snow.h, but should probably be passed in via the state file
MAX_SURFACE_SWE = 0.125
CH_ICE = 2100E+03

# Absolute tolerance under which HRU area fractions are considered zero
ZERO_AREA_FRAC_TOL = 0.00001

# Areas where snow/ice coverage exceeds GLACIER_MASK_THICKNESS_THRESHOLD
# (in units of meters) are considered to be glacier
GLACIER_THICKNESS_THRESHOLD = 2.0