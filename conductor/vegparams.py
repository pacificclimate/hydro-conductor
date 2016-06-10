"""vegparams.py

   This module provides functions for reading/writing the VIC Vegetation
   Parameter File, a delightfully custom data serialization format.

"""

from collections import OrderedDict
import csv
from conductor.cells import Band, HydroResponseUnit

def read_one_cell(f):
  """Reads all data (elevation bands/hrus) for one cell and advance the
    file pointer to the next cell.
  """
  try:
    cell_id, num_veg = f.readline().split()
  except ValueError:
    return None
  hru_dict = {}
  for _ in range(int(num_veg)):
    line = f.readline()
    split_line = line.split()
    veg_type = int(split_line[0])
    area_frac = float(split_line[1])
    root_zone_parms = [ float(x) for x in split_line[2:8] ]
    band_id = int(split_line[8])
    key = (band_id, veg_type)
    hru_dict[key] = HydroResponseUnit(area_frac, root_zone_parms, band_id, veg_type)

  return cell_id, hru_dict

def load_veg_parms(filename):
  """ Reads in VIC vegetation parameter file and creates and partially
    initializes all VIC grid cells.
  """
  cells = OrderedDict()
  # iterate over each block of cells in the file and populate the class dict
  with open(filename, 'r') as f:
    while True:
      try:
        id_, cell = read_one_cell(f)
      except TypeError:
        break
      else:
        cells[id_] = cell
  return cells

def save_veg_parms(cells, filename):
  """ Writes the vegetation parameters out to a file of the same format as the
    original vegetation parameters file.
  """
  with open(filename, 'w') as f:
    writer = csv.writer(f, delimiter=' ')
    for cell_id, cell in cells.items():
      writer.writerow( [ cell_id, sum([band.num_hrus for band in cell.bands]) ] )
      for band_id, band in enumerate(cell.bands):
        for veg_type, hru in band.hrus.items():
          line = [veg_type, hru.area_frac] + hru.root_zone_parms + [ band_id ]
          writer.writerow(line)
