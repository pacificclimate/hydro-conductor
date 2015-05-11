'''vegparams.py

   This module provides a single class for reading/writing the "Vegetation
   Parameter File", a delightfully custom data serialization format.
   One can use the module as such:

   >>> vp = VegParams(my_file)

   Then the following attributes will be accessible:

   >>> vp.cell_ids, vp.num_veg_tiles, vp.cells
'''

from collections import OrderedDict
import csv

__all__ = ['VegParams']

def leaves(d):
    '''Return the number of leaves in a multilevel dict
    '''
    if hasattr(d, 'values'):
        return sum([leaves(value) for value in d.values()])
    elif hasattr(d, '__len__'):
        return len(d)
    else:
        return 1
                    
class VegParams(object):
    glacier_id = '22'
    bare_soil_id ='19'

    def __init__(self, veg_parm_file=None):
        self.cells = OrderedDict()
        self.residual_area_fracs = {}
        if veg_parm_file:
            self.load(veg_parm_file)

    @property
    def cell_ids(self):
        '''A list of all cell ids in the vegetation parameter file
        '''
        return list(self.cells.keys())

    @property
    def num_veg_tiles(self):
        '''dictionary reporting the total number vegetation tiles over each
           cell. This number is the sum of all vegetation tiles in all bands.
        '''
        # Count is two levels deep
        return OrderedDict( (cell_id, leaves(cell)) for cell_id, cell in self.cells.items() )

    @staticmethod
    def read_one_cell(f):
        '''read all data from one cell and advance the file pointer to the next
        '''
        try:
            cell_id, num_veg = f.readline().split()
        except ValueError:
            return None
        cell = OrderedDict()
        for _ in range(int(num_veg)):
            line = f.readline()
            band_id = line.split()[-1]
            try:
                cell[band_id].append(line.split())
            except KeyError:
                cell[band_id] = [ (line.split()) ]
        return cell_id, cell

    def load(self, filename):
        '''iterate over each block of cells in the file and
           populate the class dict
        '''
        with open(filename, 'r') as f:
            while True:
                try:
                    id_, cell = self.read_one_cell(f)
                except TypeError:
                    return
                else:
                    self.cells[id_] = cell    

    def save(self, filename):
        ''' write the vegetation parameters out to a file of the same 
            format as the original vegetation parameters file
        '''
        print('VegParams.save({})...'.format(filename))
        with open(filename, 'w') as f:
            writer = csv.writer(f, delimiter=' ')
            for cell in self.cells:
                writer.writerow([cell, self.num_veg_tiles[cell]])
                for band in self.cells[cell]:
                    for line in self.cells[cell][band]:
                        writer.writerow(line)
                        print(' '.join(map(str, line)))

   def init_residual_area_fracs(self, snb_parms):
        """ Reads the initial snow band area fractions and glacier vegetation (HRU) 
            tile area fractions and calculates the initial residual area fractions 
        """
        print('VegParams.init_residual_area_fracs()...')
        for cell in self.cells:
            print('cell: {}'.format(cell))
            self.residual_area_fracs[cell] = {}
            for band in self.cells[cell]:
                print('band: {}'.format(band))
                glacier_exists = False
                for line_idx, line in enumerate(self.cells[cell][band]):
                    print('line: {}'.format(line))
                    if line[0] == self.glacier_id:
                        glacier_exists = True
                        self.residual_area_fracs[cell][band] = snb_parms[cell][0][int(band)] - float(self.cells[cell][band][line_idx][1])
                        print('Glacier exists in this band. residual_area_fracs[{}][{}] = {} - {} = {}'.format(cell, band, snb_parms[cell][0][int(band)], float(self.cells[cell][band][line_idx][1]), self.residual_area_fracs[cell][band]))
                        if self.residual_area_fracs[cell][band] < 0:
                            print('init_residual_area_fracs(): Error: Calculated a negative residual area fraction for cell {}, band {}. The sum of vegetation tile fraction areas for a given band in the Vegetation Parameter File must be equal to the area fraction for that band in the Snow Band File. Exiting.\n'.format(cell, band))
                            sys.exit(0)
                        break
                if not glacier_exists:
                    self.residual_area_fracs[cell][band] = snb_parms[cell][0][int(band)]
                    print('No glacier in this band.  residual_area_fracs[{}][{}] = {}'.format(cell, band, self.residual_area_fracs[cell][band]))

    def update(self, area_frac_bands, area_frac_glacier):
        """ Updates vegetation parameters for all VIC grid cells by applying calculated changes 
            in glacier area fractions across all elevation bands 
        """
        print('VegParams.update()...')
        for cell in self.cells:
            print('cell: {}'.format(cell))
            for band in self.cells[cell]:
                if area_frac_bands[cell][int(band)] > 0: # If we (still) have anything in this elevation band
                    new_residual_area_frac = area_frac_bands[cell][int(band)] - area_frac_glacier[cell][int(band)]
                    delta_residual = self.residual_area_fracs[cell][band] - new_residual_area_frac
                    veg_types = [] # identify and temporarily store vegetation types (HRUs) currently existing within this band
                    bare_soil_exists = False # temporary boolean to identify if bare soil vegetation type (HRU) currently exists within this band
                    glacier_exists = False # temporary boolean to identify if glacier vegetation type (HRU) currently exists within this band 
                    sum_previous_band_area_fracs = 0 # sum of band vegetation tile (HRU) area fractions in the last iteration
                    if delta_residual < 0: # glacier portion of this band has SHRUNK; update its area fraction and increase the bare soil component accordingly
                        print('\nGlacier portion of band {} has SHRUNK (delta_residual = {})'.format(band, delta_residual))
                        for line_idx, line in enumerate(self.cells[cell][band]): 
                            print('Current line in band {}:  {}'.format(band, line))
                            veg_types.append(int(line[0]))
                            print('Existing vegetation type identified in this elevation band: {}'.format(veg_types[-1]))
                            if veg_types[-1] == int(self.bare_soil_id):
                                bare_soil_exists = True
                                self.cells[cell][band][line_idx][1] = str(float(self.cells[cell][band][line_idx][1]) + abs(delta_residual))
                                print('Bare soil already exists in this band.  New area fraction: vp.cells[{}][{}][{}][1] = {}'.format(cell, band, line_idx, self.cells[cell][band][line_idx][1]))
                            elif veg_types[-1] == int(self.glacier_id):
                                self.cells[cell][band][line_idx][1] = str(area_frac_glacier[cell][int(band)])
                                print('Updated glacier area fraction: vp.cells[{}][{}][{}][1] = {}'.format(cell, band, line_idx, self.cells[cell][band][line_idx][1]))
                        if not bare_soil_exists: # insert new line for bare soil vegetation (HRU) type in the correct numerical position within this band
                            #new_line = self.bare_soil_id + ' ' + str(delta_residual) + ' ' + ''.join(str(x) for x in bare_soil_root_parms) + ' ' + str(band)
                            new_line = []
                            new_line.append(self.bare_soil_id)
                            new_line.append(str(delta_residual))
                            for num in bare_soil_root_parms:
                                new_line.append(str(num))
                            new_line.append(str(band))                        
                            bisect.insort_left(veg_types, int(self.bare_soil_id))
                            position = veg_types.index(int(self.bare_soil_id))
                            self.cells[cell][band].insert(position, new_line)
                            print('Bare soil did NOT already exist in this band. Inserted line: vp.cells[{}][{}][{}] = {}'.format(cell, band, position, self.cells[cell][band][position]))
                    elif delta_residual > 0: # glacier portion of this band has GROWN; decrease fraction of non-glacier HRUs by a total of delta_residual
                        print('\nGlacier portion of band {} has GROWN (delta_residual = {})'.format(band, delta_residual))
                        for line_idx, line in enumerate(self.cells[cell][band]): 
                            sum_previous_band_area_fracs += float(line[1])
                        print('sum_previous_band_area_fracs = {}'.format(sum_previous_band_area_fracs))
                        for line_idx, line in enumerate(self.cells[cell][band]):
                            print('Current line in band {}:  {}'.format(band, line))
                            veg_types.append(int(line[0])) # keep a record of vegetation types (HRUs) currently existing within this band
                            print('Existing vegetation type identified in this elevation band: {}'.format(veg_types[-1]))
                            if veg_types[-1] == int(self.glacier_id): # set glacier tile area fraction
                                glacier_exists = True
                                self.cells[cell][band][line_idx][1] = str(area_frac_glacier[cell][int(band)])
                                print('Glacier already exists in this band.  New glacier area fraction: vp.cells[{}][{}][{}][1] = {}'.format(cell, band, line_idx, self.cells[cell][band][line_idx][1]))
                            else:
                                # Get the change in area fraction for this vegetation type (HRU) based on its previous share of the residuals in this elevation band
                                print('Existing area fraction for this vegetation type: {}'.format(self.cells[cell][band][line_idx][1]))
                                delta_veg_area_frac = delta_residual * (float(self.cells[cell][band][line_idx][1]) / sum_previous_band_area_fracs)
                                print('Calculated delta_veg_area_frac = {}'.format(delta_veg_area_frac))
                                new_veg_area_frac = float(self.cells[cell][band][line_idx][1]) - delta_veg_area_frac
                                if new_veg_area_frac >= 0:
                                    self.cells[cell][band][line_idx][1] = str(new_veg_area_frac)
                                    print('Reduced area fraction of vegetation type {} by delta_veg_area_frac.  New area fraction: vp.cells[{}][{}][{}][1] = {}'.format(veg_types[-1], cell, band, line_idx, self.cells[cell][band][line_idx][1]))
                                else: #the existing vegetation tile (HRU) was overtaken by glacier
                                    # delete the existing vegetation tile (HRU) that was overtaken by glacier?
                                    #print('Area fraction of vegetation type {} in cell {}, band {} was reduced to {}. Removing tile (HRU) from vegetation parameters'.format(veg_types[-1], cell, band, new_veg_area_frac)
                                    #del self.cells[cell][band][line_idx]
                                    print('VegParams.update(): Error: Calculated a negative area fraction ({}) for vegetation type {} in cell {}, band {}. Exiting.\n'.format(new_veg_area_frac, veg_types[-1], cell, band))
                                    sys.exit(0)
                        if not glacier_exists: # insert new line for glacier vegetation (HRU) type in the correct numerical position within this band
                            #new_line = self.glacier_id + ' ' + str(area_frac_glacier[cell][int(band)]) + ' ' + ''.join(str(x) for x in glacier_root_parms) + ' ' + str(band)
                            new_line = []
                            new_line.append(self.glacier_id)
                            new_line.append(str(area_frac_glacier[cell][int(band)]))
                            for num in glacier_root_parms:
                                new_line.append(str(num))
                            new_line.append(str(band))
                            bisect.insort_left(veg_types, int(self.glacier_id))
                            position = veg_types.index(int(self.glacier_id))
                            self.cells[cell][band].insert(position, new_line)
                            print('Glacier did NOT already exist in this band. Inserted line with new glacier area fraction: vp.cells[{}][{}][{}] = {}'.format(cell, band, position, self.cells[cell][band][position]))
                else: # We have nothing left in this elevation band.  Set all fractional areas for all vegetation tiles (HRUs) in this band to zero
                    print('\nNothing left in elevation band {}.  Setting all fractional areas to zero.'.format(band))
                    for line_idx, line in enumerate(self.cells[cell][band]):
                        self.cells[cell][band][line_idx][1] = 0
                        new_residual_area_frac = 0
                # Set residual area fractions to the new calculated values for the next iteration
                self.residual_area_fracs[cell][band] = new_residual_area_frac

