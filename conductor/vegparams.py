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

__all__ = ['VegParams', 'VegTypeParms', 'BandInfo']

class BandInfo(object):
    def __init__(self, area_frac_glacier, area_frac_open_ground, veg_type, area_frac, root_zone_parms ):
        self.area_frac_glacier = area_frac_glacier # this is per-band, and should have an initialization function, and gets updated during update()
        self.area_frac_open_ground = area_frac_open_ground # this is per-band, and should have an initialization function, and gets updated during update()
        self.veg_tile_parms = [ VegTileParms(veg_type, area_frac, root_zone_parms) ]

class VegTileParms(object):
    def __init__(self, veg_type, area_frac, root_zone_parms):
        self.veg_type = veg_type
        self.area_frac = area_frac
        self.root_zone_parms = root_zone_parms
                   
class VegParams(object):
    glacier_id = '22'
    bare_soil_id ='19'

    def __init__(self, veg_parm_file=None):
        self.cells = OrderedDict()
        self.area_frac_non_glacier = {}
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
        count = OrderedDict()
        for cell in self.cells:
            count[cell] = 0
            for band in self.cells[cell]:
                for tile in self.cells[cell][band].veg_tile_parms:
                    count[cell]+=1
        return count

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
            split_line = line.split()
            veg_type = split_line[0]
            area_frac = split_line[1]
            root_zone_parms = split_line[2:8]
            band_id = split_line[8]
            # TODO: if glacier_root_zone_parms or open_ground_zone_parms were provided at command line, 
            # override root_zone_parms here:
            #if veg_type is glacier_id and glacier_root_zone_parms:
            #    root_zone_parms = glacier_root_zone_parms
            #elif veg_type is open_ground and open_ground_zone_parms:
            #    root_zone_parms = open_ground_zone_parms
            try:
                cell[band_id].veg_tile_parms.append(VegTileParms(veg_type, area_frac, root_zone_parms))
            except KeyError:
                cell[band_id] = BandInfo(None, None, veg_type, area_frac, root_zone_parms)
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
                    for tile_idx, tile in enumerate(self.cells[cell][band].veg_tile_parms):
                        line = [ self.cells[cell][band].veg_tile_parms[tile_idx].veg_type, \
                            self.cells[cell][band].veg_tile_parms[tile_idx].area_frac ] \
                            + self.cells[cell][band].veg_tile_parms[tile_idx].root_zone_parms \
                            + [ band ]
                        writer.writerow(line)
                        print(' '.join(map(str, line)))

#NOTE: it seems we only need to init and track glacier_area_frac for each band, and can deduce non_glacier_area_frac at each update
    def init_non_glacier_area_fracs(self, snb_parms):
        """ Reads the initial snow band area fractions and glacier vegetation (HRU) 
            tile area fractions and calculates the initial non-glacier area fractions 
        """
        print('VegParams.init_non_glacier_area_fracs()...')
        for cell in self.cells:
            print('cell: {}'.format(cell))
            area_frac_non_glacier[cell] = {}
            for band, band_idx in enumerate(self.cells[cell]):
                print('band: {}'.format(band))
                glacier_exists = False
                for line_idx, line in enumerate(self.cells[cell][band]):
                    print('line: {}'.format(line))
                    if line[0] == self.glacier_id:
                        glacier_exists = True
                        area_frac_non_glacier[cell][band] = snb_parms.cells[cell].area_fracs[band_idx] - float(self.cells[cell][band][line_idx][1])
                        print('Glacier exists in this band. area_frac_non_glacier[{}][{}] = {} - {} = {}'\
                            .format(cell, band, snb_parms.cells[cell].area_fracs[band_idx], \
                                float(self.cells[cell][band][line_idx][1], area_frac_non_glacier[cell][band_idx])))
                        break
                    # add condition here for open ground and get initial A_open?
                    # if line[0] == self.open_ground_id:
                        #open_ground_exists = True
                        #area_frac_open_ground[cell][band] = snb_parms.cells[cell].area_fracs[band_idx] - float(self.cells[cell][band][line_idx][1])

                if not glacier_exists:
                    area_frac_non_glacier[cell][band] = snb_parms.cells[cell].area_fracs[band_idx]
                    print('No glacier in this band.  area_frac_non_glacier[{}][{}] = {}'.format(cell, band, area_frac_non_glacier[cell][band_idx]))
        return area_frac_non_glacier

    def update(self, snb_parms, old_area_frac_glacier, new_area_frac_glacier, area_frac_non_glacier):
        """ Updates vegetation parameters for all VIC grid cells by applying calculated changes 
            in glacier area fractions across all elevation bands 
        """
        print('VegParams.update()...')
        for cell in self.cells:
            print('cell: {}'.format(cell))
            for band, band_idx in enumerate(self.cells[cell]):
                if snb_parms.cells[cell].area_fracs[band_idx] > 0: # If we (still) have anything in this elevation band
                    if new_area_frac_glacier is not None:
                        # TODO: need to handle case where a band has been created since last iteration, 
                        # thus old_area_frac_glacier[cell][band_idx] will not exist.  
                        # This should probably be within VegParams class, with a create_band() method like SnbParams has
                        delta_glacier_area_frac = new_area_frac_glacier[cell][band_idx] - old_area_frac_glacier[cell][band_idx]
                    # Update if there was a change in glacier area, or if this is the initialization phase
                    if (delta_glacier_area_frac != 0) or (new_area_frac_glacier is None):
                        if new_area_frac_glacier is None: # initialization case
                            new_area_frac_glacier[cell][band_idx] = 0
                        
                        new_area_frac_non_glacier = snb_parms.cells[cell].area_fracs[band_idx] - new_area_frac_glacier[cell][band_idx]
                        residual_area_frac = new_area_frac_non_glacier - area_frac_non_glacier[cell][band_idx]
                        # Set area_frac_non_glacier[cell][band_idx] to new value for next iteration
                        area_frac_non_glacier[cell][band_idx] = new_area_frac_non_glacier

                        # Now need to apply changes to open ground area fraction


                    veg_types = [] # identify and temporarily store vegetation types (HRUs) currently existing within this band
                    bare_soil_exists = False # temporary boolean to identify if bare soil vegetation type (HRU) currently exists within this band
                    glacier_exists = False # temporary boolean to identify if glacier vegetation type (HRU) currently exists within this band 
                    #sum_previous_band_area_fracs = 0 # sum of band vegetation tile (HRU) area fractions in the last iteration
      # TODO: change this to look at delta area_frac_glacier instead, not delta_residual
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

