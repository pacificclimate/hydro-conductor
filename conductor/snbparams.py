'''snbparams.py

   This module provides functions for reading/writing the VIC "Snow Band File",
   
   The format of the snow band file is one line per VIC cell:
   cell_id_0 area_frac_band_0 ... area_frac_band_N median_elev_band_0 ... median_elev_band_N
   (and optionally, Pfactor_band_0 ... Pfactor_band_N  although VIC no longer uses these)
   where N should be equal to num_snow_bands
'''

__all__ = ['load_snb_parms', 'save_snb_parms']

from collections import OrderedDict, deque
import csv

from conductor.cells import Band, HydroResponseUnit

def front_padding(l):
    ''' Count and return the number of falsy occurrances at the
        beginning of a list
    '''
    count = 0
    for x in l:
        if x:
            break
        else:
            count += 1
    return count


class PaddedDeque(deque):
    '''A Double Ended Queue (deque) of fixed length that is optionally
       None padded on both ends. Standard append operations consume
       unused padding until the deque is full. Standard pop() operations
       return the first non-None values and leave None-padding in their
       wake.
    '''
    def __init__(self, *args, left_padding=0):
        self.left_padding = left_padding
        super().__init__(*args)

    @property
    def right_padding(self):
        return self.maxlen - len(self) - self.left_padding

    def append(self, x):
        if len(self) == self.maxlen or self.right_padding == 0:
            raise IndexError("Cannot append to item to full PaddedDeque")
        super().append(x)

    def appendleft(self, x):
        if len(self) == self.maxlen or self.left_padding == 0:
            raise IndexError("Cannot append to item to full PaddedDeque")
        self.left_padding -= 1
        super().appendleft(x)

    def popleft(self):
        rv = super().popleft()
        self.left_padding += 1
        return rv

    def __getitem__(self, i):
        if i > self.maxlen - 1:
            raise IndexError("Index out of bounds for PaddedDeque")
        if i < self.left_padding or i >= (self.left_padding + super().__len__()):
            return None
        else:
            return super().__getitem__(i - self.left_padding)

    def padded_iter(self):
        yield from ([None] * self.left_padding)
        yield from super().__iter__()
        yield from ([None] * self.right_padding)

    def unpadded_enumerate(self):
        for index, item in enumerate(super().__iter__()):
            yield index + self.left_padding, item

    def peekleft(self):
        return super().__getitem__(0)

    def peekright(self):
        return super().__getitem__(-1)

    def __repr__(self):
        return '{}({}, {}, left_padding={})'.format(self.__class__.__name__,
                repr(list(super().__iter__())), self.maxlen, self.left_padding)


def load_snb_parms(snb_file, num_snow_bands):
    """ Reads in a Snow Band Parameter File and populates the median elevation
        property for each band withing an existing set of VIC cells. Creates a 
        band map to keep track of the lower bounds of each band (each spanning an
        elevation of band_size) and any zero pads provided by the user in the 
        Snow Band Parameter File (zero pads are required by VIC, to allow for glacier 
        growth/slide into previously non-existent elevations between iterations).
    """
    def assign_dummy_band_elevations(elevs):
        """ Replaces 0 pads in elevation list that is read in from the Snow Band 
            Parameter File with floor elevations for the dummy bands they are
            placeholders for.
        """
        left_pads = 0
        leftmost_floor = 0
        right_pads = 0
        rightmost_floor = 0
        for count1, elev in enumerate(elevs):
            if elev != 0:
                left_pads = count1
                leftmost_floor = elev - elev % Band.band_size
                break
        elevs.reverse()
        for count2, elev in enumerate(elevs):
            if elev != 0:
                right_pads = count2
                rightmost_floor = elev - elev % Band.band_size
                break
        elevs.reverse()
        left_fills = list(range((leftmost_floor - left_pads*Band.band_size), leftmost_floor, Band.band_size ))
        right_fills = list(range((rightmost_floor + Band.band_size), (rightmost_floor + right_pads*Band.band_size + Band.band_size), Band.band_size ))
        elevs[0:len(left_fills)] = left_fills
        elevs[elevs.index(0):] = right_fills
        return elevs

    with open(snb_file, 'r') as f:
        cells = OrderedDict()
        for line in f:
            split_line = line.split()
            cell_id = split_line[0]
            # Should have the cell_id followed by num_snow_bands columns 
            # for each of area fractions and median elevations 
            # (and NO Pfactor values, which are deprecated!)
            if len(split_line) != (num_snow_bands * 2 + 1):
                raise Exception(
                        'Number of columns ({}) in snow band file {} is '
                        'incorrect for the number of SNOW_BAND ({}) '
                        'given in the global parameter file (should be a '
                        '2 * SNOW_BAND, plus 1). Are you still including '
                        '(deprecated) Pfactor values? If so, remove them.'
                        .format(len(split_line), snb_file, num_snow_bands)
                )
            elevs = [ int(z) for z in split_line[num_snow_bands+1:] ]
            #left_padding = front_padding(elevs)
# NOTE: we are now creating placeholder Bands in the zero-pad positions (calculating the
# floor elevation for each and assigning it to the median_elev attribute)   
            #bands = [ Band(z) for z in elevs if z ]
# NOTE: we no longer need the PaddedDeque, as all possible band positions can be represented in a list
            #cell = PaddedDeque(bands, num_snow_bands, left_padding=left_padding)

            # Assign median (floor) elevations to 0-pad-derived bands
            elevs = assign_dummy_band_elevations(elevs)
            
            # Cell consists of a list of Bands (both valid and placeholders for potential Bands)
            cell = [ Band(z) for z in elevs ]

            cells[cell_id] = cell
    return cells

# TODO: update this...
def save_snb_parms(cells, filename, band_map):
    """ Assembles and writes updated snow band parameters to a new temporary
        Snow Band Parameter File for feeding back into VIC in the next iteration.
    """
    with open(filename, 'w') as f:
         writer = csv.writer(f, delimiter=' ')
         for cell_id, cell in cells.items():
            area_fracs = [ band.area_frac if map_value else 0 for map_value, band in zip(band_map, cell) ]
# TODO: only write out elevations for Bands that contain HRUs (otherwise they are placeholders and a 0 should be written)            
            elevations = [ band.median_elev if map_value else 0 for map_value, band in zip(band_map, cell) ]
            line = [cell_id] + area_fracs + elevations
            writer.writerow(line)
