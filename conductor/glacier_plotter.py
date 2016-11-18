from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

class GlacierPlotter:
  def __init__(self, current_surf_dem, glacier_mask, bed_dem, date, \
    output_trace_files, temp_files_path, glacier_thickness_threshold):
    self.output_trace_files = output_trace_files
    self.temp_files_path = temp_files_path
    self.glacier_thickness_threshold = glacier_thickness_threshold

    self.fig = plt.figure(figsize=(16,12))
    self.surf_dem_plot_title = 'Surface DEM ' + date
    self.sub1 = self.fig.add_subplot(131)
    self.sub1.set_title(self.surf_dem_plot_title)
    self.sub1.set_xticks([])
    self.sub1.set_yticks([])
    self.sub1.get_axes().set_frame_on(True)
    self.img1 = self.sub1.imshow(current_surf_dem)

    self.sub1.invert_yaxis() # makes North up
    self.divider1 = make_axes_locatable(plt.gca())
    self.cax1 = self.divider1.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(self.img1, cax=self.cax1)

    self.glacier_mask_plot_title = 'Glacier Mask ' + date
    self.sub2 = self.fig.add_subplot(132)
    self.sub2.set_title(self.glacier_mask_plot_title)
    self.sub2.set_xticks([])
    self.sub2.set_yticks([])
    self.sub2.get_axes().set_frame_on(True)
    self.img2 = self.sub2.imshow(glacier_mask)
    self.sub2.invert_yaxis()
    self.divider2 = make_axes_locatable(plt.gca())
    self.cax2 = self.divider2.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(self.img2, cax=self.cax2)

    self.glacier_thickness_plot_title = 'Glacier Thickness ' + date
    self.sub3 = self.fig.add_subplot(133)
    self.sub3.set_title(self.glacier_thickness_plot_title)
    self.sub3.set_xticks([])
    self.sub3.set_yticks([])
    self.sub3.get_axes().set_frame_on(True)
    self.img3 = self.sub3.imshow( \
      current_surf_dem - glacier_thickness_threshold - bed_dem)
    self.sub3.invert_yaxis()
    self.divider3 = make_axes_locatable(plt.gca())
    self.cax3 = self.divider3.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(self.img3, cax=self.cax3)

    self.fig.tight_layout()
    plt.show(block=False)

    if self.output_trace_files:
      self.fig.savefig(temp_files_path + 'dem_+_glacier_mask_+_thickness_' + \
        date, dpi=self.fig.dpi)

  def update_plots(self, current_surf_dem, glacier_mask, bed_dem, date):
    self.surf_dem_plot_title = 'Surface DEM ' + date
    self.sub1.set_title(self.surf_dem_plot_title)
    self.img1.set_array(current_surf_dem)
    plt.show(block=False)
    plt.pause(0.01)

    self.glacier_mask_plot_title = 'Glacier Mask ' + date
    self.sub2.set_title(self.glacier_mask_plot_title)
    self.img2.set_array(glacier_mask)
    plt.show(block=False)
    plt.pause(0.01)

    self.glacier_thickness_plot_title = 'Glacier Thickness ' + date
    self.sub3.set_title(self.glacier_thickness_plot_title)
    self.img3.set_array(\
      current_surf_dem - glacier_thickness_threshold - bed_dem)
    plt.show(block=False)
    plt.pause(0.01)

    if self.output_trace_files:
      self.fig.savefig(self.temp_files_path + 'dem_+_glacier_mask_+_thickness_' + \
      date, dpi=self.fig.dpi)
