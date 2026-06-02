import numpy as np
import numpy_compat  # noqa: F401  # NumPy 2.0 vs Landlab
import matplotlib.pyplot as plt

from landlab import RasterModelGrid, imshow_grid, imshowhs_grid # landlab grid components
from landlab.io import read_esri_ascii

def ss_fault(grid, fault_loc_y, slip_rate, method, accumulate):

    #grid: raster model grid
    #fault_loc_y: y value (row) where fault is located
    #total_slip: total slip in meters
    #slip_rate: slip rate in mm/yr taken from input file
    #method: drop or roll
    #accumulate: accumulated slip for this event (in meters), used to determine number of columns to shift.
    #faulted surface faulted_fields=[['topographic__elevation'], ['bedrock__elevation'], ['soil__depth']]

    nrows= grid.number_of_node_rows
    ncols= grid.number_of_node_columns

    z_original= grid.at_node["topographic__elevation"]
    z_original_reshaped = np.reshape(z_original, (nrows, ncols))
    bed_original = grid.at_node["bedrock__elevation"]
    bed_original_reshaped = np.reshape(bed_original, (nrows, ncols))
    soil_original = grid.at_node["soil__depth"]
    soil_original_reshaped = np.reshape(soil_original, (nrows, ncols))


    # INPUT FROM USER ABOUT TECTONICS
    #slip_rate = (total_slip / total_time)  # in m/yr
    print(' the slip rate of your fault is ' + str(slip_rate) + ' in mm/yr')
    number_cols = int(accumulate/grid.dx)   # in m assuming characteristic behavior #number of columns
    print('number of columns dropped per event is: ' + str(number_cols))

    if method == 'drop':
        # Drop/replicate on lower domain (footwall) interior columns; upper (hanging wall) fixed.
        # Right lateral: remove trailing interior columns, pad at left with duplicated first column of shortened strip.
        # Interior only (1:-1) matches roll boundary handling.

        index_columns_to_be_deleted = -1 * number_cols
        hang_z = z_original_reshaped[fault_loc_y:, :]
        foot_z = z_original_reshaped[:fault_loc_y, :]
        foot_in_z = foot_z[:, 1:-1]
        foot_drop_z = np.delete(foot_in_z, np.s_[index_columns_to_be_deleted:], axis=1)
        filling_z = np.array([foot_drop_z[:, 0], ] * number_cols).transpose()
        foot_z[:, 1:-1] = np.hstack((filling_z, foot_drop_z))
        new_z = np.vstack((foot_z, hang_z))
        z_reshaped_after_shift = np.reshape(new_z, (nrows * ncols))
        z_original[:] = z_reshaped_after_shift

        hang_s = soil_original_reshaped[fault_loc_y:, :]
        foot_s = soil_original_reshaped[:fault_loc_y, :]
        foot_in_s = foot_s[:, 1:-1]
        foot_drop_s = np.delete(foot_in_s, np.s_[index_columns_to_be_deleted:], axis=1)
        filling_s = np.array([foot_drop_s[:, 0], ] * number_cols).transpose()
        foot_s[:, 1:-1] = np.hstack((filling_s, foot_drop_s))
        new_soil = np.vstack((foot_s, hang_s))
        soil_reshaped_after_shift = np.reshape(new_soil, (nrows * ncols))
        soil_original[:] = soil_reshaped_after_shift

        hang_b = bed_original_reshaped[fault_loc_y:, :]
        foot_b = bed_original_reshaped[:fault_loc_y, :]
        foot_in_b = foot_b[:, 1:-1]
        foot_drop_b = np.delete(foot_in_b, np.s_[index_columns_to_be_deleted:], axis=1)
        filling_b = np.array([foot_drop_b[:, 0], ] * number_cols).transpose()
        foot_b[:, 1:-1] = np.hstack((filling_b, foot_drop_b))
        new_bed = np.vstack((foot_b, hang_b))
        bed_reshaped_after_shift = np.reshape(new_bed, (nrows * ncols))
        bed_original[:] = bed_reshaped_after_shift

        # imshow_grid(grid, z_original, cmap='viridis')
        # # imshow_grid(grid, bed_original, cmap='gray')
        # # imshow_grid(grid, soil_original, cmap='viridis')
        # # fig = plt.figure(figsize=[16,10])
        # plt.show()

    if method == 'roll':
        # Roll lower domain (footwall) along strike; upper domain (hanging wall) stays fixed.
        # Interior columns only (1:-1) match prior boundary handling.

        z_bottom_new = np.roll(z_original_reshaped[:fault_loc_y, 1:-1], number_cols, axis=1)
        z_original_reshaped[:fault_loc_y, 1:-1] = z_bottom_new
        z_reshaped_after_shift = np.reshape(z_original_reshaped, (nrows * ncols))
        z_original[:] = z_reshaped_after_shift

        soil_bottom_new = np.roll(soil_original_reshaped[:fault_loc_y, 1:-1], number_cols, axis=1)
        soil_original_reshaped[:fault_loc_y, 1:-1] = soil_bottom_new
        soil_reshaped_after_shift = np.reshape(soil_original_reshaped, (nrows * ncols))
        soil_original[:] = soil_reshaped_after_shift

        bed_bottom_new = np.roll(bed_original_reshaped[:fault_loc_y, 1:-1], number_cols, axis=1)
        bed_original_reshaped[:fault_loc_y, 1:-1] = bed_bottom_new
        bed_reshaped_after_shift = np.reshape(bed_original_reshaped, (nrows * ncols))
        bed_original[:] = bed_reshaped_after_shift

        # imshow_grid(grid, z_original, cmap='viridis')
        # plt.show()
    return grid 
    #displacement += slip_per_event

    # imshow_grid(grid, z_original, cmap='viridis')
    # plt.show()


# This is just an example


# #READING ORIGINAL TOPO
# (grid,z_original)=read_esri_ascii('/Users/taranguiz/Research/CSDMS_summer_2022/output_new_topo_ddd_5/finaltopo_topographic__elevation.asc', name="topographic__elevation")
# grid.set_closed_boundaries_at_grid_edges(bottom_is_closed=False, left_is_closed=True, right_is_closed=True, top_is_closed=True)
# imshow_grid(grid,z_original, cmap='viridis')
# plt.show()
#
# grid.add_zeros("node", "soil__depth") #add field to the grid
# grid.at_node["soil__depth"]=grid.at_node["soil__depth"]+2 #2
# grid.at_node["bedrock__elevation"]=grid.at_node["topographic__elevation"] - grid.at_node["soil__depth"]
#
# soil=grid.at_node['soil__depth']
# bed= grid.at_node['bedrock__elevation']
#
# imshow_grid(grid, soil)
# plt.show()
# displacement = 0
# time_track =0
# model_time=100000
# dt=100
# total_slip=400
# desired_slip_per_event=(total_slip/model_time)*dt
# print(desired_slip_per_event)
# # effective_slip= grid.dx/desired_slip_per_event
# accumulate=0
#     # np.copy(desired_slip_per_event)
#
# while time_track < model_time:
#
#     accumulate+= desired_slip_per_event
#     time_track += dt
#     print('is accumulating')
#
#     if accumulate >= grid.dx:
#         ss_fault(grid=grid, fault_loc_y=int(grid.number_of_node_rows/2.),  total_slip=total_slip, total_time=model_time, method='drop', accumulate=accumulate)
#         accumulate= accumulate%grid.dx
#         print('one slip')
#         imshow_grid(grid, z_original, cmap='viridis', shrink=0.5)
#         plt.title('Dropping method')
#         plt.show()
#         plt.close()
