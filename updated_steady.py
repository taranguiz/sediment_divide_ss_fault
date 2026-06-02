#!/usr/bin/env python3
# import time
import numpy as np
import numpy_compat  # noqa: F401  # NumPy 2.0 vs Landlab: must run before landlab
import matplotlib.pyplot as plt
import pickle
import os
from collections import defaultdict
from copy import deepcopy

#from Landlab
from landlab import RasterModelGrid, imshow_grid, imshowhs_grid
from landlab.io import read_esri_ascii
from landlab.io.netcdf import write_netcdf
from landlab.io.netcdf import write_raster_netcdf
from landlab.io.native_landlab import save_grid

#Hillslope geomorphology
from landlab.components import DepthDependentTaylorDiffuser
# from landlab.components import DepthDependentDiffuser

#Fluvial Geomorphology and Flow routing
from landlab.components import FlowDirectorMFD #trying the FlowDirectorMFD
from landlab.components import FlowAccumulator, Space, FastscapeEroder, PriorityFloodFlowRouter
from landlab.components.space import SpaceLargeScaleEroder
from landlab_compat import exponential_weatherer
from util import get_file_sequence, add_file_to_writer, save_grid_state

def init_grid(config):
    print("initializing grid")
    #Instantiate model grid
    mg= RasterModelGrid ((config.nrows,config.ncols), config.dxy)
    #add field topographic elevation.
    mg.add_zeros("node", "topographic__elevation")

    np.random.seed(seed=5000)
    #creating initial model topography
    random_noise = (np.random.rand(len(mg.node_y)))

    #add the topo to the field
    mg["node"]["topographic__elevation"] +=random_noise

    # add field 'soil__depth' to the grid
    mg.add_zeros("node", "soil__depth", clobber=True)

    # Set of initial soil depth at core nodes
    mg.at_node["soil__depth"][mg.core_nodes] = config.H0  # meters

    # Add field 'bedrock__elevation' to the grid
    mg.add_zeros("bedrock__elevation", at="node")
    # Sum 'soil__depth' and 'bedrock__elevation'
    # to yield 'topographic elevation'
    mg.at_node["bedrock__elevation"][:] = mg.at_node["topographic__elevation"]
    mg.at_node["topographic__elevation"][:] += mg.at_node["soil__depth"]
    mg.add_zeros("node", "soil_production__rate", clobber=True)
    # soil_production_rate= mg.at_node["soil_production__rate"]
    mg.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=False,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True
    )
    return mg


def build_steady_topo(config, writer):
    print("building steady topo")

    mg = init_grid(config)
    # Initialize list to store all states
    grid_states = {}

    rock = mg.at_node["bedrock__elevation"]
    soil = mg.at_node["soil__depth"]
    z = mg.at_node["topographic__elevation"]

    # instantiate components
    #Weathering
    
    #if the new version of landlab is only one _, if the old one is __
    expweath = exponential_weatherer(mg, config.P0, config.Hstar)


    # Hillslope with Taylor Diffuser
    ddtd=DepthDependentTaylorDiffuser(mg,slope_crit=config.Sc,
                                    soil_transport_velocity=config.V0,
                                    soil_transport_decay_depth=config.Hstar,
                                    nterms=2,
                                    dynamic_dt=True,
                                    if_unstable='warn', 
                                    courant_factor=0.1)

    #Flow Router
    fr=PriorityFloodFlowRouter(mg,
                            flow_metric='D8',
                            separate_hill_flow=False, #this is the default changing to test the effect of separate hill flow
                            hill_flow_metric="Quinn",
                            update_hill_flow_instantaneous=False, #this is the default changing to test the effect of separate hill flow
                            suppress_out=True, runoff_rate=config.run_off)                        
    #SPACE Large Scale
    space= SpaceLargeScaleEroder(mg,
                                K_sed=config.K_sed,
                                K_br=config.K_br,
                                F_f=config.F_f,
                                phi=config.phi,
                                H_star=config.H_star,
                                v_s=config.Vs,
                                m_sp=config.m_sp,
                                n_sp=config.n_sp,
                                sp_crit_sed=config.sp_crit_sed,
                                sp_crit_br=config.sp_crit_br)
    #tracking time
    time = 0

    while time < config.total_steady_time:

        rock[mg.core_nodes] += config.uplift_rate * config.dt_steady
        z[mg.core_nodes] += config.uplift_rate * config.dt_steady
        # z[:] = rock + soil
        
        #hillslope
        # expweath.run_one_step()
        expweath.calc_soil_prod_rate()
        ddtd.run_one_step(config.dt_steady)

        #run flow router
        fr.run_one_step()

        #run space
        space.run_one_step(config.dt_steady)
       

        if time%10000 == 0:
            # # Save grid state
            # state = save_grid_state(mg, time)
            # grid_states[time] = deepcopy(mg)
            # print(f"Saved grid state at time {time}")
            
            # Also save visualization
            imshow_grid(mg, z, cmap='coolwarm', shrink=0.5, grid_units=['m', 'm'])
            plt.title('Topography after ' + str(time) + ' years')
            loop_topo_img = f'{config.home_path}/{config.save_location}/{config.model_name}{get_file_sequence(time, config)}.png'
            plt.savefig(
                loop_topo_img, 
                dpi=300, facecolor='white'
            )
            add_file_to_writer(writer, loop_topo_img)
            plt.clf()
        if time%1000000 == 0:
            # Save grid state
            state = save_grid_state(mg, time)
            grid_states[time] = deepcopy(mg)
            print(f"Saved grid state at time {time}")
        
        print(time)
        time = time + config.dt_steady

    # Save final state into the in-memory dictionary
    final_state = save_grid_state(mg, config.total_steady_time)
    grid_states[config.total_steady_time] = deepcopy(mg)

    # Save all states (time -> grid) to a pickle file in the model's output folder
    output_filename = f'{config.home_path}/{config.save_location}/{config.model_name}_grid_states.pkl'
    with open(output_filename, 'wb') as f:
        pickle.dump(grid_states, f)
    
    print(f"\nSaved {len(grid_states)} grid states to {output_filename}")

    # Also save the final steady-state grid as a standalone pickle for faulting runs
    steady_dir = os.path.join(config.home_path, "output", "steady_state_files")
    os.makedirs(steady_dir, exist_ok=True)
    final_grid_path = os.path.join(steady_dir, f"final_state_{config.model_name}.pkl")
    with open(final_grid_path, "wb") as f:
        pickle.dump(deepcopy(mg), f)
    print(f"Saved final steady grid to {final_grid_path}")

    return mg  # Return the final grid instead of the metrics

# build_steady_topo()

# #fig = plt.figure(figsize=[8, 8])
# imshow_grid(mg, z, cmap='terrain', grid_units=['m', 'm'])
# plt.title('Topography after ' + str(int((tmax))) + ' years')
# #plt.show()        
# #imshow_grid(mg, z, cmap='terrain', grid_units=['m', 'm'])
# final_topo_img= f'/home/jupyter-taranguiz/StrikeSlip/steady/output_topo/{model_name}/final.png'
# plt.savefig(final_topo_img, dpi=300, facecolor='white')
# add_file_to_writer(final_topo_img)
# writer.close()

# print(type(Mean_elev))
# print(Mean_elev.size)
# print(type(model_time))
# print(model_time.size)

# fig2= plt.figure(figsize=[8,8])
# plt.plot(model_time, Mean_elev)
# plt.xlabel('model iterations')
# plt.ylabel('mean elevation[m]')
# plt.savefig(f'/home/jupyter-taranguiz/StrikeSlip/steady/output_topo/{model_name}/mean_elevation.jpg')

# fig3= plt.figure(figsize=[8,8])
# plt.plot(model_time, Mean_da)
# plt.xlabel('time [years]')
# plt.ylabel('mean drainage area [m2]')
# plt.savefig(f'/home/jupyter-taranguiz/StrikeSlip/steady/output_topo/{model_name}/mean_drainage.jpg')

# fig4= plt.figure(figsize=[8,8])
# plt.plot(model_time, Mean_soil)
# plt.xlabel('time [years]')
# plt.ylabel('mean soil_depth [m]')
# plt.savefig(f'/home/jupyter-taranguiz/StrikeSlip/steady/output_topo/{model_name}/mean_soil.jpg')

# np.savetxt(f'/home/jupyter-taranguiz/StrikeSlip/steady/output_topo/{model_name}/{model_name}mean_elev.txt',
#            (Mean_elev),
#            delimiter=',',
#            header='Mean_elev',
#            comments='')
# np.savetxt(f'/home/jupyter-taranguiz/StrikeSlip/steady/output_topo/{model_name}/{model_name}mean_da.txt',
#            (Mean_da),
#            delimiter=',',
#            header='Mean_da',
#            comments='')
# np.savetxt(f'/home/jupyter-taranguiz/StrikeSlip/steady/output_topo/{model_name}/{model_name}mean_soil.txt',
#            (Mean_soil),
#            delimiter=',',
#            header='Mean_soil',
#            comments='')

# print(mg.fields())
# try: 
#     write_netcdf(
#         f'/home/jupyter-taranguiz/StrikeSlip/steady/output_topo/{model_name}/steady-state.nc', 
#         mg, 
#         format='NETCDF4', 
#         names=[
#             'bedrock__elevation',
#             'drainage_area',
#             'flood_status_code',
#             'flow__link_to_receiver_node',
#             'flow__receiver_node',
#             'flow__receiver_proportions',
#             'flow__upstream_node_order',
#             'soil__depth',
#             'soil_production__rate',
#             'surface_water__discharge',
#             'topographic__elevation',
#             'topographic__steepest_slope',
#             'water__unit_flux_in'
#         ]
#     )
# except Exception as e:
#     print(str(e))

# try: 
#     mg.save(f'/home/jupyter-taranguiz/StrikeSlip/steady/output_topo/{model_name}/steady-state.nc')
# except Exception as e:
#     print(str(e))

#     # if float(diff) < float(uplift_rate):
#     #     print('steady state reached')
#     #     break


