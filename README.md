# openCFS-3D-geometry

3d_ipopt_ramp_press-for-eqn
Optimizer used: ipopt
method: RAMP
Boundary Condition at Nozzle_curve: Pressure equation

running Terminal Commands:
1) cfs -m box3d-t_0.3-nx_200-ny_20-nz_20.mesh F3D_distributed_load
2) show_density F3D_distributed_load.density.xml
3) plotviz.py F3D_distributed_load.plot.dat -x 1 -y 2 4
4) For comparison between multiple .dat files: plotviz.py *.plot.dat -x 1 -y 2 4
