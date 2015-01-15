[vtk]
|
|---- solution.pvd
|
|---- time_step_$(TIME_STEP).pvtu
|
|---- time_step_$(TIME_STEP)_$(PROCESS_RANK).vtu

For every simulation, ONE process generates "solution.pvd" (which referes to all "time_step_$(TIME_STEP).pvtu" files).
For every time step that shall be saved, ONE process generates the corresponding "time_step_$(TIME_STEP).pvtu" file (which refers to all corresponding "time_step_$(TIME_STEP)_$(PROCESS_RANK).vtu" files).
For every time step that shall be saved, EACH process generates one "time_step_$(TIME_STEP)_$(PROCESS_RANK).vtu" file (which contains the actual simulation data).

In order to visualize the sample data located in folder "vtk", start ParaView, open the file "solution.pvd", and press "Apply".
Select Filters->Alphabetical->Delaunay 2D and press "Apply". You may have to change "Coloring" from "Solid Color" to "Temperature".
You can also use:
Select Filters->Alphabetical->Glyph, chose "Glyph Type" "Box", "X Length"=0.25 / "Y Length"=0.25 / "Z Length"=0.01, check "Edit" besides "Set Scale Factor", "Set Scale Factor"=1, and press "Apply".