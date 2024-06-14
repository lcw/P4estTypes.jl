from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()
forest = XMLPartitionedUnstructuredGridReader(registrationName='forest_2d.pvtu', FileName=['forest_2d.pvtu'])
forest.CellArrayStatus = ['treeid', 'level', 'mpirank']

# get active view
view = GetActiveViewOrCreate('RenderView')

# show data in view
display = Show(forest, view, 'UnstructuredGridRepresentation')

# change representation type
display.SetRepresentationType('Surface With Edges')

# set scalar coloring
ColorBy(display, ('CELLS', 'mpirank'))

display.SetScalarBarVisibility(view, True)

# get color transfer function/color map for 'mpirank'
lut = GetColorTransferFunction('mpirank')

# Properties modified on lut
lut.InterpretValuesAsCategories = 1
lut.AnnotationsInitialized = 1
lut.Annotations = ['0', '0', '1', '1', '2', '2', '3', '3', '4', '4']
lut.IndexedColors = [1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0]
lut.IndexedOpacities = [1.0, 1.0, 1.0, 1.0, 1.0]
lut.ApplyPreset('Brewer Diverging Purple-Orange (5)', True)
colorbar = GetScalarBar(lut, view)
colorbar.Title = 'MPI Rank'

# update the view to ensure updated data information
view.Update()

view.OrientationAxesVisibility = 0
view.InteractionMode = '2D'
view.CameraPosition = [1.5, 1.0, 1.0]
view.CameraFocalPoint = [1.5, 1.0, 0.0]
view.CameraParallelScale = 1
view.ViewSize = [900, 600]

SaveScreenshot('forest_2d.png', magnification=1, quality=100, view=view)
