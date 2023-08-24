from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()
forest = XMLPartitionedUnstructuredGridReader(registrationName='forest_3d.pvtu', FileName=['forest_3d.pvtu'])
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
view.CameraPosition = [-21.21422329177474, -15.375102411690401, -8.921290188684675]
view.CameraFocalPoint = [1.5000000000000002, 0.9999999999999998, 2.9999999999999967]
view.CameraViewUp = [-0.2631790038546555, -0.205978948089547, 0.9425017161119629]
view.CameraParallelScale = 7.927230181930915
view.ResetCamera(False)

view.ViewSize = [600, 800]

SaveScreenshot('forest_3d.png', magnification=1, quality=100, view=view)
