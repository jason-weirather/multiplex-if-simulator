# multiplex-if-simulator

Simulate multiplex IF data 

##### Read the docs: https://jason-weirather.github.io/multiplex-if-simulator/

##### Quickstart

Generate a small Infiltrated mIF image

```python
from multiplexifsimulator.models.general import SlideModelInfiltrated as SMI
from multiplexifsimulator.formats.inform import FrameEmitterInForm
from plotnine import *
import matplotlib.pyplot as plt

## Set your image size
shape = (100,150)

## Set the cell coordinates and phenotypes with a model
model = SMI(shape)
## Initialize the emitter
fe = FrameEmitterInForm(shape)
fe.set_cell_coordinates(model.cells)
## Generate the base images in numpy array format
cell_image, edge_image, processed_image = fe.make_cell_image()

## Show the images
(ggplot(model.cells,aes(x='x',y='y',fill='phenotype_label'))
 + geom_point(shape='h',size=5,stroke=0.3)
 + theme_minimal()
 + theme(figure_size=(4.5,4),aspect_ratio=shape[0]/shape[1])
 + xlim(0,shape[1])
 + ylim(0,shape[0])
).draw()
plt.show()
plt.imshow(cell_image,origin='lower')
plt.show()
plt.imshow(edge_image,origin='lower')
plt.show()
plt.imshow(processed_image,origin='lower')
plt.show()
```

![Image of phenotypes](images/phenotypes.png?raw=true)

![Image of cell map](images/cell_map.png?raw=true)

![Image of edge map](images/edge_map.png?raw=true)

![Image of processed area](images/processed_image.png?raw=true)


```{python}

### Generate a binary_seg_map.tif based on a cell seg data file

```python
import matplotlib.pyplot as plt
from multiplexifsimulator.formats.inform import FrameEmitterInForm

# Read in a cell seg file
fname = 'input/MEL3_120116_2_cell_seg_data.txt'
seg = pd.read_csv(fname,sep="\t")

# Create an image emiter that will expand cells a maximum of 17 pixels
#     in each direction and will expand the processed border a further
#     20 pixels
fe = FrameEmitterInForm(shape=(1040, 1392),cell_steps=17,boundary_steps=20)
# Set the cell coordinates
fe.set_cell_coordinates(seg[['Cell X Position',
                             'Cell Y Position',
                             'Cell ID']].rename(
                        columns={'Cell X Position':'x',
                                 'Cell Y Position':'y',
                                 'Cell ID':'id'}))
# Produce numpy arrays
nuc, mem, proc = fe.make_cell_image()

# Display images
plt.imshow(nuc)
plt.show()
plt.imshow(mem)
plt.show()
plt.imshow(proc)
plt.show()

# Save the combined tif
fe.save_binary_seg_maps('output/MEL3_120116_2_binary_seg_maps.tif')
```