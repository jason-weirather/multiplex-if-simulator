# multiplex-if-simulator

Simulate multiplex IF data 

##### Read the docs: https://jason-weirather.github.io/multiplex-if-simulator/

# Quickstart

### Generate a small Infiltrated mIF image

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

![Image of phenotypes](https://github.com/jason-weirather/multiplex-if-simulator/raw/master/images/phenotypes.png?raw=true)

![Image of cell map](https://github.com/jason-weirather/multiplex-if-simulator/raw/master/images/cell_map.png?raw=true)

![Image of edge map](https://github.com/jason-weirather/multiplex-if-simulator/raw/master/images/edge_map.png?raw=true)

![Image of processed area](https://github.com/jason-weirather/multiplex-if-simulator/raw/master/images/processed_image.png?raw=true)


### Generate an dataset shaped similar to those of InForm outputs

This will have 4 samples with 3 images each from 3 different model types
It will create both a python compatible `Example` folder and an R compatible
`Example-R` folder

```python
from multiplexifsimulator.formats.inform import FrameEmitterInForm
from multiplexifsimulator.models.general import SlideModelExcluded as SME
from multiplexifsimulator.models.general import SlideModelInfiltrated as SMI
from multiplexifsimulator.models.general import SlideModelUniform as SMU

shape = (400,600)
fe = FrameEmitterInForm(shape=shape)
for n in range(1,5):
    for i in range(1,4):
        c = SME(shape=shape).cells
        fe.make_inform_frame(c,'Example','E'+str(n),str(i))
        fe.make_inform_frame(c,'Example-R','E'+str(n),str(i),r_format=True)
        c = SMI(shape=shape).cells
        fe.make_inform_frame(c,'Example','I'+str(n),str(i))
        fe.make_inform_frame(c,'Example-R','I'+str(n),str(i),r_format=True)
        c = SMU(shape=shape).cells
        fe.make_inform_frame(c,'Example','U'+str(n),str(i))
        fe.make_inform_frame(c,'Example-R','U'+str(n),str(i),r_format=True)
```

```
$ tree Example
Example
├── E1
│   ├── E1_1_binary_seg_maps.tif
│   ├── E1_1_cell_seg_data.txt
│   ├── E1_1_score_data.txt
│   ├── E1_2_binary_seg_maps.tif
│   ├── E1_2_cell_seg_data.txt
│   ├── E1_2_score_data.txt
│   ├── E1_3_binary_seg_maps.tif
│   ├── E1_3_cell_seg_data.txt
│   └── E1_3_score_data.txt
├── E2
│   ├── E2_1_binary_seg_maps.tif
│   ├── E2_1_cell_seg_data.txt
│   ├── E2_1_score_data.txt
│   ├── E2_2_binary_seg_maps.tif
│   ├── E2_2_cell_seg_data.txt
│   ├── E2_2_score_data.txt
│   ├── E2_3_binary_seg_maps.tif
│   ├── E2_3_cell_seg_data.txt
│   └── E2_3_score_data.txt
├── E3
│   ├── E3_1_binary_seg_maps.tif
│   ├── E3_1_cell_seg_data.txt
│   ├── E3_1_score_data.txt
│   ├── E3_2_binary_seg_maps.tif
│   ├── E3_2_cell_seg_data.txt
│   ├── E3_2_score_data.txt
│   ├── E3_3_binary_seg_maps.tif
│   ├── E3_3_cell_seg_data.txt
│   └── E3_3_score_data.txt
├── E4
│   ├── E4_1_binary_seg_maps.tif
│   ├── E4_1_cell_seg_data.txt
│   ├── E4_1_score_data.txt
│   ├── E4_2_binary_seg_maps.tif
│   ├── E4_2_cell_seg_data.txt
│   ├── E4_2_score_data.txt
│   ├── E4_3_binary_seg_maps.tif
│   ├── E4_3_cell_seg_data.txt
│   └── E4_3_score_data.txt
├── I1
│   ├── I1_1_binary_seg_maps.tif
│   ├── I1_1_cell_seg_data.txt
│   ├── I1_1_score_data.txt
│   ├── I1_2_binary_seg_maps.tif
│   ├── I1_2_cell_seg_data.txt
│   ├── I1_2_score_data.txt
│   ├── I1_3_binary_seg_maps.tif
│   ├── I1_3_cell_seg_data.txt
│   └── I1_3_score_data.txt
├── I2
│   ├── I2_1_binary_seg_maps.tif
│   ├── I2_1_cell_seg_data.txt
│   ├── I2_1_score_data.txt
│   ├── I2_2_binary_seg_maps.tif
│   ├── I2_2_cell_seg_data.txt
│   ├── I2_2_score_data.txt
│   ├── I2_3_binary_seg_maps.tif
│   ├── I2_3_cell_seg_data.txt
│   └── I2_3_score_data.txt
├── I3
│   ├── I3_1_binary_seg_maps.tif
│   ├── I3_1_cell_seg_data.txt
│   ├── I3_1_score_data.txt
│   ├── I3_2_binary_seg_maps.tif
│   ├── I3_2_cell_seg_data.txt
│   ├── I3_2_score_data.txt
│   ├── I3_3_binary_seg_maps.tif
│   ├── I3_3_cell_seg_data.txt
│   └── I3_3_score_data.txt
├── I4
│   ├── I4_1_binary_seg_maps.tif
│   ├── I4_1_cell_seg_data.txt
│   ├── I4_1_score_data.txt
│   ├── I4_2_binary_seg_maps.tif
│   ├── I4_2_cell_seg_data.txt
│   ├── I4_2_score_data.txt
│   ├── I4_3_binary_seg_maps.tif
│   ├── I4_3_cell_seg_data.txt
│   └── I4_3_score_data.txt
├── U1
│   ├── U1_1_binary_seg_maps.tif
│   ├── U1_1_cell_seg_data.txt
│   ├── U1_1_score_data.txt
│   ├── U1_2_binary_seg_maps.tif
│   ├── U1_2_cell_seg_data.txt
│   ├── U1_2_score_data.txt
│   ├── U1_3_binary_seg_maps.tif
│   ├── U1_3_cell_seg_data.txt
│   └── U1_3_score_data.txt
├── U2
│   ├── U2_1_binary_seg_maps.tif
│   ├── U2_1_cell_seg_data.txt
│   ├── U2_1_score_data.txt
│   ├── U2_2_binary_seg_maps.tif
│   ├── U2_2_cell_seg_data.txt
│   ├── U2_2_score_data.txt
│   ├── U2_3_binary_seg_maps.tif
│   ├── U2_3_cell_seg_data.txt
│   └── U2_3_score_data.txt
├── U3
│   ├── U3_1_binary_seg_maps.tif
│   ├── U3_1_cell_seg_data.txt
│   ├── U3_1_score_data.txt
│   ├── U3_2_binary_seg_maps.tif
│   ├── U3_2_cell_seg_data.txt
│   ├── U3_2_score_data.txt
│   ├── U3_3_binary_seg_maps.tif
│   ├── U3_3_cell_seg_data.txt
│   └── U3_3_score_data.txt
└── U4
    ├── U4_1_binary_seg_maps.tif
    ├── U4_1_cell_seg_data.txt
    ├── U4_1_score_data.txt
    ├── U4_2_binary_seg_maps.tif
    ├── U4_2_cell_seg_data.txt
    ├── U4_2_score_data.txt
    ├── U4_3_binary_seg_maps.tif
    ├── U4_3_cell_seg_data.txt
    └── U4_3_score_data.txt
```


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