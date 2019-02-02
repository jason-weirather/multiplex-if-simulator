# multiplex-if-simulator

Simulate multiplex IF data 

### Generate a binary_seg_map.tif based on a segmentation file

```{python}
import matplotlib.pyplot as plt
from multiplexifsimulator.formats.inform import FrameEmitterInForm

# Read in a cell seg file
fname = 'input/MEL3_120116_2_cell_seg_data.txt'
seg = pd.read_csv(fname,sep="\t")

# Create an image emiter of the proper size
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