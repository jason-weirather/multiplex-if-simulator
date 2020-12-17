from pythologist_image_utilities import watershed_image, map_image_ids, image_edges
import numpy as np
import pandas as pd
import sys, math
from scipy.ndimage.filters import gaussian_filter
class FrameEmitter(object):
    """
    Generic class for generating a mIF dataset

    Common properties:
        - **shape** (tuple) - the size of the image integers (y,x)
        - **cell_steps** (int) - how many pixels to try to fill out in the process of defining a cell
        - **boundary_steps** (int) - how many pixels farther than the image to fill out to define the active region

    """
    def __init__(self,*args,**kwargs):
        self.shape = kwargs['shape']
        self.cell_steps = kwargs['cell_steps']
        self.boundary_steps = kwargs['boundary_steps']
        #self.export=kwargs['export']
        #self.locations = None
        #self.phenotypes1= []
        #self.phenotypes2= []
        self.components = {} 
        self.cells = None
        return
    def set_cell_model(self,cell_model):
        """
        Set the cell positions with a cell model.

        """
        self.cell_model = cell_model
        return

    def make_cell_image(self):
        """
        Requires the 'cell_model' be set by **\*.set_cell_model(model)**

        Sets properties:
            - cell_image (numpy.array)
            - nucleus_image (numpy.array)
            - edge_image (numpy.array)
            - processed_image (numpy.array)

        Returns:
            - cell_image
            - edge_image
            - processed_image
        """
        if self.cell_model is None:
            raise ValueError("cell_model is required to be set. see the method .set_cell_model(model)")
        nuc = np.zeros(self.shape).astype(np.uint32)
        # Initialize the image to a pixel at the centroid
        for index,row in self.cell_model.cells.iterrows():
            #print((index,list(row.index)))
            nuc[row['y']][row['x']] = index
        nuc1 = nuc.copy() # work on the cell_image
        nuc2 = nuc.copy() # work on the nucleus_image
        sys.stderr.write("Making cell image\n")
        _total_count = self.cell_model.cells.shape[0]
        # Preserve a point list of the non-centroid points that we can watershed into.
        finish = map_image_ids(nuc1,remove_zero=False).query('id==0').apply(lambda x: (x['x'],x['y']),1)
        if finish.shape[0]==0: finish = []
        for index,row in self.cell_model.cells.iterrows():
            sys.stderr.write("Making cell "+str(index)+" of "+str(_total_count)+"         \r")
            start = [(row['x'],row['y'])] # Shared starting point

            # Fill in the cell_image
            nuc1 = watershed_image(nuc1,list(start),list(finish),fill_value=row.name,steps=self.cell_steps,border=0)

            # Fill in the nucleus
            nuc2 = watershed_image(nuc2,list(start),list(finish),fill_value=row.name,steps=math.ceil(self.cell_steps/4),border=0)

        sys.stderr.write("\n")
        self.nucleus_image = nuc2
        self.cell_image = nuc1
        self.edge_image = image_edges(nuc1)

        # Work on the processed image starting from the cell_image
        start = map_image_ids(nuc1,remove_zero=False).query('id!=0').apply(lambda x: (x['x'],x['y']),1)
        finish = map_image_ids(nuc1,remove_zero=False).query('id==0').apply(lambda x: (x['x'],x['y']),1)
        if start.shape[0]==0: start = []
        if finish.shape[0]==0: finish = []
        temp = watershed_image(nuc1,list(start),list(finish),steps=self.boundary_steps,border=0)
        self.processed_image = temp.astype(np.bool).astype(np.uint8)
        return

    def make_component_image(self,nucleus_width=5,membrane_width=20,DAPI=True,verbose=True,gaussian_filter_sigma=4):
        components = {}
        if DAPI:
            if verbose: sys.stderr.write("DAPI\n")
            components['DAPI'] = _generate_components(self.shape,self.cell_model.cells,
                                                      nucleus_width=nucleus_width,
                                                      membrane_width=membrane_width,gaussian_filter_sigma=gaussian_filter_sigma)

        for phenotype in self.cell_model.phenotypes_to_channels.keys():
            #if  phenotype in ignore_phenotypes: continue
            channel_name = self.cell_model.phenotypes_to_channels[phenotype]
            if verbose: sys.stderr.write(channel_name+"\n")
            if phenotype in self.cell_model.phenotypes1:
                p = self.cell_model.cells[self.cell_model.cells['phenotype_label1']==phenotype]
            elif phenotype in self.cell_model.phenotypes2:
                p = self.cell_model.cells[self.cell_model.cells['phenotype_label2']==phenotype]
            else: 
                raise ValueError("Trying to use a nonexistant phenotype.")
            #p = p.merge(self.cell_model.cells,left_index=True,right_index=True).drop(columns='phenotype')
            components[channel_name] = _generate_components(self.shape,p,
                                                         nucleus_width=nucleus_width,
                                                         membrane_width=membrane_width,gaussian_filter_sigma=gaussian_filter_sigma)
        for binary_name in self.cell_model.binary_names_to_channels.keys():
            channel_name = self.cell_model.binary_names_to_channels[binary_name]
            if verbose: sys.stderr.write(channel_name+"\n")
            p = self.cell_model.cells[self.cell_model.cells[binary_name]=='+']
            components[channel_name] = _generate_components(self.shape,p,
                                                         nucleus_width=nucleus_width,
                                                         membrane_width=membrane_width,gaussian_filter_sigma=gaussian_filter_sigma)
        self.components = components
        return components

def celldataframe_to_cell_model(cdf,frame_name=None):
    raise ValueError("Not currently supported")
    if frame_name is None and len(cdf['frame_name'].unique()) > 1:
        raise ValueError("set no more than one frame at a time")
    cdf = cdf.loc[cdf['frame_name']==frame_name].copy()
    return cdf[['cell_index','x','y','phenotype_label','scored_calls']].\
        set_index(['cell_index','x','y','phenotype_label'])['scored_calls'].\
        apply(pd.Series).applymap(lambda x: '-' if x==0 else '+').reset_index().\
        rename(columns={'cell_index':'id'})

def _generate_components(shape,locations,nucleus_width=5,membrane_width=20,gaussian_filter_sigma=4):
    y, x = np.indices(shape)
    image = np.zeros(shape)
    z = 0
    r1 = nucleus_width
    for i,r in locations.sample(frac=1).iterrows():
        z += 1
        x1 = r['x']
        y1 = r['y']
        mask_circle1 = (x - x1)**2 + (y - y1)**2 < r1**2
        image = np.add(image, mask_circle1)
        image = np.add(image, mask_circle1)
        image = np.add(image, mask_circle1)
    r1 = membrane_width
    for i,r in locations.sample(frac=1).iterrows():
        z += 1
        x1 = r['x']
        y1 = r['y']
        mask_circle1 = (x - x1)**2 + (y - y1)**2 < r1**2
        image = np.add(image, mask_circle1)
    image = np.add(image,gaussian_filter(image,membrane_width*2))
    return gaussian_filter(image,gaussian_filter_sigma)
