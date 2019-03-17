from pythologist_image_utilities import watershed_image, map_image_ids, image_edges
import numpy as np
import pandas as pd
import sys, random
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
        self.locations = None
        self.phenotypes = None
        return
    def set_cell_coordinates(self,cell_df):
        """
        Set the cell positions with a DataFrame.  *id* field is optional.  
        If not set the is will be shuffled.

        ===  ===  ===
        id    x    y
        ===  ===  ===
        1    5    5
        2    25   5
        3    45   5
        4    65   5
        ===  ===  ===

        """
        # set phenotypes if we can
        if 'phenotype_label' in cell_df:
            mc = cell_df.copy()
            phenotypes = mc['phenotype_label'].dropna().unique()
            df = pd.DataFrame(index=mc['id'].index,columns=phenotypes).fillna('-')
            for i,r in mc[['id','phenotype_label']].dropna().iterrows():
                df.loc[i,r['phenotype_label']] = '+'
            self.phenotypes = df.merge(mc.drop(columns=['x','y','phenotype_label']),left_index=True,right_on='id').set_index('id')


        if cell_df['x'].max() >= self.shape[1]:
            sys.stderr.write("Warning: shape is larger than coordinates... trimming")
            cell_df = cell_df.loc[cell_df['x']<self.shape[1]]
        if cell_df['y'].max() >= self.shape[0]:
            sys.stderr.write("Warning: shape is larger than coordinates... trimming")
            cell_df = cell_df.loc[cell_df['y']<self.shape[0]]
        if 'x' not in cell_df.columns or 'y' not in cell_df.columns: 
            raise ValueError('x and y must be in dataframe')
        if 'id' in cell_df.columns:
            self.locations = cell_df.groupby('id').first()[['x','y']]#.reset_index()
            return
        self.locations = cell_df[['x','y']].drop_duplicates().sample(frac=1).reset_index(drop=True)
        self.locations.index = [x+1 for x in self.locations.index]
        self.locations.index.name = 'id'
        # If there is phenotype data set that
    def make_cell_image(self):
        """
        Requires the 'locations' be set by **\*.set_cell_coordinates(cell_df)**

        Sets properties:
            - cell_image (numpy.array)
            - edge_image (numpy.array)
            - processed_image (numpy.array)

        Returns:
            - cell_image
            - edge_image
            - processed_image
        """
        # requires 'locations' property set
        nuc = np.zeros(self.shape).astype(np.uint32)
        for index,row in self.locations.iterrows():
            nuc[row['y']][row['x']] = row.name
        start = map_image_ids(nuc,remove_zero=False).query('id!=0').apply(lambda x: (x['x'],x['y']),1)
        finish = map_image_ids(nuc,remove_zero=False).query('id==0').apply(lambda x: (x['x'],x['y']),1)
        nuc1 = watershed_image(nuc,list(start),list(finish),steps=self.cell_steps,border=0)
        self.cell_image = nuc1
        self.edge_image = image_edges(nuc1)
        start = map_image_ids(nuc1,remove_zero=False).query('id!=0').apply(lambda x: (x['x'],x['y']),1)
        finish = map_image_ids(nuc1,remove_zero=False).query('id==0').apply(lambda x: (x['x'],x['y']),1)
        temp = watershed_image(nuc1,list(start),list(finish),steps=self.boundary_steps,border=0)
        self.processed_image = temp.astype(np.bool).astype(np.uint8)
        return self.cell_image, self.edge_image, self.processed_image
    def make_component_image(self,nucleus_width=5,membrane_width=20,DAPI=True,verbose=True,ignore_phenotypes=['OTHER'],gaussian_filter_sigma=4):
        components = {}
        if DAPI:
            if verbose: sys.stderr.write("DAPI\n")
            components['DAPI'] = _generate_components(self.shape,self.locations,
                                                      nucleus_width=nucleus_width,
                                                      membrane_width=membrane_width,gaussian_filter_sigma=gaussian_filter_sigma)
        for phenotype in self.phenotypes.columns:
            if  phenotype in ignore_phenotypes: continue
            if verbose: sys.stderr.write(phenotype+"\n")
            p = self.phenotypes[self.phenotypes[phenotype]=='+'].\
                rename(columns={phenotype:'phenotype'})[['phenotype']]
            p = p.merge(self.locations,left_index=True,right_index=True).drop(columns='phenotype')
            components[phenotype] = _generate_components(self.shape,p,
                                                         nucleus_width=nucleus_width,
                                                         membrane_width=membrane_width,gaussian_filter_sigma=gaussian_filter_sigma)
        self.components = components
        return components

def celldataframe_to_cell_model(cdf,frame_name=None):
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
