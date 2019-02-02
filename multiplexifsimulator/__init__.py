from pythologist.formats.utilities import watershed_image, map_image_ids, image_edges
import numpy as np
import sys
class FrameEmitter(object):
    def __init__(self,*args,**kwargs):
        self.shape = kwargs['shape']
        self.cell_steps = kwargs['cell_steps']
        self.boundary_steps = kwargs['boundary_steps']
        return
    def set_cell_coordinates(self,cell_df):
        # cell_df: <x> <y> <id (optional)>
        if cell_df['x'].max() >= self.shape[1]:
            sys.stderr.write("Warning: shape is larger than coordinates... trimming")
            cell_df = cell_df.loc[cell_df['x']<self.shape[1]]
        if cell_df['y'].max() >= self.shape[0]:
            sys.stderr.write("Warning: shape is larger than coordinates... trimming")
            cell_df = cell_df.loc[cell_df['y']<self.shape[0]]
        if 'x' not in cell_df.columns or 'y' not in cell_df.columns: 
            raise ValueError('x and y must be in dataframe')
        if 'id' in cell_df.columns:
            self.locations = cell_df.groupby('id').first()[['x','y']].reset_index()
            return
        self.locations = cell_df[['x','y']].drop_duplicates().sample(frac=1).reset_index(drop=True)
        self.locations.index = [x+1 for x in self.locations.index]
        self.locations.index.name = 'id'
    def make_cell_image(self):
        # requires 'locations' property set
        nuc = np.zeros(self.shape).astype(np.uint32)
        for index,row in self.locations.iterrows():
            nuc[row['y']][row['x']] = row.name
        start = map_image_ids(nuc,remove_zero=False).query('id!=0').apply(lambda x: (x['x'],x['y']),1)
        finish = map_image_ids(nuc,remove_zero=False).query('id==0').apply(lambda x: (x['x'],x['y']),1)
        nuc1 = watershed_image(nuc,list(start),list(finish),steps=self.cell_steps)
        self.cell_image = nuc1
        self.edge_image = image_edges(nuc1,seek_distance=1)
        start = map_image_ids(nuc1,remove_zero=False).query('id!=0').apply(lambda x: (x['x'],x['y']),1)
        finish = map_image_ids(nuc1,remove_zero=False).query('id==0').apply(lambda x: (x['x'],x['y']),1)
        temp = watershed_image(nuc1,list(start),list(finish),steps=self.boundary_steps)
        self.processed_image = temp.astype(np.bool).astype(np.uint8)
        return self.cell_image, self.edge_image, self.processed_image