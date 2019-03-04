from pythologist_image_utilities import watershed_image, map_image_ids, image_edges
import numpy as np
import pandas as pd
import sys, random
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
        self.edge_image = image_edges(nuc1,seek_distance=1)
        start = map_image_ids(nuc1,remove_zero=False).query('id!=0').apply(lambda x: (x['x'],x['y']),1)
        finish = map_image_ids(nuc1,remove_zero=False).query('id==0').apply(lambda x: (x['x'],x['y']),1)
        temp = watershed_image(nuc1,list(start),list(finish),steps=self.boundary_steps,border=0)
        self.processed_image = temp.astype(np.bool).astype(np.uint8)
        return self.cell_image, self.edge_image, self.processed_image
    def gradient_model1(self):
        xdim = self.shape[1]
        ydim = self.shape[0]
        cells = []
        def _make_gradient1(x,xdim,label):
            # we want 1/3 to 2/3 drop off of tumor
            fill = 'TUMOR'
            if x < xdim/3: return fill
            if x > 2*xdim/3: return label
            if (x-(xdim/3)) > random.random()*xdim/3: return label
            return fill
        def _make_gradient2(x,xdim,label):
            # we want a 1 to 0 gradient of t-cell
            fill = 'T-CELL'
            if (1-x/xdim) < random.random() and random.random() < 0.2: return fill
            return label
        def _make_gradient3(x,y,xdim,ydim,label):
            # we want a 1 to 0 gradient of t-cell
            fill = '+'
            if x/xdim < random.random() and label == 'T-CELL' and y/ydim > 0.5: return '+'
            return '-'
        def _make_gradient4(x,y,xdim,ydim,label):
            # we want a 1 to 0 gradient of t-cell
            fill = '+'
            if x/xdim < random.random() and label == 'OTHER' and y/ydim > 0.3: return '+'
            return '-'

        for y in range(5,ydim,10):
            for x in range(5,xdim,10):
                cells.append([x,y])
        cells = pd.DataFrame(cells,columns=['x','y']).reset_index().\
            rename(columns={'index':'id'})
        ## Lets start everything as OTHER
        cells['phenotype_label'] = 'OTHER'
        ## Lets start PD1 as negative on everything
        cells['PD1'] = '-'
        ## Lets start PDL1 as negative on everything
        cells['PDL1'] = '-'
        cells['phenotype_label'] = cells.apply(lambda x:
            _make_gradient1(x['x'],xdim,x['phenotype_label'])
        ,1)
        cells['phenotype_label'] = cells.apply(lambda x:
            _make_gradient2(x['x'],xdim,x['phenotype_label'])
        ,1)
        # Make a hard PDL1 gradient on the top
        cells['PDL1'] = cells.apply(lambda x: '+' if x['y'] < 2*ydim/3 and x['phenotype_label']=='TUMOR'  and random.random() < 0.9 else '-',1)
        cells['PD1'] = cells.apply(lambda x: 
            _make_gradient3(x['x'],x['y'],xdim,ydim,x['phenotype_label'])
        ,1)
        cells['PDL1'] = cells.apply(lambda x: 
            _make_gradient4(x['x'],x['y'],xdim,ydim,x['phenotype_label'])
        ,1)
        cells = cells.sort_values('phenotype_label').set_index('id').reset_index(drop=True).reset_index().\
            rename(columns={'index':'id'})
        return cells
