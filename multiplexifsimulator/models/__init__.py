import pandas as pd
from random import random
class SlideModelGeneric(object):
    """
    A Generic class to spawn cell frames based on different models

    General properties:
        - **shape** (tuple) the integers (y,x) dimensions
        - **cell_width** (int)
        - **offset** (int) distnace to offset from side and every other row
        - **cells** (int) the pandas.DataFrame with cell and phenotype information
    """
    def __init__(self,shape,cell_width=10):
        self.shape = shape
        self.cell_width = cell_width
        self.offset=5
        self.cells = self._initialize()
        return

    def _initialize(self):
        cells = []
        for y in range(self.offset,self.shape[0],self.cell_width):
            for x in range(self.offset,self.shape[1],self.cell_width):
                cells.append([x,y])
        cells = pd.DataFrame(cells,columns=['x','y']).reset_index().\
            rename(columns={'index':'id'})
        cells['x'] = cells.apply(lambda x: x['x'] if (x['y']-self.offset)%(self.cell_width*2)==0 else x['x']-self.offset,1)
        cells['phenotype_label'] = 'OTHER'
        return cells

    def expanded_cells(self):
        """
        Returns:
            pandas.DataFrame with all phenotypes expanded in the phenotype_label
        """
        ec = self.cells.copy()
        scored = [x for x in self.cells.columns if x not in ['id','x','y','phenotype_label']]
        #print(scored)
        for name in scored:
            ec['phenotype_label'] = ec.apply(lambda x:
                x['phenotype_label']+' '+name+x[name]
            ,1)
        ec = ec.drop(columns=scored)
        return ec


    def fill_uniform(self,
        column_name = 'phenotype_label',
        fill_label='TUMOR',
        fill_probability=0.5,
        condition=lambda x: 1==1):
        """
        Updates the **cells** property with added cell phenotypes uniformally randomly distributed
        """
        self.cells[column_name] = self.cells.apply(lambda x: 
                fill_label if condition(x) and \
                     random() < fill_probability else x[column_name]
            ,1)


    def fill_gradient_margin(self,
                        column_name = 'phenotype_label',
                        fill_label='TUMOR',
                        axis='x',
                        breaks=[1/2,3/4],
                        fill_probability=[1,0],
                        condition=lambda x: 1==1):
        """
        Set a tme with a loose margin between the breaks
            set_fraction is how frequently to apply the tumor label
        """
        iwidth = self.shape[0] if axis == 'y' else self.shape[1]
        ## Do the left
        self.cells[column_name] = self.cells.apply(lambda x: 
                fill_label if x[axis] < iwidth*breaks[0] and \
                     condition(x) and \
                     random() < fill_probability[0] else x[column_name]
            ,1)
        ## Do the right
        self.cells[column_name] = self.cells.apply(lambda x: 
                fill_label if x[axis] >= iwidth*breaks[1] and \
                     condition(x) and \
                     random() < fill_probability[1] else x[column_name]
            ,1)

        ## Now do the gradient
        def _prob(iwidth,value,breaks):
            return ((value/iwidth)-breaks[0])/(breaks[1]-breaks[0])
        p1 = 1-fill_probability[0]
        p2 = 1-fill_probability[1]
        self.cells[column_name] = self.cells.apply(lambda v: 
                fill_label if v[axis] >= iwidth*breaks[0] and \
                     v[axis] < iwidth*breaks[1] and \
                     condition(v) and \
                     p1+_prob(iwidth,v[axis],breaks)*(p2-p1) < random() else v[column_name]
            ,1)


