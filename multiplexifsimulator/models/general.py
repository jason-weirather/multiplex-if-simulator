from multiplexifsimulator.models import SlideModelGeneric
class SlideModelExcluded(SlideModelGeneric):
    """
    A slide where there are more T cells outside of the tumor than inside of the
    tumor and an immunosuppressive phenotypes inside the tumor

    .. figure::  ../images/exclusion.png

    """
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self.phenotypes = ['OTHER','T-CELL','TUMOR']
        self.cells['PDL1'] = '-'
        self.cells['PD1'] = '-'
        self.fill_gradient_margin(
                column_name = 'phenotype_label',
                fill_label = 'TUMOR',
                axis='x',
                breaks=[6/13,7/13],
                fill_probability=[1,0])
        self.fill_uniform(
                column_name = 'phenotype_label',
                fill_label = 'OTHER',
                fill_probability=0.05)
        self.fill_gradient_margin(
                column_name = 'phenotype_label',
                fill_label = 'T-CELL',
                axis='x',
                breaks=[2/5,3/5],
                fill_probability=[0.02,0.11])
        self.fill_gradient_margin(
                column_name = 'PDL1',
                fill_label = '+',
                axis='x',
                breaks=[0,1],
                fill_probability=[1,0],
                condition=lambda x: x['phenotype_label']=='OTHER')
        self.fill_gradient_margin(
                column_name = 'PD1',
                fill_label = '+',
                axis='x',
                breaks=[2/5,3/5],
                fill_probability=[0.9,0.05],
                condition=lambda x: x['phenotype_label']=='T-CELL')
        self.fill_uniform(
                column_name = 'PDL1',
                fill_label = '+',
                fill_probability=0.7,
                condition=lambda x: x['phenotype_label']=='TUMOR')

class SlideModelInfiltrated(SlideModelGeneric):
    """
    A slide where there are more T cells outside of the tumor than inside of the
    tumor and an immunosuppressive phenotypes inside the tumor

    .. figure::  ../images/infiltration.png

    """
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self.phenotypes = ['OTHER','T-CELL','TUMOR']
        self.cells['PDL1'] = '-'
        self.cells['PD1'] = '-'
        self.fill_gradient_margin(
                column_name = 'phenotype_label',
                fill_label = 'TUMOR',
                axis='x',
                breaks=[6/13,8/13],
                fill_probability=[1,0])
        self.fill_uniform(
                column_name = 'phenotype_label',
                fill_label = 'OTHER',
                fill_probability=0.05)
        self.fill_gradient_margin(
                column_name = 'phenotype_label',
                fill_label = 'T-CELL',
                axis='x',
                breaks=[2/5,3/5],
                fill_probability=[0.11,0.02])
        self.fill_gradient_margin(
                column_name = 'PDL1',
                fill_label = '+',
                axis='x',
                breaks=[0,1],
                fill_probability=[1,0.02],
                condition=lambda x: x['phenotype_label']=='OTHER')
        self.fill_gradient_margin(
                column_name = 'PD1',
                fill_label = '+',
                axis='x',
                breaks=[2/5,3/5],
                fill_probability=[0.06,0.9],
                condition=lambda x: x['phenotype_label']=='T-CELL')
        self.fill_uniform(
                column_name = 'PDL1',
                fill_label = '+',
                fill_probability=0.71,
                condition=lambda x: x['phenotype_label']=='TUMOR')

class SlideModelUniform(SlideModelGeneric):
    """
    A slide where all cells are uniformally distributed

    .. figure::  ../images/uniform.png

    """
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self.phenotypes = ['OTHER','T-CELL','TUMOR']
        self.cells['PDL1'] ='-'
        self.cells['PD1'] = '-'
        self.fill_uniform(
                column_name = 'phenotype_label',
                fill_label = 'TUMOR',
                fill_probability=0.49)
        self.fill_uniform(
                column_name = 'phenotype_label',
                fill_label = 'T-CELL',
                fill_probability=0.065)
        self.fill_uniform(
                column_name = 'PD1',
                fill_label = '+',
                fill_probability=0.20,
                condition=lambda x: x['phenotype_label']=='T-CELL')
        self.fill_uniform(
                column_name = 'PDL1',
                fill_label = '+',
                fill_probability=0.7,
                condition=lambda x: x['phenotype_label']=='TUMOR')
        self.fill_uniform(
                column_name = 'PDL1',
                fill_label = '+',
                fill_probability=0.28,
                condition=lambda x: x['phenotype_label']=='OTHER')
