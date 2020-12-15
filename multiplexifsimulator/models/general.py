from multiplexifsimulator.models import SlideModelGeneric

phenotypes_to_channels = {
    'TUMOR':'CYTOK (Opal 520)',
    'T-CELL':'CD3 (Opal 540)',
    'MACROPHAGE':'CD68 (Opal 650)'
}
binary_names_to_channels = {
    'PDL1':'PDL1 (Opal 570)',
    'PD1':'PD1 (Opal 620)',
    'FOXP3':'FOXP3 (Opal 690)'
}

class SlideModelExcluded(SlideModelGeneric):
    """
    A slide where there are more T cells outside of the tumor than inside of the
    tumor and an immunosuppressive phenotypes inside the tumor

    .. figure::  ./images/exclusion.png

    """
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self.phenotypes_to_channels = phenotypes_to_channels.copy()
        self.binary_names_to_channels = binary_names_to_channels.copy()
        self.phenotypes1 = ['OTHER','T-CELL','TUMOR']
        self.phenotypes2 = ['OTHER','MACROPHAGE']
        self.cells['PDL1'] = '-'
        self.cells['PD1'] = '-'
        self.cells['FOXP3'] = '-'
        self.fill_gradient_margin(
                column_name = 'phenotype_label1',
                fill_label = 'TUMOR',
                axis='x',
                breaks=[6/13,7/13],
                fill_probability=[1,0])
        self.fill_uniform(
                column_name = 'phenotype_label1',
                fill_label = 'OTHER',
                fill_probability=0.05)
        self.fill_uniform(
                column_name = 'phenotype_label2',
                fill_label = 'MACROPHAGE',
                fill_probability=0.2)
        self.fill_gradient_margin(
                column_name = 'phenotype_label1',
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
                condition=lambda x: x['phenotype_label1']=='OTHER' and x['phenotype_label2']!='MACROPHAGE')
        self.fill_gradient_margin(
                column_name = 'PDL1',
                fill_label = '+',
                axis='x',
                breaks=[0,1],
                fill_probability=[1,0],
                condition=lambda x: x['phenotype_label2']=='MACROPHAGE')
        self.fill_gradient_margin(
                column_name = 'PD1',
                fill_label = '+',
                axis='x',
                breaks=[2/5,3/5],
                fill_probability=[0.9,0.05],
                condition=lambda x: x['phenotype_label1']=='T-CELL')
        self.fill_gradient_margin(
                column_name = 'FOXP3',
                fill_label = '+',
                axis='x',
                breaks=[2/5,3/5],
                fill_probability=[0.9,0.05],
                condition=lambda x: x['phenotype_label1']=='T-CELL')
        self.fill_uniform(
                column_name = 'PDL1',
                fill_label = '+',
                fill_probability=0.7,
                condition=lambda x: x['phenotype_label1']=='TUMOR')

class SlideModelInfiltrated(SlideModelGeneric):
    """
    A slide where there are more T cells outside of the tumor than inside of the
    tumor and an immunosuppressive phenotypes inside the tumor

    .. figure::  ./images/infiltration.png

    """



    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self.phenotypes_to_channels = phenotypes_to_channels.copy()
        self.binary_names_to_channels = binary_names_to_channels.copy()
        self.phenotypes1= ['OTHER','T-CELL','TUMOR']
        self.phenotypes2 = ['OTHER','MACROPHAGE']
        self.cells['PDL1'] = '-'
        self.cells['PD1'] = '-'
        self.cells['FOXP3'] = '-'
        self.fill_gradient_margin(
                column_name = 'phenotype_label1',
                fill_label = 'TUMOR',
                axis='x',
                breaks=[6/13,8/13],
                fill_probability=[1,0])
        self.fill_uniform(
                column_name = 'phenotype_label1',
                fill_label = 'OTHER',
                fill_probability=0.05)
        self.fill_uniform(
                column_name = 'phenotype_label2',
                fill_label = 'MACROPHAGE',
                fill_probability=0.2)
        self.fill_gradient_margin(
                column_name = 'phenotype_label1',
                fill_label = 'T-CELL',
                axis='x',
                breaks=[2/5,3/5],
                fill_probability=[0.2,0.02])
        self.fill_gradient_margin(
                column_name = 'PDL1',
                fill_label = '+',
                axis='x',
                breaks=[0,1],
                fill_probability=[1,0.02],
                condition=lambda x: x['phenotype_label1']=='OTHER' and x['phenotype_label2']!='MACROPHAGE')
        self.fill_gradient_margin(
                column_name = 'PDL1',
                fill_label = '+',
                axis='x',
                breaks=[0,1],
                fill_probability=[1,0.02],
                condition=lambda x: x['phenotype_label2']=='MACROPHAGE')
        self.fill_gradient_margin(
                column_name = 'PD1',
                fill_label = '+',
                axis='x',
                breaks=[2/5,3/5],
                fill_probability=[0.06,0.9],
                condition=lambda x: x['phenotype_label1']=='T-CELL')
        self.fill_gradient_margin(
                column_name = 'FOXP3',
                fill_label = '+',
                axis='x',
                breaks=[2/5,3/5],
                fill_probability=[0.06,0.9],
                condition=lambda x: x['phenotype_label1']=='T-CELL')
        self.fill_uniform(
                column_name = 'PDL1',
                fill_label = '+',
                fill_probability=0.71,
                condition=lambda x: x['phenotype_label1']=='TUMOR')

class SlideModelUniform(SlideModelGeneric):
    """
    A slide where all cells are uniformally distributed

    .. figure::  ./images/uniform.png

    """
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self.phenotypes_to_channels = phenotypes_to_channels.copy()
        self.binary_names_to_channels = binary_names_to_channels.copy()
        self.phenotypes1 = ['OTHER','T-CELL','TUMOR']
        self.phenotypes2 = ['OTHER','MACROPHAGE']
        self.cells['PDL1'] ='-'
        self.cells['PD1'] = '-'
        self.cells['FOXP3'] = '-'
        self.fill_uniform(
                column_name = 'phenotype_label1',
                fill_label = 'TUMOR',
                fill_probability=0.49)
        self.fill_uniform(
                column_name = 'phenotype_label2',
                fill_label = 'MACROPHAGE',
                fill_probability=0.2)
        self.fill_uniform(
                column_name = 'phenotype_label1',
                fill_label = 'T-CELL',
                fill_probability=0.1)
        self.fill_uniform(
                column_name = 'PD1',
                fill_label = '+',
                fill_probability=0.20,
                condition=lambda x: x['phenotype_label1']=='T-CELL')
        self.fill_uniform(
                column_name = 'FOXP3',
                fill_label = '+',
                fill_probability=0.2,
                condition=lambda x: x['phenotype_label1']=='T-CELL')
        self.fill_uniform(
                column_name = 'PDL1',
                fill_label = '+',
                fill_probability=0.7,
                condition=lambda x: x['phenotype_label1']=='TUMOR')
        self.fill_uniform(
                column_name = 'PDL1',
                fill_label = '+',
                fill_probability=0.28,
                condition=lambda x: x['phenotype_label1']=='OTHER' and x['phenotype_label2']!='MACROPHAGE')
        self.fill_uniform(
                column_name = 'PDL1',
                fill_label = '+',
                fill_probability=0.6,
                condition=lambda x: x['phenotype_label2']=='MACROPHAGE')
