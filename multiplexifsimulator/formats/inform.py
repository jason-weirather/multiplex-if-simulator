from multiplexifsimulator import FrameEmitter
from skimage.external.tifffile import TiffWriter
import numpy as np
import pandas as pd
import os
#_long_names = {
#    'TUMOR':'CYTOK (Opal 520)',
#    'T-CELL':'CD3 (Opal 540)',
#    'PDL1':'PDL1 (Opal 570)',
#    'PD1':'PD1 (Opal 620)',
#    'MACROPHAGE':'CD68 (Opal 650)',
#    'FOXP3':'FOXP3 (Opal 690)'
#}
tif_tag_integer = {'ImageDescription':'270'}
_nuc_xml = '<?xml version="1.0" encoding="utf-16"?>\r\n<SegmentationImage>\r\n <Version>1</Version>\r\n <CompartmentType>Nucleus</CompartmentType>\r\n</SegmentationImage>'
_mem_xml = '<?xml version="1.0" encoding="utf-16"?>\r\n<SegmentationImage>\r\n <Version>1</Version>\r\n <CompartmentType>Membrane</CompartmentType>\r\n</SegmentationImage>'
_pro_xml = '<?xml version="1.0" encoding="utf-16"?>\r\n<ProcessRegionImage>\r\n <Version>1</Version>\r\n</ProcessRegionImage>'
class FrameEmitterInForm(FrameEmitter):
    """
    Generate mIF data similar in shape to InForm Exports.

    Extends **FrameEmitter**
    """
    def __init__(self,shape=(1040,1392),
                      cell_steps=17,
                      boundary_steps=10,
                      ):
        super(FrameEmitterInForm, self).__init__(shape=shape,cell_steps=cell_steps,boundary_steps=boundary_steps)
        return        
    def save_custom_mask(self,path,x_intercept_fraction=0.667):
      """
      Save a Tumor file 
      """
      _v = np.zeros((self.shape[0],self.shape[1],4)).astype(np.uint8)
      span_range = int(_v.shape[1]*x_intercept_fraction)
      _v[:,:span_range,] = (255,221,221,255)
      with TiffWriter(path) as tif:
        tif.save(_v,compress=9)
    def save_invasive_margin_mask(self,path,x_intercept_fraction=0.667,width_pixels=10):
      """
      Save a Invsive Margin file 
      """
      _v = np.zeros((self.shape[0],self.shape[1],4)).astype(np.uint8)
      span_range = int(_v.shape[1]*x_intercept_fraction)
      left_add = int(width_pixels/2)
      right_add = width_pixels-left_add
      _v[:,span_range-left_add:span_range+right_add,] = (255,255,255,255)
      with TiffWriter(path) as tif:
        tif.save(_v,compress=9)
    def save_binary_seg_maps(self,path,processed_image=True,regions=False,x_intercept_fraction=0.667):
        """
        Save a binary_seg_map image based on the **.make_cell_image_()** images in a python-readable format
        """
        region_image = np.zeros((self.shape[0],self.shape[1])).astype(int)
        span_range = int(region_image.shape[1]*x_intercept_fraction)
        region_image[:,:span_range] = 1
        region_image[:,span_range:] = 2
        _region_xml = '''<?xml version="1.0" encoding="utf-16"?>\r\n<TissueClassMap>\r\n <Version>1</Version>\r\n <Entry>\r\n  <Name>Tumor</Name>\r\n  <Color>255,255,0,0</Color>\r\n  <ID>746136c0-5450-440a-aff7-0969f493f5e3</ID>\r\n </Entry>\r\n <Entry>\r\n  <Name>Stroma</Name>\r\n  <Color>255,0,255,0</Color>\r\n  <ID>6c77d103-7122-4239-b588-5287ddfd028f</ID>\r\n </Entry>\r\n</TissueClassMap>'''

        with TiffWriter(path) as tif:
            if regions:
              tif.save(region_image.astype(np.uint32),
                compress=9, extratags=[(tif_tag_integer['ImageDescription'],'s',0,_region_xml,True)]
                )
            tif.save(self.nucleus_image.astype(np.uint32), 
                     compress=9, extratags=[(tif_tag_integer['ImageDescription'],'s',0,_nuc_xml,True)])
            tif.save(self.edge_image.astype(np.uint32), 
                     compress=9, extratags=[(tif_tag_integer['ImageDescription'],'s',0,_mem_xml,True)])
            if processed_image:
                tif.save(self.processed_image.astype(np.uint32), 
                         compress=9, extratags=[(tif_tag_integer['ImageDescription'],'s',0,_pro_xml,True)])
    def save_component(self,path,frame_name):
        """
        Save a component image
        """
        with TiffWriter(path) as tif:
            for channel in self.components:     
                tif.save(self.components[channel].astype(np.float32), 
                     compress=9, extratags=[(tif_tag_integer['ImageDescription'],
                                             's',
                                             0,
                                             _get_description(channel,frame_name),
                                             True)])

    def make_inform_frame(self,
                          image_folder,
                          frame_name,
                          annotations=None,
                          gimp_tsi=('Tumor','Invasive_Margin'),
                          custom_label=None,
                          binary_names=['PDL1','PD1'],
                          phenotype_strategy=1):
        """
        Save the inform 'cell_seg_data.txt', 'score_data.txt' and 'binary_seg_maps.tif' 
        to  **<image_folder>/**

        annotations can be either 'GIMP TSI', or 'GIMP CUSTOM' or 'InForm'
        gimp_tsi is which images to output.  can be either Tumor or Tumor and Invsive_Margin
        custom_label can be anything
        """
        if self.cell_model is None:
          raise ValueError("Need to set the model cells to make calls")
        cell_seg = _construct_cell_seg(self.cell_model,phenotype_strategy)
        score = _construct_score(self.cell_model,binary_names,annotations)
        score.loc[:,'Sample Name'] = frame_name
        if not os.path.exists(image_folder):
            os.makedirs(image_folder)
        cell_seg.to_csv(os.path.join(image_folder,frame_name+'_cell_seg_data.txt'),
            index=False,sep="\t")
        score.to_csv(os.path.join(image_folder,frame_name+'_score_data.txt'),
            index=False,sep="\t")
        self.make_cell_image()
        bfile = os.path.join(image_folder,frame_name+'_binary_seg_maps.tif')
        self.save_binary_seg_maps(bfile,processed_image=True,regions=True if annotations=='InForm' else False)
        cfile = os.path.join(image_folder,frame_name+'_component_data.tif')
        self.make_component_image(nucleus_width=3,
                                  membrane_width=6,
                                  DAPI=True,
                                  verbose=True,
                                  gaussian_filter_sigma=4,)
        self.save_component(cfile,frame_name)
        if annotations is not None:
          if annotations == 'GIMP TSI':
            if 'Tumor' not in gimp_tsi:
               raise ValueError("need a Tumor mask with this annotation strategy")
            self.save_custom_mask(os.path.join(image_folder,frame_name+'_Tumor.tif'))
            if 'Invasive_Margin' in gimp_tsi:
              self.save_invasive_margin_mask(os.path.join(image_folder,frame_name+'_Invasive_Margin.tif'))
          elif annotations == 'GIMP CUSTOM':
            if custom_label is None:
              raise ValueError("need a custom label if we are setting a custom mask")
            self.save_custom_mask(os.path.join(image_fold,frame_name+'_'+str(custom_label)+'.tif'))
          elif annotations == 'InForm':
            1==1 # we are good.  This case is taken care of in binary seg map creation
          elif annotations != 'None':
            raise ValueError("unknown annotation strategy")
        return #path,cell_seg,score
        

def _construct_cell_seg(cell_model,phenotype_strategy):
    def _fill_chunk(name):
      return [
          'Entire Cell '+name+' Mean (Normalized Counts, Total Weighting)',
          'Nucleus '+name+' Mean (Normalized Counts, Total Weighting)',
          'Membrane '+name+' Mean (Normalized Counts, Total Weighting)',
      ]
    fill = []
    headings_by_marker1 = {}
    for binary_name in cell_model.binary_names_to_channels.keys(): 
      headings_by_marker1[binary_name] = _fill_chunk(cell_model.binary_names_to_channels[binary_name])
      fill+=_fill_chunk(binary_name)
    headings_by_marker2 = {}
    for phenotype_name in cell_model.phenotypes_to_channels.keys():
      headings_by_marker2[phenotype_name] = _fill_chunk(cell_model.phenotypes_to_channels[phenotype_name])
      fill+=_fill_chunk(cell_model.phenotypes_to_channels[phenotype_name])

    header = ['Cell ID','Cell X Position','Cell Y Position','Nucleus Area (pixels)',
          'Nucleus Area (percent)','Nucleus Compactness','Nucleus Minor Axis',
          'Nucleus Major Axis']+fill+[
          'Membrane Area (pixels)','Membrane Area (percent)','Membrane Compactness',
          'Membrane Minor Axis','Membrane Major Axis',
          'Entire Cell Area (pixels)',
          'Entire Cell Area (percent)','Entire Cell Compactness',
          'Entire Cell Minor Axis','Entire Cell Major Axis','Phenotype','Confidence']
    cs = cell_model.cells.copy().rename(columns={'id':'Cell ID',
                                  'x':'Cell X Position',
                                  'y':'Cell Y Position'})

    fill_in = [x for x in header if x not in cs.columns]
    for col in fill_in: cs[col] = 1
    for marker in headings_by_marker1:
      for column_name in headings_by_marker1[marker]:
          cs.loc[cs[marker]=='-',column_name] = 0
    for phenotype_name in cell_model.phenotypes_to_channels.keys():
      for column_name in headings_by_marker2[phenotype_name]:
        if phenotype_name in cs['phenotype_label1']:
          cs.loc[cs['phenotypes_label1']==phenotype_name,headings_by_marker2[phenotype_name]] = 0
        if phenotype_name in cs['phenotype_label2']:
          cs.loc[cs['phenotypes_label2']==phenotype_name,headings_by_marker2[phenotype_name]] = 0
    if phenotype_strategy==1:
      cs = cs.rename(columns={'phenotype_label1':'Phenotype'}).drop(columns=['phenotype_label2'])
    else:
      cs = cs.rename(columns={'phenotype_label2':'Phenotype'}).drop(columns=['phenotype_label1'])

    output = cs.drop(columns=list(cell_model.binary_names_to_channels.keys()))
    return output

def _construct_score(cell_model,binary_names,annotations):
  if len(binary_names) == 0: return None
  if len(binary_names) > 2: 
    raise ValueError("Inform exports dont support more than 2 channels for thresholding to binary phenotypes")
  if len(binary_names) ==2:
    header = ['Path','Sample Name','Tissue Category',
          'First Cell Compartment','First Stain Component',
          'Second Cell Compartment','Second Stain Component',
          'Double Negative',
          'Single '+cell_model.binary_names_to_channels[binary_names[0]],
          'Single '+cell_model.binary_names_to_channels[binary_names[1]],
          'Double Positive',
          'Tissue Category Area (Percent)',
          'Number of Cells',
          cell_model.binary_names_to_channels[binary_names[0]]+' Threshold',
          cell_model.binary_names_to_channels[binary_names[1]]+' Threshold','Lab ID','Slide ID',
          'TMA Sector','TMA Row','TMA Column','TMA Field',
          'inForm 2.1.5430.24864']
    if annotations=='InForm':
      score = pd.DataFrame([len(header)*[0.5],len(header)*[0.5]],columns=header)
    else:
      score = pd.DataFrame([len(header)*[0.5]],columns=header)
    score['Path'] = '/location'
    score['Sample Name'] = 'sample_name'
    score['Tissue Category'] = ''
    score['First Stain Component'] = cell_model.binary_names_to_channels[binary_names[0]]
    score['Second Stain Component'] = cell_model.binary_names_to_channels[binary_names[1]]
    score['First Cell Compartment'] = 'Membrane'
    score['Second Cell Compartment'] = 'Nucleus'
    if annotations=='InForm':
      score.iloc[0,2] = 'Tumor'
      score.iloc[1,2] = 'Stroma'
    return score
  if len(binary_names) == 1: # The single channel case
    header = [
      'Path'
      'Sample Name',
      'Tissue Category',
      'Cell Compartment',
      'Stain Component',
      'Positivity',
      'Tissue Category Area (Percent)',
      'Number of Cells',
      'Positivity Threshold',
      'Lab ID',
      'Slide ID',
      'Annotation ID',
      'TMA Sector',
      'TMA Row',
      'TMA Column',
      'TMA Field',
      'inForm 2.4.6921.16061'
    ]
    if annotations=='InForm':
      score = pd.DataFrame([len(header)*[0.5],len(header)*[0.5]],columns=header)
    else:
      score = pd.DataFrame([len(header)*[0.5]],columns=header)
    score['Path'] = '/location'
    score['Sample Name'] = 'sample_name'
    score['Tissue Category'] = ''
    score['Stain Component'] = cell_model.binary_names_to_channels[binary_names[1]]
    score['Cell Compartment'] = 'Nucleus'
    if annotations=='InForm':
      score.iloc[0,2] = 'Tumor'
      score.iloc[1,2] = 'Stroma'
    return score
def _get_description(channel_name,frame_name):
    istr = '<?xml version="1.0" encoding="utf-16"?>\r\n<PerkinElmer-QPI-ImageDescription>\r\n  <DescriptionVersion>2</DescriptionVersion>\r\n  <AcquisitionSoftware>Mantra</AcquisitionSoftware>\r\n  <ImageType>FullResolution</ImageType>\r\n  <Identifier>d512eea2-f3fd-4ba5-92ed-69b9c9aa9bf2</Identifier>\r\n  <SlideID>'+\
        frame_name+'</SlideID>\r\n  <Barcode />\r\n  <ComputerName>LAPTOP</ComputerName>\r\n  <IsUnmixedComponent>True</IsUnmixedComponent>\r\n  <ExposureTime>16969</ExposureTime>\r\n  <SignalUnits>33</SignalUnits>\r\n  <Name>'+\
        channel_name+'</Name>\r\n  <Color>255,255,0</Color>\r\n  <Objective>20x [20x]</Objective>\r\n  <ValidationCode>80D81BB521952BDCC3ECA60A218A344B5FA618CF97F07F9D7E462CFAC6F50BA0</ValidationCode>\r\n</PerkinElmer-QPI-ImageDescription>'
    return(istr)