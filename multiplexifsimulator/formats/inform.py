from multiplexifsimulator import FrameEmitter
from tifffile import TiffWriter
import numpy as np
import pandas as pd
import os
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
                      boundary_steps=10):
        super(FrameEmitterInForm, self).__init__(shape=shape,cell_steps=cell_steps,boundary_steps=boundary_steps)
        return        
    def save_binary_seg_maps(self,path,processed_image=True):
        """
        Save a binary_seg_map image based on the **.make_cell_image_()** images in a python-readable format
        """
        with TiffWriter(path) as tif:
            tif.save(self.cell_image.astype(np.uint16), 
                     compress=9, extratags=[('ImageDescription','s',0,_nuc_xml,True)])
            tif.save(self.edge_image.astype(np.uint16), 
                     compress=9, extratags=[('ImageDescription','s',0,_mem_xml,True)])
            if processed_image:
                tif.save(self.processed_image.astype(np.uint16), 
                         compress=9, extratags=[('ImageDescription','s',0,_pro_xml,True)])
    def save_binary_seg_maps_r(self,path,processed_image=True):
        """
        Save a binary_seg_map image based on the **.make_cell_image_()** images in an R-readable format
        """
        with TiffWriter(path) as tif:
            tif.save(self.cell_image.astype(np.uint16), 
                     compress=9, description=_nuc_xml)
            tif.save(self.edge_image.astype(np.uint16), 
                     compress=9, description=_mem_xml)
            if processed_image:
                tif.save(self.processed_image.astype(np.uint16), 
                         compress=9, description=_pro_xml)
    def save_component(self,path,frame_name):
        """
        Save a component image
        """
        with TiffWriter(path) as tif:
            for channel in self.components:     
                tif.save(self.components[channel].astype(np.float16), 
                     compress=9, extratags=[('ImageDescription','s',0,_get_description(channel,frame_name),True)])

    def make_inform_frame(self,model_cells,base_path,sample_name,frame_name,r_format=False):
        """
        Save the inform 'cell_seg_data.txt', 'score_data.txt' and 'binary_seg_maps.tif' 
        to  **basepath/sample_name/**
        """
        self.set_cell_coordinates(model_cells)
        cell_seg = _construct_cell_seg(model_cells)
        score = _construct_score(model_cells)
        score.loc[:,'Sample Name'] = frame_name
        path = os.path.join(base_path,sample_name)
        if not os.path.exists(path):
            os.makedirs(path)
        cell_seg.to_csv(os.path.join(path,sample_name+'_'+frame_name+'_cell_seg_data.txt'),
            index=False,sep="\t")
        score.to_csv(os.path.join(path,sample_name+'_'+frame_name+'_score_data.txt'),
            index=False,sep="\t")
        self.make_cell_image()
        bfile = os.path.join(path,sample_name+'_'+frame_name+'_binary_seg_maps.tif')
        if r_format:
            self.save_binary_seg_maps_r(bfile,processed_image=True)
        else:
            self.save_binary_seg_maps(bfile,processed_image=True)
        return #path,cell_seg,score
        

def _construct_cell_seg(cells):
    header = ['Cell ID','Cell X Position','Cell Y Position','Nucleus Area (pixels)',
          'Nucleus Area (percent)','Nucleus Compactness','Nucleus Minor Axis',
          'Nucleus Major Axis','Entire Cell PD-1 (Opal 540) Mean (Normalized Counts, Total Weighting)',
          'Nucleus PD-1 (Opal 540) Mean (Normalized Counts, Total Weighting)',
          'Membrane PD-1 (Opal 540) Mean (Normalized Counts, Total Weighting)',
          'Membrane Area (pixels)','Membrane Area (percent)','Membrane Compactness',
          'Membrane Minor Axis','Membrane Major Axis',
          'Nucleus PD-Ligand-1 (Opal 690) Mean (Normalized Counts, Total Weighting)',
          'Entire Cell PD-Ligand-1 (Opal 690) Mean (Normalized Counts, Total Weighting)',
          'Membrane PD-Ligand-1 (Opal 690) Mean (Normalized Counts, Total Weighting)',
          'Entire Cell Area (pixels)',
          'Entire Cell Area (percent)','Entire Cell Compactness',
          'Entire Cell Minor Axis','Entire Cell Major Axis','Phenotype','Confidence']
    cs = cells.copy().rename(columns={'id':'Cell ID',
                                  'x':'Cell X Position',
                                  'y':'Cell Y Position',
                                  'phenotype_label':'Phenotype'})
    fill_in = [x for x in header if x not in cs.columns]
    for col in fill_in: cs[col] = 1
    cs.loc[cs['PD1']=='-',cs.columns.str.contains('PD-1')] = 0
    cs.loc[cs['PDL1']=='-',cs.columns.str.contains('PD-Ligand-1')] = 0
    return cs.drop(columns=['PD1','PDL1'])

def _construct_score(cells):
    header = ['Path','Sample Name','Tissue Category',
          'First Cell Compartment','First Stain Component',
          'Second Cell Compartment','Second Stain Component',
          'Double Negative','Single PD-Ligand-1 (Opal 690)',
          'Single PD-1 (Opal 540)','Double Positive',
          'Tissue Category Area (Percent)','Number of Cells',
          'PD-Ligand-1 (Opal 690) Threshold',
          'PD-1 (Opal 540) Threshold','Lab ID','Slide ID',
          'TMA Sector','TMA Row','TMA Column','TMA Field',
          'inForm 2.1.5430.24864']
    score = pd.DataFrame([len(header)*[0.5]],columns=header)
    score['Path'] = '/location'
    score['Sample Name'] = 'sample_name'
    score['Tissue Category'] = ''
    score['First Stain Component'] = 'PD-Ligand-1 (Opal 690)'
    score['Second Stain Component'] = 'PD-1 (Opal 540)'
    score['First Cell Compartment'] = 'Membrane'
    score['Second Cell Compartment'] = 'Nucleus'
    return score
def _get_description(channel_name,frame_name):
    istr = '<?xml version="1.0" encoding="utf-16"?>\r\n<PerkinElmer-QPI-ImageDescription>\r\n  <DescriptionVersion>2</DescriptionVersion>\r\n  <AcquisitionSoftware>Mantra</AcquisitionSoftware>\r\n  <ImageType>FullResolution</ImageType>\r\n  <Identifier>d512eea2-f3fd-4ba5-92ed-69b9c9aa9bf2</Identifier>\r\n  <SlideID>'+\
        frame_name+'</SlideID>\r\n  <Barcode />\r\n  <ComputerName>LAPTOP</ComputerName>\r\n  <IsUnmixedComponent>True</IsUnmixedComponent>\r\n  <ExposureTime>16969</ExposureTime>\r\n  <SignalUnits>33</SignalUnits>\r\n  <Name>'+\
        channel_name+'</Name>\r\n  <Color>255,255,0</Color>\r\n  <Objective>20x [20x]</Objective>\r\n  <ValidationCode>80D81BB521952BDCC3ECA60A218A344B5FA618CF97F07F9D7E462CFAC6F50BA0</ValidationCode>\r\n</PerkinElmer-QPI-ImageDescription>'
    return(istr)