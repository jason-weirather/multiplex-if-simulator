from multiplexifsimulator import FrameEmitter
from tifffile import TiffWriter
import numpy as np
_nuc_xml = '<?xml version="1.0" encoding="utf-16"?>\r\n<SegmentationImage>\r\n <Version>1</Version>\r\n <CompartmentType>Nucleus</CompartmentType>\r\n</SegmentationImage>'
_mem_xml = '<?xml version="1.0" encoding="utf-16"?>\r\n<SegmentationImage>\r\n <Version>1</Version>\r\n <CompartmentType>Membrane</CompartmentType>\r\n</SegmentationImage>'
_pro_xml = '<?xml version="1.0" encoding="utf-16"?>\r\n<ProcessRegionImage>\r\n <Version>1</Version>\r\n</ProcessRegionImage>'
class FrameEmitterInForm(FrameEmitter):
    def __init__(self,shape=(1040,1392),
                      cell_steps=17,
                      boundary_steps=10):
        super(FrameEmitterInForm, self).__init__(shape=shape,cell_steps=cell_steps,boundary_steps=boundary_steps)
        return        
    def save_binary_seg_maps(self,path,processed_image=True):
        with TiffWriter(path) as tif:
            tif.save(self.cell_image.astype(np.uint16), 
                     compress=9, extratags=[(270,'s',1,_nuc_xml,True)])
            tif.save(self.edge_image.astype(np.uint16), 
                     compress=9, extratags=[(270,'s',1,_mem_xml,True)])
            if processed_image:
                tif.save(self.processed_image.astype(np.uint16), 
                         compress=9, extratags=[(270,'s',1,_pro_xml,True)])
