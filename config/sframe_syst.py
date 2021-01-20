#!/usr/bin/env python

import os, sys, itertools

# sys.path.append('/nfs/dust/cms/user/tholenhe/installs/varial-stable/Varial')
sys.path.append('/nfs/dust/cms/user/froehlia/CMSSW_10_2_10/src/UHH2/Varial')

##################################### definition of UserConfig item changes ###

sys_uncerts = {
    # 'name' : {'item name': 'item value', ...},
     # 'NOMINAL'               : {'b_TopPtReweight':'true', 'b_PDFVariation':'true'}, #{...} is a dummy
     # 'SCALE_upup'            : {'ScaleVariationMuR':'up','ScaleVariationMuF':'up'},
     # 'SCALE_upnone'          : {'ScaleVariationMuR':'up','ScaleVariationMuF':'none'},
     # 'SCALE_noneup'          : {'ScaleVariationMuR':'none','ScaleVariationMuF':'up'},
     # 'SCALE_nonedown'        : {'ScaleVariationMuR':'none','ScaleVariationMuF':'down'},
     # 'SCALE_downnone'        : {'ScaleVariationMuR':'down','ScaleVariationMuF':'none'},
     # 'SCALE_downdown'        : {'ScaleVariationMuR':'down','ScaleVariationMuF':'down'},
     'PU_up'                 : {'Systematic_PU':'up'},
     'PU_down'               : {'Systematic_PU':'down'},
     # 'MUOID_up'               : {'Systematic_MuonID':'up'},
     # 'MUOID_down'             : {'Systematic_MuonID':'down'},
     # 'MUOTR_up'               : {'Systematic_MuonTrigger':'up'},
     # 'MUOTR_down'             : {'Systematic_MuonTrigger':'down'},
     # 'MUOISO_up'              : {'Systematic_MuonIso':'up'},
     # 'MUOISO_down'            : {'Systematic_MuonIso':'down'},
     'ELEID_up'               : {'Systematic_ElectronID':'up'},
     'ELEID_down'             : {'Systematic_ElectronID':'down'},
     'ELETR_up'               : {'Systematic_ElectronTr':'up'},
     'ELETR_down'             : {'Systematic_ElectronTr':'down'},
     'ELEREC_up'              : {'Systematic_ElectronReco':'up'},
     'ELEREC_down'            : {'Systematic_ElectronReco':'down'},
     'BTAG_bc_up'             : {'Systematic_BTag':'up_bc'},
     'BTAG_bc_down'           : {'Systematic_BTag':'down_bc'},
     'BTAG_udsg_up'           : {'Systematic_BTag':'up_udsg'},
     'BTAG_udsg_down'         : {'Systematic_BTag':'down_udsg'},
    'TOPTAG_merged_up'         : {'Systematic_TopTag':'up_merged'},
    'TOPTAG_merged_down'       : {'Systematic_TopTag':'down_merged'},
    'TOPTAG_semi_up'           : {'Systematic_TopTag':'up_semi'},
    'TOPTAG_semi_down'         : {'Systematic_TopTag':'down_semi'},
    'TOPTAG_non_up'            : {'Systematic_TopTag':'up_non'},
    'TOPTAG_non_down'          : {'Systematic_TopTag':'down_non'},
     'TOPPT_alpha_up'         : {'Systematic_TopPt_a':'up'},
     'TOPPT_alpha_down'       : {'Systematic_TopPt_a':'down'},
     'TOPPT_beta_up'          : {'Systematic_TopPt_b':'up'},
     'TOPPT_beta_down'        : {'Systematic_TopPt_b':'down'},
     'Prefiring_up'           : {'Systematic_Prefiring':'up'},
     'Prefiring_down'         : {'Systematic_Prefiring':'down'},

     # 'JEC_up'               : {'b_TopPtReweight':'true', 'b_PDFVariation':'false'},
     # 'JEC_down'             : {'b_TopPtReweight':'true', 'b_PDFVariation':'false'},
     # 'JER_up'               : {'b_TopPtReweight':'true', 'b_PDFVariation':'false'},
     # 'JER_down'             : {'b_TopPtReweight':'true', 'b_PDFVariation':'false'},
     }
start_all_parallel = False

############################################################### script code ###
import varial
import sys
import os

if len(sys.argv) != 2:
    print 'Plz. give me da name of da sframe-config! ... dude!'
    exit(-1)


def set_uncert_func(uncert_name):
    uncert = sys_uncerts[uncert_name]
    def do_set_uncert(element_tree):
        cycle = element_tree.getroot().find('Cycle')
        user_config = cycle.find('UserConfig')
        output_dir = cycle.get('OutputDirectory')

        cycle.set('OutputDirectory', os.path.join(output_dir, uncert_name+'/'))

        for name, value in uncert.iteritems():
            uc_item = list(i for i in user_config if i.get('Name') == name)
            assert uc_item, 'could not find item with name: %s' % name
            uc_item[0].set('Value', value)

    return do_set_uncert


from varial.extensions.sframe import SFrame
from varial import tools
if start_all_parallel:
    ToolChain = tools.ToolChainParallel
else:
    ToolChain = tools.ToolChain


class MySFrameBatch(SFrame):

    def configure(self):
        self.xml_doctype = self.xml_doctype + """
<!--
   <ConfigParse NEventsBreak="100000" FileSplit="0" AutoResubmit="0" />
   <ConfigSGE RAM ="2" DISK ="2" Mail="alexander.froehlich@desy.de" Notification="as" Workdir="workdir"/>
-->
"""
        if os.path.exists(self.cwd + 'workdir'):
            opt = ' -f'
        else:
            opt = ' -s'

        self.exe = 'sframe_batch.py' + opt


sframe_tools = ToolChain(
    'SFrameUncertsSR_Electron_18',
    list(
        MySFrameBatch(
            cfg_filename=sys.argv[1],
            xml_tree_callback=set_uncert_func(uncert),
            name='SFrame_' + uncert,
            halt_on_exception=False,
        )
        for uncert in sys_uncerts
    )
)


if __name__ == '__main__':
    varial.tools.Runner(sframe_tools)
