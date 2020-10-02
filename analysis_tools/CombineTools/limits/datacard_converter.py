import argparse
import subprocess


class DatacardConverter(object):
    def __init__(self, datacard_name, keywords = {}):
        self.datacard_name = datacard_name
        self.keywords = keywords

    def convert_datacard(self, masspoint, workspace_name = ""):
        if not workspace_name:
            workspace_name = datacard_name.replace("datacard","workspace").replace(".txt","")
        keywords_string = " ".join(["--keyword-value {}={}".format(k,v) for k,v in self.keywords.iteritems()])
        command = "text2workspace.py {datacard} -m {mass} {keywords} -o {workspace}_{mass}.root".format(**{"datacard":self.datacard_name, "mass":masspoint, "workspace":workspace_name, "keywords":keywords_string})
        subprocess.Popen(command,shell=True)
        return workspace_name

