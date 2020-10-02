import argparse
import subprocess

def convert_datacard(datacard_name, masspoint):
    workspace_name = datacard_name.replace("datacard","workspace").replace(".txt","")
    command = "text2workspace.py {datacard} -m {mass} -o {workspace}_{mass}.root".format({"datacard":datacard_name, "mass":masspoint, "workspace":workspace_name})
    subprocess.Popen(command,shell=True)
