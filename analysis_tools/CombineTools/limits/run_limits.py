import datacard_converter.DatacardConverter as DatacardConverter

if __name__ == "__main__":
    """
    run datacards for all masspoints with Combine
    """
    parser = argparse.ArgumentParser(description="run combine for datacards in given directory")
    parser.add_argument("-j", "--max_jobs", default=16, type=int,
                        help="maximum number of parallel jobs")

    args = parser.parse_args()
    max_jobs = args.max_jobs

    # Move this to JSON or YAML in the future
    # Settings
    datacard = "datacard_BstarToTW.txt"
    keywords = {"SHAPES":"shapes_BstarToTW.root"}
    combine_setup = "-M HybridNew --rMin -2 --rMax 2 -H AsymptoticLimits"
    masspoints = [700, 800, 900, 1000, 1100, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200]
    

    jobs = []
    converter = DatacardConverter(datacard, keywords)    
    for masspoint in masspoints:
        workspace = converter.convert_datacard(masspoint)
        combine_command = "combine {workspace}_{mass}.root {combine_setup}".format(**{"workspace":workspace,"mass":mass,"combine_setup":combine_setup})
        f = open("log/{}.log".format(workspace),'w')
        jobs.append(combine_command

    runner = CombineRunner(workspaces)

    # Plotter
    predicted = [32.8, 17.06, 9.46, 5.93, 3.23, 1.98, 0.81, 0.35, 0.17, 0.081, 0.041, 0.022, 0.012, 0.0067, 0.0038, 0.001861, 0.001064, 0.0006169, 0.0003632, 0.0002179, 0.0001324]

