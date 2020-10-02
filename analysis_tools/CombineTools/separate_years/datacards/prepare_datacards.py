import datacard2workspace

masspoints = [700, 800, 900, 1000, 1100, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200]
datacard_name = "datacard_BstarToTW_all.txt"

if __name__ == "__main__":
    for masspoint in masspoints:
        datacard2workspace.convert_datacard(datacard_name, masspoint)
        
