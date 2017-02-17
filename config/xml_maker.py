import csv, fileinput

datasetPath = '/nfs/dust/cms/user/froehlia/CMSSW_8_0_24_patch1/src/UHH2/common/datasets/RunII_80X_v3/'

with open('Background.csv', 'r') as f:
  reader = csv.reader(f)
  samplelist = list(reader)


with open('Base.xml', 'r') as xmlfile:
    with open('New.xml', 'w') as newfile:
        for line in xmlfile:
            newfile.write(line)
            if '<!-- Underground -->' in line:
                for row in samplelist:
                    s = '<!ENTITY ' + row[0] + ' SYSTEM "' + datasetPath + row[0] + '.xml'
                    s = s + '">\n'
                    newfile.write(s)
            if '<!-- InputData for Underground -->' in line:
                for row in samplelist:
                    s1 = '<InputData Lumi="' + str(int(float(row[1]))) + '" NEventsMax="&NEvents;" Type="MC" Version="' + row[0] + '" Cacheable="&cacheable;">\n'
                    s2 = '&' + row[0] + ';\n'
                    s3 = '<InputTree Name="AnalysisTree" /> \n <OutputTree Name="AnalysisTree" /> \n </InputData> \n'
                    newfile.write(s1+s2+s3+'\n')
