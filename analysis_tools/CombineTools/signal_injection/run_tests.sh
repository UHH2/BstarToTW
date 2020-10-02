
# m = 1000; r = 1.0
mass=1000
r=1.0
text2workspace.py datacard_BstarToTW.txt -m ${mass} --keyword-value SHAPES=shapes_BstarToTW.root -o workspace_BstarToTW_${mass}.root
combine workspace_BstarToTW_${mass}.root -M FitDiagnostics -t 500 --expectSignal ${r} --toysFrequentist --bypassFrequentistFit --rMin -2 --rMax 2 --cminDefaultMinimizerStrategy 0
mv fitDiagnostics.root results/singal_injection_m_${mass}_r_injected.root
combine workspace_BstarToTW_${mass}.root -M FitDiagnostics -t 500 --expectSignal 0.0 --toysFrequentist --bypassFrequentistFit --rMin -2 --rMax 2 --cminDefaultMinimizerStrategy 0
mv fitDiagnostics.root results/singal_injection_m_${mass}_no_r.root

mass=2000
r=0.1
text2workspace.py datacard_BstarToTW.txt -m ${mass} --keyword-value SHAPES=shapes_BstarToTW.root -o workspace_BstarToTW_${mass}.root
combine workspace_BstarToTW_${mass}.root -M FitDiagnostics -t 500 --expectSignal ${r} --toysFrequentist --bypassFrequentistFit --rMin -2 --rMax 2 --cminDefaultMinimizerStrategy 0
mv fitDiagnostics.root results/singal_injection_m_${mass}_r_injected.root
combine workspace_BstarToTW_${mass}.root -M FitDiagnostics -t 500 --expectSignal 0.0 --toysFrequentist --bypassFrequentistFit --rMin -2 --rMax 2 --cminDefaultMinimizerStrategy 0
mv fitDiagnostics.root results/singal_injection_m_${mass}_no_r.root

mass=3000
r=0.05
text2workspace.py datacard_BstarToTW.txt -m ${mass} --keyword-value SHAPES=shapes_BstarToTW.root -o workspace_BstarToTW_${mass}.root
combine workspace_BstarToTW_${mass}.root -M FitDiagnostics -t 500 --expectSignal ${r} --toysFrequentist --bypassFrequentistFit --rMin -2 --rMax 2 --cminDefaultMinimizerStrategy 0
mv fitDiagnostics.root results/singal_injection_m_${mass}_r_injected.root
combine workspace_BstarToTW_${mass}.root -M FitDiagnostics -t 500 --expectSignal 0.0 --toysFrequentist --bypassFrequentistFit --rMin -2 --rMax 2 --cminDefaultMinimizerStrategy 0
mv fitDiagnostics.root results/singal_injection_m_${mass}_no_r.root

