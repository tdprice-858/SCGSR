import yaml

temps = [temp for temp in range(300,625,25)]
methyl_coverages = []
oxygen_coverages = []
file_dir = '/Users/tdprice/Desktop/Cantera test_catalytic combustion example'

for temp in temps:
    if temp == 600:
        file = f"{file_dir}/MeOH_ox_Ag_catalytic_combustion_{temp}.0.yaml"
    else:
        file = f"{file_dir}/MeOH_ox_Ag_catalytic_combustion_{temp}.yaml"
    with open(file, 'r') as file:
        sim_data = yaml.safe_load(file)
    print(sim_data['soln1']['surface']['coverages'])
    methyl_coverages.append(sim_data['soln1']['surface']['coverages']['C[Pt](15)'])
    oxygen_coverages.append(sim_data['soln1']['surface']['coverages']['OX(46)'])
print(sim_data['soln1'])
print(methyl_coverages)
print(oxygen_coverages)
