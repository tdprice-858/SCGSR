
path_to_species_dict = '/Users/tdprice/Desktop/Desktop/\
RMG-Py/examples/rmg/catalysis/MeOH_ox_Ag/chemkin/species_dictionary.txt'
gas_species = []
ad_species = []
with open(path_to_species_dict, 'r') as f:
    list = f.read().split('\n')
    length = len(list)
    for count,line in enumerate(list):
        if count==0:
            if 'X' in list[count + 1]:
                ad_species.append(list[count])
            else:
                gas_species.append(list[count])
        elif line=='' and (count+1)<length:
            if list[count+1]!='':
                if 'X' in list[count+1]:
                    ad_species.append(list[count+1])
                else:
                    gas_species.append(list[count+1])
        else:
            continue
print(gas_species)
print(ad_species)



