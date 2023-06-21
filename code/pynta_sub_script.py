#!/bin/bash -l
#SBATCH --ntasks=24
#SBATCH --time=21-00:00:00
#SBATCH --partition=month-long-cpu
#SBATCH --account=shikim
#SBATCH --job-name=FW_run
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.error

from pynta.main import Pynta

pyn = Pynta(path="/home/shikim/pynta-production/han",launchpad_path="/home/shikim/local_launchpad.yaml",
                fworker_path="/home/shikim/my_fworker.yaml",
                queue_adapter_path="/home/shikim/local_qadapter.yaml",
                rxns_file="/home/shikim/pynta-production/han/han.yaml",
                surface_type="fcc111",metal="Cu",socket=False,queue=True,njobs_queue=5,
                repeats=(3,3,4),label="han",num_jobs=1,
                software_kwargs={'kpts': (3, 3, 4), 'tprnfor': True, 'occupations': 'smearing',
                            'smearing':  'marzari-vanderbilt', 'input_dft': 'BEEF-vdW',
                            'degauss': 0.01, 'ecutwfc': 40, 'nosym': True,
                            'conv_thr': 1e-6, 'mixing_mode': 'local-TF',
                            "pseudopotentials": {"Cu": 'Cu.pbe-spn-kjpaw_psl.1.0.0.UPF',
                            "H": 'H.pbe-kjpaw_psl.1.0.0.UPF',
                            "O": 'O.pbe-n-kjpaw_psl.1.0.0.UPF',
                            'N': 'N.pbe-n-kjpaw_psl.1.0.0.UPF'
                                },
                            "command": 'srun /opt/custom/espresso/6.6_nostress/bin/pw.x < PREFIX.pwi > PREFIX.pwo'},
                software_kwargs_gas={'kpts': 'gamma', 'tprnfor': True, 'occupations': 'smearing',
                            'smearing':  'gauss', 'input_dft': 'BEEF-vdW',
                            'degauss': 0.005, 'ecutwfc': 40, 'nosym': True,
                            'conv_thr': 1e-6, 'mixing_mode': 'local-TF',
                            'mixing_beta': 0.2, 'mixing_ndim': 10,
                            "pseudopotentials": {"Cu": 'Cu.pbe-spn-kjpaw_psl.1.0.0.UPF',
                            "H": 'H.pbe-kjpaw_psl.1.0.0.UPF',
                            "O": 'O.pbe-n-kjpaw_psl.1.0.0.UPF',
                            'N': 'N.pbe-n-kjpaw_psl.1.0.0.UPF'
                                },
                            "command": 'mpirun -np 1 /home/shikim/qe-7.1/bin/pw.x -inp < PREFIX.pwi > PREFIX.pwo'},
                slab_path="/home/shikim/pynta-production/han/slab.xyz",
               TS_opt_software_kwargs={"conv_thr":1e-11},
               )
pyn.generate_slab()