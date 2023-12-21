from pathlib import Path
from bladegen_files.run_case import main

if __name__ == "__main__":
    # Input
    bladegen_working_dir = Path.cwd() / 'bladegen_files'
    sc90c_working_dir = Path.cwd() / 'sc90c_files'
    turbogrid_working_dir = Path.cwd() / 'cfd/turbogrid_files'
    name = 'dt3a'
    n_stages = 1
    flg_igv = 1
    axis = 'x'  # select axis of rotation

    main(
        bladegen_working_dir=bladegen_working_dir, 
        sc90c_working_dir=sc90c_working_dir,
        turbogrid_working_dir=turbogrid_working_dir,
        name=name,
        n_stages=n_stages,
        flg_igv=flg_igv,
        log=False, # for debug parsing of SC90
        plot=True, # generate plots
        axis=axis,
    )
