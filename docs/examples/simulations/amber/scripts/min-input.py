import os

os.chdir(os.path.dirname(os.path.abspath(__file__)))

from atomea.containers.workflow.amber import Amber22Inputs

amber_inputs = Amber22Inputs()
amber_inputs.imin = 1
amber_inputs.ntx = 1
amber_inputs.irest = 0
amber_inputs.ntmin = 1
amber_inputs.maxcyc = 5000
amber_inputs.ncyc = 1000
amber_inputs.ntr = 1
amber_inputs.restraint_wt = 2.0
amber_inputs.restraintmask = "!(:WAT) & (@C,CA,N,O,O5',P,O3',C3',C4',C5')"
amber_inputs.ntb = 1
amber_inputs.ntf = 1
amber_inputs.ntc = 1
amber_inputs.cut = 10.0
amber_inputs.ntxo = 2
amber_inputs.ntwr = 200
amber_inputs.ntpr = 1
amber_inputs.ntwx = 200
amber_inputs.ioutfm = 1
amber_inputs.iwrap = 1

amber_inputs.write_render("min.in")
