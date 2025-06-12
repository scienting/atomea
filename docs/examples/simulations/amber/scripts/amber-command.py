import os

os.chdir(os.path.dirname(os.path.abspath(__file__)))

from atomea.project.workflow.amber import Amber22CLI

amber_cli = Amber22CLI(
    mdin="min.in",
    mdout="min.mdout",
    mdinfo="min.mdinfo",
    prmtop="mol.prmtop",
    inpcrd="mol.inpcrd",
)
amber_cli.module = "pmemd.cuda"
amber_command = amber_cli.render()

print(amber_command)
