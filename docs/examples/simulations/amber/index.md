# Amber

[Amber](https://ambermd.org/) is a widely used biomolecular simulation package for decades.
Much of the simulation workflows are routine that can be automated and rendered on demand.

## Input files

```python
--8<-- "docs/examples/simulations/amber/scripts/min-input.py:5"
```

```text
--8<-- "docs/examples/simulations/amber/scripts/min.in:1:9"

    ...

--8<-- "docs/examples/simulations/amber/scripts/min.in:55"
```

## Command

```python
--8<-- "docs/examples/simulations/amber/scripts/amber-command.py:5"
```

```python
["pmemd.cuda -i min.in -o min.mdout -inf min.mdinfo -p mol.prmtop -c mol.inpcrd"]
```
