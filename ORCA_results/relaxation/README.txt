Relaxation of octylsilanetriol

Gas-phase geometry optimization and vibrational frequencies of octylsilanetriol, Si(OH)3-C8H17, at the
HF/6-31G* level with ORCA. The start structure is a relaxed 9-carbon molecule without its last CH2 group.

Files
octylsilanetriol_start.xyz       start structure
octylsilanetriol.inp             ORCA input
octylsilanetriol.out             ORCA output
octylsilanetriol.xyz             relaxed structure
octylsilanetriol_trj.xyz         structures of the optimization steps
octylsilanetriol.hess            Hessian of the frequency calculation
octylsilanetriol.engrad          energy and gradient of the last optimization step
octylsilanetriol_property.txt    properties, as written by ORCA
octylsilanetriol.opt             data of the optimizer, as written by ORCA
octylsilanetriol.densities       electron densities, as written by ORCA
octylsilanetriol.err             error messages of ORCA (empty)
extract_bonds_and_angles.py      writes the bond lengths and angles of the relaxed structure to
                                 bonds_and_angles.txt and prints the Si-C1 bond
bonds_and_angles.txt             bond lengths (Angstrom) and angles (degrees) of the relaxed structure,
                                 as printed by ORCA

The .gbw file, the binary wavefunction of ORCA, is not included.

Run
orca octylsilanetriol.inp > octylsilanetriol.out
python3 extract_bonds_and_angles.py
