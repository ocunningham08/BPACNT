# Bisphenol A (BPA) on Carbon Nanotubes (CNTs)

* ## cnt

A directory containing all of the LAMMPS MD files needed for simulating BPA on a CNT.

Individiual molecule files are provided for both BPA and Water, an input.data file containing the CNT geometry, a params.lmp file containing all the force field paramters, and the cnt_in.lmp file to insert the BPA and water and run the simulation.

* ## graphene

A directory containg the LAMMPS MD files for 60 BPA on a graphene sheet.

This is the same as the cnt directory except the BPA and water are already present in the input.data file with the graphene.

* ## MetaD

This contains the plumed files to run both Well-Tempered Metadynamics (-wt.dat) and Parallel Bias Metadynamics with Partitioned families (-pb.dat) as well as for a 2D Well Tempered run (-2d.dat).

These can be implementd by adding the line below to the LAMMPS input file

"fix pl all plumed plumedfile plumed.dat outfile p.log"

Relevant reweighting files can also be found in this directory.

# External Codes
The following codes have been used in this work...
* [LAMMPS](https://www.lammps.org/#gsc.tab=0)
* [PLUMED](https://www.plumed.org/doc-v2.9/user-doc/html/index.html)
* [moltemplate](https://moltemplate.org/)



