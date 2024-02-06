#!/bin/sh

graphlan_annotate.py --annot Figures/Figure3/annot.txt Data/gtdbtk_bacteria.phy Data/bacterial_tree_annotated.xml

graphlan.py Data/bacterial_tree_annotated.xml Figures/Figure3/SciData_fig3.svg --dpi 300 --size 7