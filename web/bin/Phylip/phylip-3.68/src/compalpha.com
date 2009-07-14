$ define lnk$library sys$library:vaxcrtl
$ write sys$output "Building PHYLIP version 3.6"
$ write sys$output "phylip"
$ cc phylip.c
$ write sys$output "seq"
$ cc seq.c
$ write sys$output "disc"
$ cc disc.c
$ write sys$output "dist"
$ cc dist.c
$ write sys$output "cont"
$ cc cont.c
$ write sys$output "moves"
$ cc moves.c
$ write sys$output "contml"
$ cc contml.c
$ link contml.obj,cont.obj,phylip.obj
$ write sys$output "clique"
$ cc clique.c
$ link clique.obj,disc.obj,phylip.obj
$ write sys$output "contrast"
$ cc contrast.c
$ link contrast.obj,cont.obj,phylip.obj
$ write sys$output "dnadist"
$ cc dnadist.c
$ link dnadist.obj,seq.obj,phylip.obj
$ write sys$output "dnainvar"
$ cc dnainvar.c
$ link dnainvar.obj,seq.obj,phylip.obj
$ write sys$output "dnaml"
$ cc dnaml.c
$ link dnaml.obj,seq.obj,phylip.obj
$ write sys$output "dnapars"
$ cc dnapars.c
$ link dnapars.obj,seq.obj,phylip.obj
$ write sys$output "draw"
$ cc draw.c
$ write sys$output "draw2"
$ cc draw2.c
$ write sys$output "drawgram"
$ cc drawgram.c
$ link draw.obj,draw2.obj,drawgram.obj,phylip.obj
$ write sys$output "drawtree"
$ cc drawtree.c
$ link draw.obj,draw2.obj,drawtree.obj,phylip.obj
$ write sys$output "factor"
$ cc factor.c
$ link factor.obj,phylip.obj
$ write sys$output "fitch"
$ cc fitch.c
$ link fitch.obj,dist.obj,phylip.obj
$ write sys$output "gendist"
$ cc gendist.c
$ link gendist.obj,phylip.obj
$ write sys$output "kitsch"
$ cc kitsch.c
$ link kitsch.obj,dist.obj,phylip.obj
$ write sys$output "neighbor"
$ cc neighbor.c
$ link neighbor.obj,dist.obj,phylip.obj
$ write sys$output "protdist"
$ cc protdist.c
$ link protdist.obj,seq.obj,phylip.obj
$ write sys$output "restdist"
$ cc restdist.c
$ link restdist.obj,seq.obj,phylip.obj
$ write sys$output "restml"
$ cc restml.c
$ link restml.obj,seq.obj,phylip.obj
$ write sys$output "retree"
$ cc retree.c
$ link retree.obj,moves.obj,phylip.obj
$ write sys$output "seqboot"
$ cc seqboot.c
$ link seqboot.obj,seq.obj,phylip.obj
$ write sys$output "Installing PHYLIP v3.6 binaries in -.-.exe"
$ copy *.exe;* [-.-.exe]
$ write sys$output "Installing font files in -.-.exe"
$ copy font*.;* [-.-.exe]
$ write sys$output "Installing documentation files in -.-.doc"
$ copy *.doc;* [-.-.doc]
$ write sys$output "Removing object files to save space"
$ delete *.obj;*
$ write sys$output "Removing executables from this directory"
$ delete *.exe;*
$ write sys$output "Finished cleanup."
$ write sys$output "Finished installation.
