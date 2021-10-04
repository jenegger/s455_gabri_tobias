# s455_gabri_tobias

Data taken from:
/lustre/r3b/202103_s455/lmd

Reference files:
main0273_0010.lmd or main0273_0001.lmd

Philipp's script to haecksle the CALIFA readout blocks and stitch the events (4 us):
/u/land/lynx.landexp/202105_s494/scripts/califa_merge_stitch.bash

Upexps unpacker used:
/u/land/fake_cvmfs/9.13/upexps/202103_s455_jentob/202103_s455

Latest R3BRoot version (dev branch) used, commit 2ca8be537a39b47d265fc3eaed23fc52b006b85c (Thu Sep 9 00:09:26 2021 +0200)



(macro Philipp used to create mapping file: /u/land/klenze/r3bmap/run/califa/mapping/ucesb_mapping_gen_CALIFA_TOT.py)


The calibration file dummy_sofia.par is the same as /u/land/r3broot/202106_testing/R3BRoot_20210726/sofia/macros/s455Up2p/parameters/CalibParam.par


For further data checking we look at following lmd files:
-> main0208_0008.lmd (as it is one of the beginnings)
-> main0220_0004.lmd (this was after the lxana01 crash, do we see some shifts in WRTS?? this has to be investigated....) 
