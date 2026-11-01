#!/bin/zsh
# E5: lighter load (1024 B every 20 ms) to isolate the source-buffer flush component
cd ~/ns-3.48-nr/results/sys/raw
B=~/ns-3.48-nr/build/scratch/ns3.48-nr-udn-handover-v2-optimized
for r in 1 2 3 4 5; do
  for pol in a3 og; do
    pre=E5_isd50_tm40_mix_$( [ $pol = a3 ] && echo a3t100 || echo og-5.5 )_pkt20_r$r
    [ -f $pre-summary.csv ] && continue
    $B --gridN=4 --simTime=20 --warm=2 --run=$r --isd=50 --tm=40 --speeds=3,30,60,120 --policy=$pol --ttt=100 --guardThr=-5.5 --pktInterval=0.02 --prefix=$pre > /dev/null 2>&1
  done
done
