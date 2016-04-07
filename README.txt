Prommute and Prommute Spacer program written by Jean O'Micks, April 5, 2016
The program was adapted from the C# version described in 
Cserháti M (2012) Prommute – A Promoter Mutation Simulation for Modeling the Evolution of Genetic Regulatory Elements. J Comput Sci Syst Biol 5:074-080.

The basic version can be run this way:

perl prommute.pl

Usage: prommute.pl
        -l promoter length
        -norg number of organisms
        -ncyc number of cycles
        -m matrix file
        -c motif cutoff
        -i list of tfs
        -o output file

Here the matrix file is tfmx_phiSITE.dat which contains all the matrix values for the 24 TFBS from the SCPD.
The TFBS are in order:

1 RLM1
2 SMP1
3 ABF1
4 CSRE
5 SCB
6 GAL4
7 GCN4
8 GCR1
9 MCB
10 MCM1
11 MATalpha2
12 MIG1
13 PHO2
14 PDR1_PDR3
15 PHO4
16 REB1
17 ROX1
18 RAP1
19 repressor_of_CAR1
20 SWI5
21 STE12
22 TBP
23 UASPHR
24 XBP1

This is important when listing the individual TFBS for the list of tfs

In order to run the program with spacers, use the prommute_spacer.pl variant:

perl prommute_spacer.pl

Usage: prommute_spacer.pl
        -l promoter length
        -norg number of organisms
        -ncyc number of cycles
        -m matrix file
        -c motif cutoff
        -s spacer length
        -i list of tfs
        -o output file

