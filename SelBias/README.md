## Minimum doc on execution

usage: selbias.py [-h] -t GENETREE [-w WDIR] [--gcode GCODE]
                  (--mxmo | --hyptest H0|H1 | --freemodel)
                  (--codaln CODALN | --seq protaln nucseq) [--codeml CODEML]
                  [--debug] [--alnformat ALNFORMAT] [--pal2nal_args PAL2NAL]
                  [--cpu CORES]

SelBias

optional arguments:
  -h, --help            show this help message and exit
  -t GENETREE, --gtree GENETREE
                        Either the filename or the newick of the genetree
  -w WDIR, --wdir WDIR  Working directory to use. Default will be a tmp folder
                        created in the current directory
  --gcode GCODE         Genetic code to use
  --mxmo                Compute MxMo as described in the paper
  --hyptest H0|H1       The two models to test. List of valid pair are:
                        bsA1|bsA, M0|b_free, b_neut|b_free, M3|bsD, M1|bsA.
                        Just use 'b_neut|b_free') or avoid this option.
  --freemodel           Free ratio on each branch of the tree
  --codaln CODALN, -c CODALN
                        Codon alignment
  --seq protaln nucseq, -s protaln nucseq
                        nucleotide sequence and protein alignment as an
                        alternative to codon alignment
  --codeml CODEML       Location of the codeml binary. It is easier to add it
                        to your home path
  --debug               Print debug information
  --alnformat ALNFORMAT, -f ALNFORMAT
                        Alignment format to read. Default is paml, which is
                        the format returned by pal2nal
  --pal2nal_args PAL2NAL
                        Supplementary arguments for pal2nal
  --cpu CORES           Maximum number of CPU cores available in the execution
                        host. If higher than 1, multi-threading will be
                        enabled (if 0 all available cores will be used)
