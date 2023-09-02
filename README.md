# srtree-opt - a CLI tool to (re)optimize the numeric parameters of symbolic regression expressions

Optimize SR expressions using Levenberg-Marquardt.

## Instalation

```bash
stack install
```

## Command line help:

```bash
srtree-opt -h
```

## Example usage:

Parse the Heuristic Lab expressions in `test/feynman-III-17-37_expressions_heuristiclab.txt` file using the training data in `test/feynman_III_17_37_train.csv` (with variable names in the header) taking the first $50$ rows as the training set and the next rows as validation, select the columns $0,1,2$ as the input variables and runs a maximum of $10$ iterations. It also simplifies the expressions using equality saturation prior to optimizing.

```bash
srtree-opt -f hl -i test/feynman-III-17-37_expressions_heuristiclab.txt -d test/feynman_III_17_37_train.csv -r 50 -c 0,1,2 --niter 10 --hasheader --simplify
```

Same as above but using column $4$ as the target variable.

```bash
srtree-opt -f hl -i test/feynman-III-17-37_expressions_heuristiclab.txt -d test/feynman_III_17_37_train.csv -r 50 -c 0,1,2 -t 4 --niter 10 --hasheader --simplify
```

By default the fitted expressions will be printed in stdout and the stats (sse of training and validation prior to optimization and after optimization) will be printed to stderr. You can specify the output files with the `-o, -s` flags.
