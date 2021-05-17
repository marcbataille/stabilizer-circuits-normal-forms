# stabilizer-circuits-normal-forms, graph states branch

This source code is the C implementation with a simple text-based user interface of an algorithm that writes stabilizer quantum circuits in normal form.

To compile the code just use the 'make' command : this creates the 'stabnf' command.

The 'manual' mode of the command reproduces the induction steps of the  algorithm.

The 'statistics' mode of the command gives simple statistics to test the efficiency of a new implementation of graph states.

Let n be the number of qubits of the circuits (2 <= n <= 1000).

To launch the manual mode : ./stabnf n man

To launch the statistics mode : ./stabnf n stat





