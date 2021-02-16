# stabilizer-circuits-normal-forms

This source code is the C implementation with a simple text-based user interface of an algorithm that writes stabilizer quantum circuits under normal form.

The process is described in a paper at : https://arxiv.org/abs/2012.09224.

To compile the code just use the 'make' command : this creates the 'stabnf' command.

The 'manual' mode of the command reproduces the induction steps of the  algorithm.

The 'statistics' mode of the command allows to test the efficiency of the algorithm to reduce the gate count of a circuit.

Let n be the number of qubits of the circuits (2 <= n <= 1000).

To launch the manual mode : ./stabnf n man

To launch the statistics mode : ./stabnf n stat





