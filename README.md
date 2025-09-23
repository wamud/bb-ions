# bb-ions

This repository is for simulating an ion trap realisation of a Bivariate Bicycle (BB) code [[1]](https://www.nature.com/articles/s41586-024-07107-7.pdf).
The ion trap architecture comprises data qubit modules and check qubit modules as per [[2]](https://arxiv.org/pdf/2508.01879).
When a data qubit module is aligned with a check qubit module, forming a pair, their qubits can interact.
The data qubit modules are fixed while the check qubit modules are equipped with a cyclic shift.
This enables all required interactions for the syndrome extraction circuits of the BB code, which we implement using algorithm 2 of [[2]](https://arxiv.org/pdf/2508.01879).
This implements X checks first then Z checks.

## Bivariate Bicycle Code

 A BB code [[1]](https://www.nature.com/articles/s41586-024-07107-7.pdf) is defined by the parity check matrices

$$H_x = [A|B] \\
H_z = [B^T|A^T]$$

where $A$ and $B$ are (matrix) polynomials of the variables

$$x = S_l \otimes \mathbb{1}_m\\
y = \mathbb{1}_l \otimes S_m$$

where $S_j$ is the cyclic matrix of dimension $j\times j$, e.g.

$$
S_3 =
\begin{bmatrix}
0 & 1 & 0 \\
0 & 0 & 1 \\
1 & 0 & 0
\end{bmatrix}
$$
