# Decomposição LU (Gaxpy)

Esse repositório apresenta o método Decomposição LU, versão Gaxpy, apresentado na página 100, 3.2.5, de [1].


### Convenção 1

Matrizes são representadas por letras maiúsculas em negrito; vetores, por minúsculas.


### Definição 1 (Simetria)

Seja uma matriz quadrada $\mathbf{A} \in \mathbb{R}^{n \times n}$ com entradas $a_{i j}$. $\mathbf{A}$ é dita simétrica se e somente se satisfizer

$$
\begin{equation}
    a_{i j} = a_{j i} \textrm{    } \forall i, j = 0, 1, \ldots, n-1
\end{equation}
$$


### Definição 2 (Positiva-Definida)

Seja uma matriz $\mathbf{A} \in \mathbb{R}^{m \times n}$. $\mathbf{A}$ é dita positiva-definida se e somente se satisfizer


$$
\begin{equation}
    v^{T} \mathbf{A} v > 0 \textrm{    } \forall v \in \mathbb{R}^{m \times 1}
\end{equation}
$$


### Decomposição LU


Seja uma matriz $\mathbf{A}$ tal que $\mathbf{A} \in \mathbb{R}^{n \times n}$. Se $\textrm{det}( \mathbf{A}(1:k, 1:k)) \neq 0 \forall k \in \{1, \hdots, n-1 \}$, então $\mathbf{A}$ possui fatoração LU. Adicionalmente, se $\mathbb{A}$ for não singular, então a fatoração LU é única e $\textrm{det}(A) = u_{11} \hdots u_{nn}$.


$$
\begin{equation}
    x_i^{(k+1)} = \dfrac{b_i + \displaystyle\sum_{j=1}^{i-1} \left( a_{i j} x_j^{(k+1)} \right) - \displaystyle\sum_{j=i+1}^{n} \left( a_{i j} x_j^{(k)} \right)}{a_{i i}}
\end{equation}
$$


**Demonstração:**


Teorema 3.2.1, página 97 de [1].


# Referências

[1] Matrix Computations; GOLUB, Gene H., VAN LOAN, Charles F.; 3ed, 1996. ISBN: 0-8018-5414-8