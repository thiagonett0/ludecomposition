# Decomposição LU (Gaxpy)

Esse repositório apresenta o método Decomposição LU, versão Gaxpy, apresentado na página 100, 3.2.5, de [1].


### Convenção 1

Matrizes são representadas por letras maiúsculas em negrito; vetores, por minúsculas.


### Convenção 2 (triangularidade inferior)

Representa-se matrizes triangulares inferiores por $\mathbf{L}$.


### Convenção 3 (triangularidade superior)

Representa-se matrizes triangulares superiores por $\mathbf{U}$.


### Decomposição LU


Seja uma matriz $\mathbf{A}$ tal que $\mathbf{A} \in \mathbb{R}^{n \times n}$. Se $\textrm{det}( \mathbf{A}(1:k, 1:k)) \neq 0 \,\,\, \forall k \in \{1, ... , n-1 \}$, então $\mathbf{A}$ possui fatoração LU. Adicionalmente, se $\mathbf{A}$ for não singular, então a fatoração LU é única e $\textrm{det}(A) = u_{11} ... u_{nn}$.


$$
\begin{equation}
    \begin{pmatrix}
        a_{11} & ... & a_{1n} \\
        \vdots & \ddots & \vdots \\
        a_{n1} & ... & a_{nn}
    \end{pmatrix}
    =
    \begin{pmatrix}
        l_{11} & ... & 0 \\
        \vdots & \ddots & \vdots \\
        l_{n1} & ... & l_{nn}
    \end{pmatrix}
    \begin{pmatrix}
        u_{11} & ... & u_{1n} \\
        \vdots & \ddots & \vdots \\
        0 & ... & u_{nn}
    \end{pmatrix}
\end{equation}
$$


**Demonstração:**


Teorema 3.2.1, página 97 de [1].


# Referências

[1] Matrix Computations; GOLUB, Gene H., VAN LOAN, Charles F.; 3ed, 1996. ISBN: 0-8018-5414-8