"""
Descrição: código que executa soluciona qualquer sistema Ax=b para A quadrada e não singular
    por meio da decomposição LU.
"""


def zeros(n):

    """
    Descrição: função que cria uma matriz quadrada com entradas nulas;

    Entrada(s):
                i) n (int): tamanho da matriz;
        
    Saída(s):
                i) Z (list): matriz nula.
    """

    Z = [[0 for _ in range(n)] for _ in range(n)]
    return Z


def identidade(n):

    """
    Descrição: função que cria uma matriz identidade quadrada;

    Entrada(s):
                i) n (int): tamanho da matriz;
        
    Saída(s):
                i) I (list): matriz identidade.
    """

    I = [[0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        I[i][i] = 1
    return I


def produtoescalar(a, u):

    """
    Descrição: opera produto escalar entre um escalar (número real) e um vetor real. Para tal,
                basta acessar cada entrada do vetor, multiplicar pelo escalar e adicionar o 
                resultado em outra lista;

    Entrada(s):
                i) a (float): escalar real;
                ii) u (list): vetor real;
    
    Saída(s):
                i) p (list): vetor resultado.
    """

    p = list()
    for i in u:
        p.append(a*i)
    return p


def produtoponto(u, v):

    """
    Descrição: opera produto escalar entre u e v. Por meio de zip, acessa as entradas de ambos os vetores
                e adiciona o produto entre ambas entradas na variável pe;

    Entrada(s):
                i) u (list): vetor operador;
                ii) v (list): vetor operado;

    Saída(s):
                i) pe (float): produto escalar entre u e v.
    """

    pe = 0
    for i, j in zip(u, v):
        pe += i*j
    return pe


def triForward(L, b):

    """
    Descrição: função que reescreve a solução do sistema triangular inferior no vetor independente;

    Entrada(s):
                i) L (list): matriz triangular inferior;
                ii) b (list): vetor independente;

    Saída(s):
                i) b (list): vetor solução.
    """

    b[0] = b[0]/L[0][0]
    for i in range(1, len(L)):
        b[i] = (b[i] - produtoponto(L[i][:i], b[:i]))/L[i][i]
    return b


def triBack(U, b):

    """
    Descrição: função que reescreve a solução do sistema triangular superior no vetor independente;

    Entrada(s):
                i) U (list): matriz triangular superior;
                ii) b (list): vetor independente;

    Saída(s):
                i) b (list): vetor solução.
    """

    n = len(U)-1
    b[n] = b[n]/U[n][n]
    for i in range(1, n+1):
        b[n-i] = (b[n-i] - produtoponto(U[n-i][n-i+1:], b[n-i+1:]))/U[n-i][n-i]
    return b


def transposta(A):

    """
    Descrição: função que calcula a matriz transposta;

    Entrada(s):
                i) A (list): matriz de entrada;

    Saída(s):
                i) At (list): matriz transposta.
    """

    At = [[None for _ in range(len(A))] for _ in range(len(A[0]))]
    for i in range(len(A[0])):
        for j in range(len(A)):
            At[i][j] = A[j][i]
    print(At)
    return At


def luDecomp(A):

    """
    Descrição: função que calcula a decomposição LU de qualquer matriz quadrada;

    Entrada(s):
                i) A (list): matriz de entrada;

    Saída(s):
                i) L (list): matriz triangular inferior;
                ii) U (list): matriz triangular superior.
    """

    L, U = identidade(len(A)), zeros(len(A))
    v = [0 for _ in range(len(A))]
    for i in range(len(A)):
        if i == 0:
            for alpha in range(i, len(A)):
                v[alpha] = A[alpha][i]
        else:
            z = [0 for _ in range(i)]
            for j in range(i):
                z[j] = (A[j][i] - produtoponto(L[j][:j], z[:i]))/L[i][i]
                U[j][i] = z[j]
            for beta in range(i, len(A)):
                v[beta] = A[beta][i] - produtoponto(L[beta][:i], z)
        if i < len(A)-1:
            for gamma in range(i+1, len(A)):
                L[gamma][i] = v[gamma]/v[i]
        U[i][i] = v[i]
    return L, U


def solucaoLU(A, b):

    """
    Descrição: função que calcula a solução de qualquer sistema quadrado Ax=b;

    Entrada(s):
                i) A (list): matriz do sistema;
                ii) b (list): vetor independente;

    Saída(s):
                i) x (list): solução do sistema.
    """

    lu = luDecomp(A)
    L, U = lu[0], lu[1]
    y = triForward(L, b)
    x = triBack(U, y)
    return x
