"""
Descrição: código que executa soluciona qualquer sistema Ax=b para A quadrada e não singular
    por meio da decomposição LU.
"""

from random import uniform



def zeros(n):
    Z = [[0 for _ in range(n)] for _ in range(n)]
    return Z


def identidade(n):
    I = [[0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        I[i][i] = 1
    return I


def produtoescalar(a, u):
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


def triBack(U, b):

    """
    Descrição: ..

    Entrada(s):
                i) ..

    Saída(s):
                i) ..
    """

    n = len(U)-1
    b[n] = b[n]/U[n][n]
    for i in range(1, n+1):
        b[n-i] = (b[n-i] - produtoponto(U[n-i][n-i+1:], b[n-i+1:]))/U[n-i][n-i]
    return b


def transposta(A):
    At = [[None for _ in range(len(A[0]))] for _ in range(len(A))]
    for i in range(len(A)):
        for j in range(len(A[0])):
            At[i][j] = A[j][i]
    return At


def luDecomp(A):
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


A = [[1, 0, 0], [0, 4, 0], [0, -2, 0]]
B = identidade(3)
luDecomp(A)


def produtoMatricial(A, B):
    P = list()
    B = transposta(B)
    for linhas in A:
        Plinhas = list()
        for colunas in B:
            Plinhas.append(produtoponto(linhas, colunas))
        P.append(Plinhas)
    return P


def verificaLU(A):
    lu = luDecomp(A)
    L, U = lu[0], lu[1]
    P = produtoMatricial(L, U)
    S = zeros(len(A))
    for i in range(len(A)):
        for j in range(len(A[0])):
            S[i][j] = abs(P[i][j]-A[i][j])
            if abs(S[i][j]) > 1e-3:
                return False
    return True


def matrizAleatoria(tamanho):
    A = zeros(tamanho)
    for i in range(tamanho):
        for j in range(tamanho):
            A[i][j] = uniform(0, 1)
    return A


def verificador(tamanho):
    for n in range(2000):
        A = matrizAleatoria(tamanho)
        v = verificaLU(A)
        if v == False:
            return False
    return True


print(verificador(10))
