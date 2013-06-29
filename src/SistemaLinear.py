# -*- coding: iso-8859-1 -*-
"""
Autores: Bismarck Gomes Souza J�nior e Fernando Vizeu Santos
Vers�o: 1.0
Data: 27/07/2011

    Esse script cont�m quatro fun��es para resolver um sistema linear. Uma pelo
m�todo direto e outras tr�s por m�todos indiretos, recebendo somente as diago-
nais da matriz esparsa penta-diagonal.
"""
from numpy import array, append
from numpy import zeros, empty

def M5D2M(D1, D2, D3, D4, D5):
    """M5D2M (Matriz Penta-Diagonal para Matriz) recebe as 5 diagonais e cons-
    troi a matriz espar�a penta-diagonal. Retorna a matriz (array).
    
    Sistema:    A*x = b,   A[n][n] e b[n**2]
        
    Argumentos:
        D1:    menor diagonal inferior      (array: n**2-n)
        D2:    diagonal inferior            (array: n**2-1)
        D3:    diagonal principal           (array: n**2)
        D4:    diagonal superior            (array: n**2-1)
        D5:    maior diagonal superior      (array: n**2-n)
        
        *n: n�mero de linhas e colunas da matriz espar�a
    """
    N = D3.size
    n = int(N**0.5)
    M = zeros((N,N))
    M[0][0], M[0][1], M[0][n] = D3[0], D4[0], D5[0]
    for i in range(1, n):
        M[i][i-1], M[i][i], M[i][i+1], M[i][i+n] = D2[i-1], D3[i], D4[i], D5[i]
    for i in range(n, N-n):
        M[i][i-n], M[i][i-1], M[i][i], M[i][i+1], M[i][i+n] = D1[i-n], D2[i-1], D3[i], D4[i], D5[i]
    for i in range(N-n, N-1):
        M[i][i-n], M[i][i-1], M[i][i], M[i][i+1] = D1[i-n], D2[i-1], D3[i], D4[i]
    M[N-1][N-n], M[N-1][N-2], M[N-1][N-1] = D1[N-n-1], D2[N-2], D3[N-1]
    return M

def SSLD(A, b):
    """Solu��o do Sistema Linear Diretamente pelo m�todo da Elimina��o Gaussia-
    na. Retorna um array com a solu��o (x).
    
    Sistema:    A*x = b,   A[n][n] e b[n**2]
    
    Argumentos:
        A:     matriz dos coeficientes (array)
        b:     vetor dos valores       (array)    
    """
    try:
        b.resize(b.size, 1)
        C = append(A, b, 1)
    except MemoryError:
        print 'Mem�ria insuficiente!'
        return empty(A.shape)
        
    i = 0
    j = 0
    print "EG:",
    while i < C.shape[0] and j < C.shape[1]:
        maxi = i
        for k in range(i+1, C.shape[0]):
            if abs(C[k][j]) > abs(C[maxi,j]): maxi = k
        if C[maxi, j]:
            C[[i,maxi]] = C[[maxi,i]]
            C[i] /= C[i][j]
            for u in range(i+1, C.shape[0]): C[u] -= C[i]*C[u,j]
            i += 1
        j += 1
    print j
    x = array([])
    for l in range(C.shape[0]):
        x_ = C[-1-l][-1]
        for m in range(l): x_ -= C[-1-l][-2-m]*x[m]
        x = append(x, x_)
    return x[::-1]
    
def SSLD5D(D1, D2, D3, D4, D5, b):
    """Solu��o do Sistema Linear Diretamente de uma matrix Penta-Diagonal pelo
    m�todo de Elimina��o Gaussiana. Retorna um array com a solu��o.
    
    Sistema:    A*x = b,   A[n][n] e b[n**2]
    
    Argumentos:
        D1:    menor diagonal inferior      (array: n**2-n)
        D2:    diagonal inferior            (array: n**2-1)
        D3:    diagonal principal           (array: n**2)
        D4:    diagonal superior            (array: n**2-1)
        D5:    maior diagonal superior      (array: n**2-n)
        b:     vetor dos valores            (array: n**2)
        
        *n: n�mero de linhas e colunas da matriz espar�a
    """
    return SSLD(M5D2M(D1, D2, D3, D4, D5), b)

def SSLI5D_SOR(D1, D2, D3, D4, D5, b, w=0.5, x_=None, erro=1.0E-6):
    """Solu��o do Sistema Linear Indiretamente de uma matrix Penta-Diagonal pelo
    m�todo de SOR(Successive Over-Relaxation). Retorna um array com a solu��o.
    
    Sistema:  A*x = b,   A[n][n] e b[n**2]
    
    Argumentos:
        D1:    menor diagonal inferior      (array: n**2-n)
        D2:    diagonal inferior            (array: n**2-1)
        D3:    diagonal principal           (array: n**2)
        D4:    diagonal superior            (array: n**2-1)
        D5:    maior diagonal superior      (array: n**2-n)
        b:     vetor dos valores            (array: n**2) 
        w:     coeficiente do m�todo de SOR (velocidade de converg�ncia)
        x_:    valores iniciais para x      (array)
        erro:  erro m�ximo admitido
        
        *n: n�mero de linhas e colunas da matriz espar�a
    """
    error = True      # Valor inicial para erro entrar no "while"
    N = b.size        # N�mero de termos no sistema
    n = int(N**0.5)   # N�mero de linhas e colunas
    j = 0             # Vari�vel para itera��o
    
    # Se x_ n�o for definido, defina-o
    if x_ is None: x_ = array([0.0 for i in b])
    
    while(error > erro):
        # Na primeira linha, h� 3 elementos
        x = array([(1-w)*x_[0] + w*(b[0] - D4[0]*x_[1] - D5[0]*x_[n])/D3[0]])
        # Da segunda linha at� a n-�ssima linha, h� 4 elementos
        for i in range(1, n):
            x = append(x, (1-w)*x_[i] + w*(b[i]- D2[i-1]*x[i-1] - D4[i]*x_[i+1] - D5[i]*x_[i+n])/D3[i])
        # Da (n+1)-�ssima linha at� a (n**2-n)-�ssima linha, h� 5 elementos
        for i in range(n, N-n):
            x = append(x, (1-w)*x_[i] + w*(b[i] - D1[i-n]*x[i-n] - D2[i-1]*x[i-1] - D4[i]*x_[i+1] - D5[i]*x_[i+n])/D3[i])
        # Da(n**2-n+1)-�ssima linha at� (n**2-1)-�ssima linha, h� 4 elementos
        for i in range(N-n, N-1):
            x = append(x, (1-w)*x_[i] + w*(b[i] - D1[i-n]*x[i-n] - D2[i-1]*x[i-1] - D4[i]*x_[i+1])/D3[i])
        # Na �ltima linha, h� 3 elementos
        x = append(x, (1-w)*x_[N-1] + w*(b[N-1]- D1[N-n-1]*x[N-n-1] - D2[N-2]*x[N-2])/D3[N-1])
        # C�lculo do novo erro        
        error = abs(x-x_).max()
        x_ = x
        j+=1
    print 'SOR:', j
    return x
    
def SSLI5D_GJ(D1, D2, D3, D4, D5, b, x_=None, erro=1.0E-6):
    """Solu��o do Sistema Linear Indiretamente de uma matrix Penta-Diagonal pelo
    m�todo de Gauss-Jacobi. Retorna um array com a solu��o.
    
    Sistema:    A*x = b,   A[n][n] e b[n**2]
    
    Argumentos:
        D1:    menor diagonal inferior      (array: n**2-n)
        D2:    diagonal inferior            (array: n**2-1)
        D3:    diagonal principal           (array: n**2)
        D4:    diagonal superior            (array: n**2-1)
        D5:    maior diagonal superior      (array: n**2-n)
        b:     vetor dos valores            (array: n**2) 
        w:     coeficiente do m�todo de SOR (velocidade de converg�ncia)
        x_:    valores iniciais para x      (array)
        erro:  erro m�ximo admitido
        
        *n: n�mero de linhas e colunas da matriz espar�a
    """
    error = True      # Valor inicial para erro entrar no "while"
    N = b.size        # N�mero de termos no sistema
    n = int(N**0.5)   # N�mero de linhas e colunas
    j = 0             # Vari�vel para itera��o
    
    # Se x_ n�o for definido, defina-o
    if x_ is None: x_ = array([0.0 for i in b])
    
    while(error > erro):
        # Na primeira linha, h� 3 elementos
        x = array([(b[0] - D4[0]*x_[1] - D5[0]*x_[n])/D3[0]])
        # Da segunda linha at� a n-�ssima linha, h� 4 elementos
        for i in range(1, n):
            x = append(x, (b[i]- D2[i-1]*x_[i-1] - D4[i]*x_[i+1] - D5[i]*x_[i+n])/D3[i])
        # Da (n+1)-�ssima linha at� a (n**2-n)-�ssima linha, h� 5 elementos
        for i in range(n, N-n):
            x = append(x, (b[i] - D1[i-n]*x_[i-n] - D2[i-1]*x_[i-1] - D4[i]*x_[i+1] - D5[i]*x_[i+n])/D3[i])
        # Da(n**2-n+1)-�ssima linha at� (n**2-1)-�ssima linha, h� 4 elementos
        for i in range(N-n, N-1):
            x = append(x, (b[i] - D1[i-n]*x_[i-n] - D2[i-1]*x_[i-1] - D4[i]*x_[i+1])/D3[i])
        # Na �ltima linha, h� 3 elementos
        x = append(x, (b[N-1]- D1[N-n-1]*x_[N-n-1] - D2[N-2]*x_[N-2])/D3[N-1])
        # C�lculo do novo erro     
        error = abs(x-x_).max()
        x_ = x
        j+=1
    print 'GJ:', j
    return x
    
def SSLI5D_GS(D1, D2, D3, D4, D5, b, x_=None, erro=1.0E-6):
    """Solu��o do Sistema Linear Indiretamente de uma matrix Penta-Diagonal pelo
    m�todo de Gauss-Seidel. Retorna um array com a solu��o.
    
    Sistema:    A*x = b,   A[n][n] e b[n**2]
    
    Argumentos:
        D1:    menor diagonal inferior      (array: n**2-n)
        D2:    diagonal inferior            (array: n**2-1)
        D3:    diagonal principal           (array: n**2)
        D4:    diagonal superior            (array: n**2-1)
        D5:    maior diagonal superior      (array: n**2-n)
        b:     vetor dos valores            (array: n**2) 
        w:     coeficiente do m�todo de SOR (velocidade de converg�ncia)
        x_:    valores iniciais para x      (array)
        erro:  erro m�ximo admitido
        
        *n: n�mero de linhas e colunas da matriz espar�a
    """
    error = True      # Valor inicial para erro entrar no "while"
    N = b.size        # N�mero de termos no sistema
    n = int(N**0.5)   # N�mero de linhas e colunas
    j = 0             # Vari�vel para itera��o
    
    # Se x_ n�o for definido, defina-o
    if x_ is None: x_ = array([0.0 for i in b])
    
    while(error > erro):
        # Na primeira linha, h� 3 elementos
        x = array([(b[0] - D4[0]*x_[1] - D5[0]*x_[n])/D3[0]])
        # Da segunda linha at� a n-�ssima linha, h� 4 elementos
        for i in range(1, n):
            x = append(x, (b[i]- D2[i-1]*x[i-1] - D4[i]*x_[i+1] - D5[i]*x_[i+n])/D3[i])
        # Da (n+1)-�ssima linha at� a (n**2-n)-�ssima linha, h� 5 elementos
        for i in range(n, N-n):
            x = append(x, (b[i] - D1[i-n]*x[i-n] - D2[i-1]*x[i-1] - D4[i]*x_[i+1] - D5[i]*x_[i+n])/D3[i])
        # Da(n**2-n+1)-�ssima linha at� (n**2-1)-�ssima linha, h� 4 elementos
        for i in range(N-n, N-1):
            x = append(x, (b[i] - D1[i-n]*x[i-n] - D2[i-1]*x[i-1] - D4[i]*x_[i+1])/D3[i])
        # Na �ltima linha, h� 3 elementos
        x = append(x, (b[N-1]- D1[N-n-1]*x[N-n-1] - D2[N-2]*x[N-2])/D3[N-1])
        # C�lculo do novo erro     
        error = abs(x-x_).max()
        x_ = x
        j+=1
    print 'GS:', j
    return x