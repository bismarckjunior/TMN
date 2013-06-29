# -*- coding: iso-8859-1 -*-
"""
Autores: Bismarck Gomes Souza Júnior e Fernando Vizeu Santos
Versão: 1.0
Data: 28/07/2011

   Esse script gera três gráficos (soluções pelos métodos: SOR, Gauss-Jacobi e 
Gauss-Seidel) para a pressão do reservatório nxn. Salva as figuras e gera três 
arquivos de textos com os dados de pressão para cada bloco pelos três métodos.
"""
from time import time       # Função para calcular o tempo
import SistemaLinear as SL  # Funções que resolvem sistemas lineares 
from matplotlib import cm   # Função usada pra plotar colorido (3D)
import matplotlib.pyplot as plt          # Plotar graficos
from mpl_toolkits.mplot3d import axes3d  # Plotar gráficos com visualização 3D
from numpy import array, ones, zeros, pi, log, meshgrid, linspace


def Reservatorio(n=10, p0=None, inj=None, prod=None, erro=1E-5, fun = SL.SSLI5D_GS):
    """Função que calcula a pressão num reservatório de acordo com a equação 
    abaixo. Utiliza a condição de contorno tal que dp/dx=dp/dy=0 para as fron-
    teiras. Os pontos estão compreendidos entre 0 e 1 formando uma matriz nxn.
    Deve-se passar os dados dos poços produtores e injetores. Para representar 
    a matriz esparsa  dos coeficientes, armazena-se apenas os valores das diago-
    nais. Retorna uma matriz (array) dos valores de pressões em cada bloco.    
    
    Equação:
        d2p/dx2 + d2p/dy2 = - 2*pi/log(0.2*dx/r_w)*(p_bh - p),   p = p(x,y)
        
    Argumentos:
        n:       número de divisões entre 0 e 1 
        p0:      valores iniciais para a pressão (array)
        inj:     dados do poço injetor           [(i,j, valor da pressão)]
        prod:    dados do poço produtor          [(i,j, valor da pressão)]
        erro:    erro da iteração aceitável
        fun:     função para solucionar o sistema linear:
                   SL.SSLD: Eliminação gaussiana
                   SL.SSLI5D_GJ: Gauss-Jacobi
                   SL.SSLI5D_GS: Gauss-Seidel
                   SL.SSLI5D_SOR: Successive Over-Relaxation 
     
    Exemplo:
        Reservatorio(n=100, inj=(20,20,-1), prod=(70,70,1))
    """
    dx = 1./n         # Tamanho dos blocos
    r_e = 2**(-0.5)   # Raio externo (Blocos num quadrado 1x1)
    r_w = 0.2*dx      # Raio do poço (Peaceman)
    b  = zeros(n**2)  # Vetor de valores do Sistema Linear
    if inj is None: inj = (n/5, n/5, 1)          # Dados do poço injetor
    if prod is None: prod = (7*n/10, 7*n/10, -1) # Dados do poço produtor
    
    #Inserindo os dados dos poços
    b[inj[0]*n + inj[1]] = -inj[2]*2.*dx**2*pi/log(r_e/r_w)
    b[prod[0]*n + prod[1]] = -prod[2]*2.*dx**2*pi/log(r_e/r_w)
    
    #Diagonal principal
    D3 = ones(n**2)*( -4. - 2.*dx**2*pi/log(r_e/r_w) )
           
    D1 = ones(n**2-n)  # Blocos ao Sul
    D5 = ones(n**2-n)  # Blocos ao Norte
    D2 = ones(n**2-1)  # Blocos ao Oeste
    D4 = ones(n**2-1)  # Blocos ao Leste
    
    D4[0] += 1.
    D2[-1] += 1.
    
    for i in range(n):  
        D1[-(i+1)] += 1.      # Condição de contorno na parte inferior
        D5[i] += 1.           # Condição de contorno na parte superior
        
    for i in range(n, n**2, n): 
        D4[i] += 1.           # Condição de contorno na parte lateral esquerda
        D2[-i] -= 1.          # Condição de contorno na parte lateral esquerda        
        D4[i-1] -= 1.         # Condição de contorno na parte lateral direita
        D2[-(i+1)] += 1.      # Condição de contorno na parte lateral direita
    
    if fun == SL.SSLD5D:
        # O método direto não recebe valores iniciais para p, nem o erro
        p = array(fun(D1, D2, D3, D4, D5, b))
    else:
        p = array(fun(D1, D2, D3, D4, D5, b, x_=p0, erro=erro))
    p.resize(n,n)
    return p
    
if __name__ == "__main__":
    # Número de divisões entre 0 e 1
    n_ = 10
    
    # Criação da malha
    x, y = meshgrid(linspace(0,1,n_), linspace(0,1, n_))
    
    # Reservatório pelo método da eliminação gaussiana
    t0 = time()
    p_EG = Reservatorio(n=n_, fun=SL.SSLI5D_GJ)
    
    # Reservatório pelo método Gauss-Jacobi
    t1 = time()
    p_GJ = Reservatorio(n=n_, fun=SL.SSLI5D_GJ)#,inj=(6,6,1), prod=(21,21,-1),)
    
    # Reservatório pelo método Gauss-Seidel
    t2 = time()
    p_GS = Reservatorio(n=n_, fun=SL.SSLI5D_GS)#, inj=(6,6,1), prod=(21,21,-1))
    
    # Reservatório pelo método SOR
    t3 = time()    
    p_SOR = Reservatorio(n=n_, fun=SL.SSLI5D_SOR)#, inj=(6,6,1), prod=(21,21,-1)
    t4 = time()
    
    # Imprimindo os tempos de execução
    print '\nEliminacao Gaussiana:', t1-t0, 's'
    print 'Metodo Gauss-Jacobi: ', t2-t1, 's'
    print 'Metodo Gauss-Seidel: ', t3-t2, 's'
    print 'Metodo SOR:          ', t4-t3, 's'
    
    # Salvando em arquivo de texto
    for p, name in zip((p_EG, p_GJ, p_GS, p_SOR),('Eliminacao-Gaussiana.txt', 'Gauss-Jacobi.txt', 'Gauss-Seidel.txt', 'SOR.txt')):
        f = file(name, 'w')
        f.write('Metodo de '+name[:-4]+'\n\n')
        f.write(' i\j'+' '*8 + (' '*12).join([str(a) for a in range(n_)]))
        for i in range(p.shape[0]):
            f.write('\n%2i  ' % (i))
            for j in range(p.shape[1]):
                f.write('%13.7f' % (p[i][j]))
        f.close()
     
    # Plotando e salvando
    fig_3d = [plt.figure(i) for i in range(4)]
    for fig_, p, name in zip(fig_3d ,(p_EG, p_GJ, p_GS, p_SOR), ('Eliminacao-Gaussiana', 'Gauss-Jacobi', 'Gauss-Seidel', 'SOR')):
        ax = fig_.add_axes((0,0,1,1),projection='3d')
        N = (p-p.min())/(p.max()-p.min()) # Normalização da pressão (p)
        ax.plot_surface(x,y,p, rstride=1, cstride=1, alpha=0.6, facecolors=cm.jet(N))
        ax.set_title(name)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('p(x,y)')        
        fig_.savefig(name+'.png')
    plt.show()
    
    fig_2d = [plt.figure(i) for i in range(4)]
    for fig_, p, name in zip(fig_2d ,(p_EG, p_GJ, p_GS, p_SOR), ('Eliminacao-Gaussiana', 'Gauss-Jacobi', 'Gauss-Seidel', 'SOR')):
        ax = fig_.add_axes((0,0,1,1))
        cset = ax.contour(x,y, p, 20)
        ax.set_title(name)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.clabel(cset, inline=1, fotsize=10)
        fig_.savefig(name+'-2d.png')
        #ax.set_zlabel('p(x,y)')"""
    plt.show()