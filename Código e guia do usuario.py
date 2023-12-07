from vpython import *
from random import *
import numpy as np

winw = 800
winh = winw

# Configuracoes da camera do Vpython
scene = canvas(title="Ciclo de Carnot", width=winw,
               height=winh, x=0, y=0, align='left')
scene.center = vector(0, 1.6, 0)
scene.forward = vector(0, -.2, -1)
scene.range = 3.5
scene.background = vector(0.788,0.788,0.788)

# Dimensoes da caixa:
# lado da base
l = 1
# altura
h = 2
# espessura
es = 0.1

# Criacao da biela
corC = color.blue

rCilindro = 0.8
fundo = cylinder(pos=vector(0, -es/2, 0),
                 axis=vector(0, es, 0), radius=rCilindro, color=corC, opacity=0.25)
circpath = paths.arc(angle1=0, angle2=2*pi, radius=rCilindro + 0.055)
rect = shapes.rectangle(width=0.1, height=2)
corpoC = extrusion(pos=vector(0, h/2, (rCilindro+0.055)/100),
                   path=circpath, shape=rect, color=color.blue, opacity=0.25)

# Criacao do pistao:
pistao = cylinder(pos=vector(0, 0.2*h + es/2, 0),
                  axis=vector(0, es, 0), radius=0.8)
pinP = cylinder(pos=pistao.pos + vector(0, es/2, -es/2),
                axis=vector(0, 0, es), radius=es/2, color=corC)


# Criacao da roda
rRoda = 0.8*l
roda = cylinder(pos=vector(0, 0, 0), axis=vector(
    0, 0, es), radius=rRoda, color=vector(0.5,0.5,0.5))
eixoR = cylinder(pos=vector(0, 0, 0), axis=vector(0,
                 0, 3*es), radius=es/2, color=color.red)

# Criacao do pino de conexao da vareta com a roda
Rp = 0.8*rRoda  # O pino de conexao esta a uma distancia de 0.8R do centro da roda
theta = -0.9*pi/2  # Angulo inicial desse pino
pino = cylinder(pos=Rp*vector(cos(theta), sin(theta), 0),
                axis=vector(0, 0, 3*es), radius=es/2, color=color.red)

conj = compound([roda, eixoR, pino])
conj.pos = vector(0, 2*h, 0)

# Salvando a localizacao desse pino, sera atualizada mais a frente
pinLoc = vector(conj.pos.x + Rp*cos(theta),
                conj.pos.y + Rp*sin(theta), pinP.pos.z)

# Criando a haste de ligacao:
haste = cylinder(pos=vector(pinP.pos.x, pinP.pos.y, 0),
                 axis=pinLoc - pinP.pos, radius=es/5, color=color.white)
Lh = mag(haste.axis)  # Comprimento da haste

# Criacao de um termometro de mercurio:

Tq = 400  # Temperatura do reservatorio quente
Tf = 339  # Temperatura do reservatorio frio
T = T0 = Tq

e = 1 - (Tf / Tq) # Rendimento
R = 8.31

merc = cylinder(pos=vector(fundo.pos.x + 1 + es/2, fundo.pos.y +
                es/2, l/2), radius=es/5, axis=vector(0, 0.7*h, 0), color=color.red)
sphere(pos=merc.pos, radius=es/2, color=color.red)
tubo = cylinder(pos=merc.pos, radius=es/3, axis=vector(0,
                0.7*h, 0), color=color.white, opacity=0.25)
bulbo = sphere(pos=merc.pos, radius=es, color=color.white, opacity=0.25)

Trot = label(pos=vector(fundo.pos.x+1/2+es/2+es,fundo.pos.y+2/2,1/2), xoffset=60, 
                line=0, box=10, opacity=0, color=color.black, text=str(int(round(T)))+' K')

tempQ = 0.7*h
tempF = 0.2*h


def Escala():  # Escala de temperatura arbitraria (so serve para variar a altura do mercurio no termometro)
    merc.axis.y = T*0.7*h/Tq

# Inicializacao do grafico

Wg = graph(xtitle='Ciclo', ytitle='W', xmin=0, xmax=10, ymin=0, ymax=100,
      x=2*winw, width=800, height=scene.height, align='left')
W = gcurve(graph = Wg, color=color.red, dot=True, dot_color=color.black)

PVg = graph(xtitle='V', ytitle='P', xmin=0.2, xmax=1.8, ymin=0, ymax=90,
      x=winw, width=600, height=scene.height, align='left')
PV = gcurve(graph = PVg, color=color.red, dot=True, dot_color=color.black)


P = P0 = 80  # Pressao inicial

Rpar = 0.03  # Raio de uma particula
dt = 0.01  # Passo

desvio = 1.1*Rp

Natomos = 50  # Numero de atomos
Matom = 4E-3/6E23  # Massa de um atomo
Ratom = 0.03  # Raio atomico
k = 1.4E-23  # Constante de Boltzmann

Atomos = []  # Lista com as particulas
cores = [color.red, color.yellow]
plist = []  # Lista com os momentos lineares
mlist = []  # Lista com as massas
desvio = 1.1*Rpar  # Margem de seguranca da posicao das particulas

# Sabendo que a enrgia cinetica media é p**2/(2*M) = (3/2)kT
pavg = sqrt(2*Matom*1.5*k*Tq)*(5E-5/dt)

for i in range(Natomos):
    # Definindo limmites de posição das particulas
    Lmin = -1/2+desvio
    Lmax = 1/2-desvio
    x = Lmin+(Lmax-Lmin)*random()
    Lmin = fundo.pos.y+es/2+desvio
    Lmax = pistao.pos.y-es/2-desvio
    y = Lmin+(Lmax-Lmin)*random()
    Lmin = -1/2+desvio
    Lmax = 1/2-desvio
    z = Lmin+(Lmax-Lmin)*random()
    Atomos.append(sphere(pos=vector(x, y, z),
                         radius=Ratom, color=cores[i % 2]))
    angle = pi*random()
    phi = 2*pi*random()
    px = pavg*sin(angle)*cos(phi)
    py = pavg*sin(angle)*sin(phi)
    pz = pavg*cos(angle)
    plist.append(vector(px, py, pz))
    mlist.append(Matom)

# Momento angular
L = vector(0, 0, 0)
# Expoente de Poisson gas monoatomico
gamma = 1.4
# Momento de Inercia
I = 10

# Determina qual parte do processo do ciclo esta sendo executado
processo = 0  # Expansao isotermica
yinicial = pistao.pos.y - es/2

processoRot = label(pos=vector(1.5, 3, 1), text = "Fase " + str(processo) + ": Fase Inicial", line=0, box=10, opacity=0, color=color.black)

u = 0.5
a = 0

while True:
    rate(150)

    # Rotacao da roda e atualizacao da posicao dos pinos e da haste
    # Forca aproximada para a pressao, atuando na direcao da haste
    torque = cross(pinLoc - conj.pos, P*norm(haste.axis))
    L += torque*dt  # Atualizacao do momento angular
    omega = L.z/I
    dtheta = omega*dt
    if omega >= 4 and T <= 340:
        break
    theta += dtheta

    conj.rotate(angle=dtheta, axis=vector(0, 0, 1))
    pinLoc = vector(conj.pos.x + Rp*cos(theta),
                    conj.pos.y + Rp*sin(theta), pinP.pos.z)

    # Diferenca de altura entre o pino e o pistao (o quanto a roda girou pra cima)
    dy = sqrt(Lh**2 - pinLoc.x**2)
    # "Volume" inicial
    y0 = pistao.pos.y - es/2
    # Atualizacao da altura do pistao
    pinP.pos.y = pinLoc.y - dy
    pistao.pos.y = pinP.pos.y - es/2
    # "Volume" final
    y = pistao.pos.y - es/2
    # Razao inversa dos volumes
    Vr = y0/y
    # Atualizacao da posicao da haste
    haste.pos = pinP.pos
    haste.axis = pinLoc - haste.pos

    if processo == 1 or processo == 3:  # Isotermas
        # Pressao e volume inversamente proporcionais
        P = P*Vr
    else:  # Processo adibatico
        P = P*Vr**gamma
    # Equacao geral dos gases
    T = T0*(P*y)/(P0*yinicial)
    
    Trot.text = str(int(round(T)))+' K'

    PV.plot(pos=(y, P))

    if processo == 0:  # Inicio da simulacao
        y0 = y = yinicial = pistao.pos.y - es/2
        PV.plot(pos=(y, P))
        plotcolor = color.red
        processo = 1
        processoRot.text = "Fase " + str(0) + ": Fase Inicial"

    if processo == 1:  # Expansao isotermica
        if y >= 1.1:
            processo = 2
            PV = gcurve(color=color.blue, dot=True, dot_color=color.black)
            PV.plot(pos=(y, P))
        processoRot.text = "Fase " + str(1) + ": Expansão Isotérmica"

    elif processo == 2:  # Expansao adiabatica
        yb = pistao.pos.y - es/2
        if merc.axis.y > tempF:
            merc.axis.y -= 0.03
        if pinLoc.x <= 0:
            processo = 3
            PV = gcurve(color=color.red, dot=True, dot_color=color.black)
            PV.plot(pos=(y, P))
        processoRot.text = "Fase " + str(2) + ": Expansão Adiabática"

    elif processo == 3:  # Compressao isotermica
        # OBS.: Onde aparece y, pode ser substituido por V, pois a área é constante
        if y < yinicial*(T0/T)**(1/(gamma - 1)):
            processo = 4
            PV = gcurve(color=color.blue, dot=True, dot_color=color.black)
            PV.plot(pos=(y, P))
        processoRot.text = "Fase " + str(3) + ": Compressão Isotérmica"
            
    else:  # Compressao adibaticaa
        if merc.axis.y < tempQ:
            merc.axis.y += 0.03
        # Voltando para a expansao isotermica
        if T >= Tq:
            processo = 1
            PV = gcurve(color=color.red, dot=True, dot_color=color.black)
            PV.plot(pos=(y, P))
        processoRot.text = "Fase " + str(4) + ": Compressão Adiabática"


    # Colisao com as paredes:
    for i in range(Natomos):
        p = plist[i]
        # Posicao futura da particula
        pos = Atomos[i].pos + (p/mlist[i])*dt
        colisao = False
        # Colisao em x:
        if not (-1/2 < pos.x < 1/2):
            plist[i].x *= -1
            colisao = True
        # Colisao em y:
        if not (fundo.pos.y+es/2+desvio < pos.y < pistao.pos.y-es/2-desvio):
            plist[i].y *= -1
            # Garantia de que nenhum atomo vai escapar do pistao
            if pos.y >= pistao.pos.y-es/2-desvio:
                Atomos[i].pos.y = pistao.pos.y-es/2-desvio
            colisao = True
        # Colisao em z:
        if not (-1/2 < pos.z < 1/2):
            plist[i].z *= -1
            colisao = True
        # Se nao há colisao, atualiza a posicao da particula
        if not colisao:
            Atomos[i].pos = pos
