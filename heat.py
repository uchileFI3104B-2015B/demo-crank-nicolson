#!/usr/bin/env python
# -*- coding: utf-8 -*-


'''
Este script resuelve un problema simple de diffusion en 1D.
La ecuación a resover es:
    dT/dt = d2T/dx2; T(0,x) = sin(pi * x); T(t, 0) = T(t, 1) = 0
'''

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


def inicializa_T(T, N_steps, h):
    '''
    Rellena T con las condiciones iniciales del problema.
    Se asegura que las condiciones en los bordes sean cero.
    '''
    for i in range(N_steps):
        x = i * h
        T[i] = np.sin(np.pi * x)
    T[0] = 0
    T[-1] = 0


def calcula_b(b, N_steps, r):
    for j in range(1, N_steps - 1):
        b[j] = r * T[j+1] + (1-2*r) * T[j] + r * T[j-1]


def calcula_alpha_y_beta(alhpa, beta, b, r, N_Steps):
    Aplus = -1 * r
    Acero = (1 + 2 * r)
    Aminus = -1 * r
    alpha[0] = 0
    beta[0] = 0  # viene de la condicion de borde T(t, 0) = 0
    for i in range(1, N_steps):
        alpha[i] = -Aplus / (Acero + Aminus*alpha[i-1])
        beta[i] = (b[i] - Aminus*beta[i-1]) / (Aminus*alpha[i-1] + Acero)


def avanza_paso_temporal(T, T_next, alpha, beta, N_steps):
    T_next[0] = 0
    T_next[-1] = 0
    for i in range(N_steps - 2, 0, -1):
        T_next[i] = alpha[i] * T_next[i+1] + beta[i]
    # T = T_next.copy() # Esto no funciona, hacerlo fuera de la funcion

# Main

# setup
N_steps = 41
N_pasos_temporales = 100

h = 1 / (N_steps - 1)
# dt = h**2 / 2 # Este es el máximo teórico para el metodo explicito
dt = 0.01
r = dt / 2 / h**2

T = np.zeros(N_steps)
T_next = np.zeros(N_steps)

b = np.zeros(N_steps)
alpha = np.zeros(N_steps)
beta = np.zeros(N_steps)

inicializa_T(T, N_steps, h)

# Queremos guardar las soluciones en cada paso
T_solucion = np.zeros((N_pasos_temporales, N_steps))
T_solucion[0, :] = T.copy()

for i in range(1, N_pasos_temporales):
    calcula_b(b, N_steps, r)
    calcula_alpha_y_beta(alpha, beta, b, r, N_steps)
    avanza_paso_temporal(T, T_next, alpha, beta, N_steps)
    T = T_next.copy()
    T_solucion[i, :] = T.copy()


# Plots

# ejemplo 1

x = np.linspace(0, 1, N_steps)

fig = plt.figure(1)
fig.clf()
ax = fig.add_subplot(111)

for i in range(0, N_pasos_temporales, 10):
    ax.plot(x, T_solucion[i, :])
ax.set_ylim(0, 1)

# ejemplo 2
# usar el plano x, t y plotear T en la 3a dimension
fig2 = plt.figure(2)
fig2.clf()
ax2 = fig2.add_subplot(111)
y = np.arange(0, N_pasos_temporales) * dt
X, Y = np.meshgrid(x, y)
ax2.pcolormesh(X, Y, T_solucion)

plt.show()
plt.draw()
