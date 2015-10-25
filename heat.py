#!/usr/bin/env python
# -*- coding: utf-8 -*-


'''
Este script resuelve un problema simple de diffusion en 1D.
La ecuación a resover es:
    dT/dt = d2T/dx2; T(0,x) = sin(pi * x); T(t, 0) = T(t, 1) = 0
'''

from __future__ import division
import numpy as np


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
    Aplus = -1  * r
    Acero = (1+2 * r)
    Aminus = -1 * r
    alpha[0] = 0
    beta[0] = 0 # viene de la condicion de borde T(t, 0) = 0
    for i in range(1, N_steps):
        alpha[i] = -Aplus / (Acero + Aminus*alpha[i-1])
        beta[i] = (b[i] - Aminus*beta[i-1]) / (Aminus*alpha[i-1] + Acero)

# Main

# setup
N_steps = 5
h = 1 / (N_steps - 1)
dt = h**2 / 2 # Este es el máximo teórico para el metodo explicito
r = dt / 2 / h**2

T = np.zeros(N_steps)
T_next = np.zeros(N_steps)
b = np.zeros(N_steps)
alpha = np.zeros(N_steps)
beta = np.zeros(N_steps)

inicializa_T(T, N_steps, h)

calcula_b(b, N_steps, r)
calcula_alpha_y_beta(alpha, beta, b, r, N_steps)

# Avanza T en el tiempo
T_next[0] = 0
T_next[-1] = 0
for i in range(N_steps - 2, 0, -1):
    T_next[i] = alpha[i] * T[i+1] + beta[i]
