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



# Main

# setup
N_steps = 5
h = 1 / (N_steps - 1)

T = np.zeros(N_steps)

inicializa_T(T, N_steps, h)
print T
