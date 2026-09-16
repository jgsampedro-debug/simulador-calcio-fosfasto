import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# --- CONFIGURACIÓN DE LA PÁGINA ---
st.set_page_title("Simulador de Equilibrio de Fosfatos", layout="wide")
st.title("🧪 Equilibrio Químico y Precipitación de Fosfato Dicálcico")
st.markdown("Ajusta los parámetros en la barra lateral para ver cómo cambia la concentración de Calcio libre en función del pH.")

# --- BARRA LATERAL (CONTROLES) ---
st.sidebar.header("Parámetros de Entrada")

# Sliders interactivos
P_total_inicial = st.sidebar.slider("Concentración total de Fósforo (P_total) [M]", 0.001, 0.1, 0.01, step=0.001, format="%.3f")
Kps_exponente = st.sidebar.slider("Exponente del Kps (10^-x)", 5.0, 10.0, 6.6, step=0.1) # 2.5e-7 es aprox 10^-6.6
Kps_dcpd = 10**(-Kps_exponente)

st.sidebar.markdown("### Rango de pH")
pH_min, pH_max = st.sidebar.slider("Selecciona el rango de pH", 2.0, 12.0, (4.0, 10.0))

# --- CONSTANTES ---
pKa1, pKa2, pKa3 = 2.15, 7.20, 12.35
K1 = 10**(-pKa1)
K2 = 10**(-pKa2)
K3 = 10**(-pKa3)

# --- FUNCIONES DE CÁLCULO ---
def calcular_especies_fosfato(pH):
    H = 10**(-pH)
    alfa0 = H**3 / (H**3 + H**2*K1 + H*K1*K2 + K1*K2*K3)
    alfa1 = (H**2*K1) / (H**3 + H**2*K1 + H*K1*K2 + K1*K2*K3)
    alfa2 = (H*K1*K2) / (H**3 + H**2*K1 + H*K1*K2 + K1*K2*K3)
    alfa3 = (K1*K2*K3) / (H**3 + H**2*K1 + H*K1*K2 + K1*K2*K3)
    return alfa0, alfa1, alfa2, alfa3

def sistema_equilibrio(vars, pH, P_total):
    Ca2_libre, HPO4_libre = vars
    res1 = Ca2_libre * HPO4_libre - Kps_dcpd
    PO4_total_calc = HPO4_libre * (1 + 10**(-pH)/K2) 
    res2 = PO4_total_calc - P_total
    return [res1, res2]

# --- PROCESAMIENTO ---
pH_valores = np.linspace(pH_min, pH_max, 100)
calcio_libre = []

for pH in pH_valores:
    guess = [1e-3, 1e-4]
    solucion = fsolve(sistema_equilibrio, guess, args=(pH, P_total_inicial))
    calcio_libre.append(solucion[0])

# --- RENDERIZADO DE LA INTERFAZ ---
col1, col2 = st.columns([2, 1])

with col1:
    st.subheader("Gráfico de Concentración de Calcio")
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(pH_valores, calcio_libre, label='Calcio libre disuelto (M)', color='#1f77b4', linewidth=2.5)
    ax.set_xlabel('pH', fontsize=12)
    ax.set_ylabel('Concentración de Ca2+ (mol/L)', fontsize=12)
    ax.set_yscale('log')
    ax.grid(True, which="both", ls="--", alpha=0.6)
    ax.legend(fontsize=10)
    st.pyplot(fig)

with col2:
    st.subheader("Valores Actuales")
    st.metric(label="Kps Calculado", value=f"{Kps_dcpd:.2e}")
    st.metric(label="P Total Inicial", value=f"{P_total_inicial} M")
    st.metric(label="Rango de pH", value=f"{pH_min} - {pH_max}")
