# WingBox Structural Analysis

Este projeto realiza a simulação estrutural de uma **caixa de asa (wingbox)** de uma aeronave sujeita a cargas aerodinâmicas. O código utiliza programação orientada a objetos em Python para organizar os módulos de cálculo, entrada de dados e visualização dos resultados.

---

## 📁 Estrutura dos Arquivos

- `structural.py`  
  Contém a **estrutura principal de classes** para cálculo dos esforços na asa. As principais classes são:
  - `structure`: define os dados da aeronave e da seção estrutural.
  - `aero_struct`: herda de `structure` e implementa a lógica de cargas aerodinâmicas, momentos fletores, esforços cortantes, deflexões e tensões normais.
  - `torsion`: herda de `aero_struct` e implementa o cálculo do **torque**, **ângulo de torção** e **fluxo de cisalhamento interno**.

- `visualization.py`  
  Contém funções auxiliares para geração de **plots**:
  - Distribuições de carga e deflexão.
  - Esforço cortante e momento fletor.
  - Campo de tensões σₓ(y,z).
  - Geometria da seção transversal com centroide e eixo neutro.

- `sim_setup.py`  
  Script principal que:
  - Define os **dados de entrada** da aeronave e seção estrutural.
  - Inicializa objetos de cálculo (`aero_struct`, `torsion`).
  - Calcula esforços e deflexões em quatro pontos críticos do diagrama V-n.
  - Gera os gráficos e resultados utilizando as funções de `visualization.py`.

---

## 📌 Funcionalidades

- 📐 Cálculo da geometria da seção H (multi-células).
- 🧮 Cálculo de:
  - Cargas distribuídas (sustentação e arrasto).
  - Esforços cortantes e momentos fletores.
  - Deflexões transversais.
  - Tensões normais σₓ(y,z).
  - Torque interno e ângulo de torção.
- 📊 Visualizações gráficas com `matplotlib`, incluindo:
  - Diagramas de força e momento.
  - Campo de tensões σₓ.
  - Ângulo de torção.
  - Visualização 3D de vetores de carga.

---

## ▶️ Como usar

1. Certifique-se de ter o Python 3 instalado com as bibliotecas:
   ```bash
   pip install numpy matplotlib
