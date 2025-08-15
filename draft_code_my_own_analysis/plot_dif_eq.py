# -----------------------------------------------------------------------------
# Copyright (C) 2025 Salvador de la Torre Gonzalez
# Co-author: Luciana Melina Luque
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# -----------------------------------------------------------------------------

# # ------------------------------------------------------------------------------
# # Script to simulate and plot the time evolution of cell volume as done in 
# # Luciana's paper in order to study necrosis
# #
# # Solves ODEs for fluid and solid volume changes over time:
# #   - dF/dt = k_f * (F_target - F)
# #   - dNs/dt = k_n * (Ns_target - Ns)
# #   - dCs/dt = k_c * (Cs_target - Cs)
# #
# # Total and compartment volumes are updated accordingly.
# # ------------------------------------------------------------------------------



# import numpy as np
# import matplotlib.pyplot as plt

# # === Initial Conditions ===
# # Based on your Volumen constructor
# total_initial = 2494.0
# total=total_initial
# fraccion_de_fluido = 0.75
# fluido = fraccion_de_fluido * total
# solido = total - fluido

# nuclear = 540.0
# nuclear_fluido = fraccion_de_fluido * nuclear
# nuclear_solido = nuclear - nuclear_fluido

# citoplasmatico = total - nuclear
# citoplasmatico_fluido = fraccion_de_fluido * citoplasmatico
# citoplasmatico_solido = citoplasmatico - citoplasmatico_fluido

# # Parameters
# citoplasma_tasa_de_cambio = 0.0032 /60.0
# nucleo_tasa_de_cambio = 0.013 / 60.0
# fluido_tasa_de_cambio = 0.050 /60.0
# tasa_de_calcificacion = 0.0042 / 60.0

# relacion_citoplasma_nucleo = citoplasmatico / ( 1e-16 + nuclear)

# #TARGET
# target_nucleo_solido = 0
# target_relacion_citoplasma_nucleo = 0.0
# target_citoplasma_solido = 0
# target_fraccion_fluido = 1.0

# # For plotting
# time = []
# total_volume = []
# fluid_volume = []
# solid_volume = []
# fluid_fraction = []

# # Simulation settings
# dt = 6 # min
# t_max = 50000  # total time in minutes
# steps = int(t_max / dt)

# # Initialize dynamic variables
# fraccion_calcificada = 0.0

# t_real_t_1= None

# for step in range(steps):
#     t = step * dt

#     # Update fluid volume
#     target_fluido = target_fraccion_fluido * total
#     fluido += dt * fluido_tasa_de_cambio * (target_fluido - fluido)
#     if fluido < 0.0:
#         fluido = 0.0

#     # Update fluid distribution
#     nuclear_fluido = (nuclear / total) * fluido
#     citoplasmatico_fluido = fluido - nuclear_fluido

#     # Update solid volumes
#     nuclear_solido += dt * nucleo_tasa_de_cambio * (target_nucleo_solido - nuclear_solido)
#     nuclear_solido = max(0.0, nuclear_solido)

#     target_citoplasma_solido = target_relacion_citoplasma_nucleo * target_nucleo_solido
#     citoplasmatico_solido += dt * citoplasma_tasa_de_cambio * (target_citoplasma_solido - citoplasmatico_solido)
#     citoplasmatico_solido = max(0.0, citoplasmatico_solido)

#     # Update total components
#     solido = nuclear_solido + citoplasmatico_solido
#     nuclear = nuclear_solido + nuclear_fluido
#     citoplasmatico = citoplasmatico_solido + citoplasmatico_fluido
#     total = nuclear + citoplasmatico

#     fraccion_calcificada += dt * tasa_de_calcificacion * (1 - fraccion_calcificada)


#     fraccion_de_fluido = fluido / (1e-16 + total)

#     # Store for plotting
#     time.append(t)
#     total_volume.append(total)
#     fluid_volume.append(fluido)
#     solid_volume.append(solido)
#     fluid_fraction.append(fraccion_de_fluido)
#     if total > 2*total_initial and t_real_t_1==None:
#         print(t, "min")
#         t_real_t_1= t
# # === Plotting ===
# # plt.figure(figsize=(10, 6))

# # plt.plot(time, total_volume, label='Total Volume')
# # plt.plot(time, fluid_volume, label='Fluid Volume')
# # plt.plot(time, solid_volume, label='Solid Volume')
# # plt.plot(time, fluid_fraction, label='Fluid Fraction')

# # if t_real_t_1 is not None:
# #     plt.axvline(x=t_real_t_1, color='red', linestyle='--', label=f'real_t total Volume ({t_real_t_1} min, step={t_real_t_1//dt})')


# # plt.xlabel("Time (minutes)")
# # plt.ylabel("Volume (μm³) or Fraction")
# # plt.title("Cell Volume Dynamics Over Time")
# # plt.legend()
# # plt.grid(True)
# # plt.tight_layout()
# # plt.show()

# # #_-----------------------------------this is the equation I am using in my symplification to have just one volume type-------------------------------------------
# import matplotlib.pyplot as plt

# # Parameters
# initial_volume = 2494.0
# current_volume = initial_volume
# target_volume = 4*initial_volume         
# relaxation_rate = 0.000059#0.054*(0.35/60)+0.75*(0.05/60)+0.196*(1.0/60)          
# dt = 6                           
# t_max = 50000                      
# steps = int(t_max / dt)

# # Lists for plotting
# time2 = []
# volumes = []

# # To track when volume real_ts
# t_real_t_2= None

# # Simulation loop
# for step in range(steps):
#     t = step * dt

#     # Update the volume using the formula
#     current_volume += (target_volume - current_volume) * relaxation_rate * dt
#     volumes.append(current_volume)
#     time2.append(t)

#     # Check if volume has real_td
#     if current_volume >= 2 * initial_volume and t_real_t_2 is None:
#         t_real_t_2 = t
#         print(f"Volume real_td at {t_real_t_2} minutes")

# # # Plotting
# # plt.figure(figsize=(10, 6))
# # plt.plot(time, volumes, label='Volume')

# # if t_real_t_2 is not None:
# #     plt.axvline(x=t_real_t, color='red', linestyle='--', label=f'real_t Volume ({t_real_t} min, step={t_real_t//dt})')

# # plt.xlabel("Time (minutes)")
# # plt.ylabel("Volume")
# # plt.title("Volume Evolution Over Time")
# # plt.grid(True)
# # plt.legend()
# # plt.tight_layout()
# # plt.show()

# # === Plotting Side by Side ===
# fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), sharex=True, sharey=True)

# # Plot first model
# ax1.plot(time, total_volume, label='Total Volume')
# ax1.plot(time, fluid_volume, label='Fluid Volume')
# ax1.plot(time, solid_volume, label='Solid Volume')
# ax1.plot(time, fluid_fraction, label='Fluid Fraction')
# if t_real_t_1:
#     ax1.axvline(t_real_t_1, color='red', linestyle='--', label=f'real_t at {t_real_t_1} min, step={t_real_t_1//dt}')
# ax1.set_title("Model 1: Detailed Volume Dynamics")
# ax1.set_xlabel("Time (minutes)")
# ax1.set_ylabel("Volume (μm³) or Fraction")
# ax1.grid(True)
# ax1.legend()

# # Plot second model
# ax2.plot(time2, volumes, label='Volume', color='purple')
# if t_real_t_2:
#     ax2.axvline(t_real_t_2, color='red', linestyle='--', label=f'real_t at {t_real_t_2} min, step={t_real_t_2//dt}')
# ax2.set_title("Model 2: Simplified Growth Model")
# ax2.set_xlabel("Time (minutes)")
# ax2.grid(True)
# ax2.legend()

# plt.tight_layout()
# plt.show()


# ------------------------------------------------------------------------------
# Script to simulate and plot the time evolution of cell volume as done in 
# Luciana's paper in order to study necrosis
#
# Solves ODEs for fluid and solid volume changes over time:
#   - dF/dt = k_f * (F_target - F)
#   - dNs/dt = k_n * (Ns_target - Ns)
#   - dCs/dt = k_c * (Cs_target - Cs)
#
# Total and compartment volumes are updated accordingly.
# ------------------------------------------------------------------------------



import numpy as np
import matplotlib.pyplot as plt

# === Initial Conditions ===
# Based on your Volumen constructor
total_initial = 2494.0
total=total_initial
fraccion_de_fluido = 0.75
fluido = fraccion_de_fluido * total
solido = total - fluido

nuclear = 540.0
nuclear_fluido = fraccion_de_fluido * nuclear
nuclear_solido = nuclear - nuclear_fluido

citoplasmatico = total - nuclear
citoplasmatico_fluido = fraccion_de_fluido * citoplasmatico
citoplasmatico_solido = citoplasmatico - citoplasmatico_fluido

# Parameters
citoplasma_tasa_de_cambio = 0.27/60.0
nucleo_tasa_de_cambio = 0.33/60.0
fluido_tasa_de_cambio = 3.0 / 60.0
    
tasa_de_calcificacion = 0.0

relacion_citoplasma_nucleo = citoplasmatico / ( 1e-16 + nuclear)

#TARGET
target_nucleo_solido = 2*nuclear_solido
target_relacion_citoplasma_nucleo = relacion_citoplasma_nucleo
target_citoplasma_solido = 2*citoplasmatico_solido
target_fraccion_fluido = fraccion_de_fluido

# For plotting
time = []
total_volume = []
fluid_volume = []
solid_volume = []
fluid_fraction = []

# Simulation settings
dt = 6 # min
t_max = 500  # total time in minutes
steps = int(t_max / dt)

# Initialize dynamic variables
fraccion_calcificada = 0.0

t_real_t_1= None

for step in range(steps):
    t = step * dt

    # Update fluid volume
    target_fluido = target_fraccion_fluido * total
    fluido += dt * fluido_tasa_de_cambio * (target_fluido - fluido)
    if fluido < 0.0:
        fluido = 0.0

    # Update fluid distribution
    nuclear_fluido = (nuclear / total) * fluido
    citoplasmatico_fluido = fluido - nuclear_fluido

    # Update solid volumes
    nuclear_solido += dt * nucleo_tasa_de_cambio * (target_nucleo_solido - nuclear_solido)
    nuclear_solido = max(0.0, nuclear_solido)

    target_citoplasma_solido = target_relacion_citoplasma_nucleo * target_nucleo_solido
    citoplasmatico_solido += dt * citoplasma_tasa_de_cambio * (target_citoplasma_solido - citoplasmatico_solido)
    citoplasmatico_solido = max(0.0, citoplasmatico_solido)

    # Update total components
    solido = nuclear_solido + citoplasmatico_solido
    nuclear = nuclear_solido + nuclear_fluido
    citoplasmatico = citoplasmatico_solido + citoplasmatico_fluido
    total = nuclear + citoplasmatico

    fraccion_calcificada += dt * tasa_de_calcificacion * (1 - fraccion_calcificada)


    fraccion_de_fluido = fluido / (1e-16 + total)

    # Store for plotting
    time.append(t)
    total_volume.append(total)
    fluid_volume.append(fluido)
    solid_volume.append(solido)
    fluid_fraction.append(fraccion_de_fluido)
    if total > 2*total_initial and t_real_t_1==None:
        print(t, "min")
        t_real_t_1= t
# === Plotting ===
# plt.figure(figsize=(10, 6))

# plt.plot(time, total_volume, label='Total Volume')
# plt.plot(time, fluid_volume, label='Fluid Volume')
# plt.plot(time, solid_volume, label='Solid Volume')
# plt.plot(time, fluid_fraction, label='Fluid Fraction')

# if t_real_t_1 is not None:
#     plt.axvline(x=t_real_t_1, color='red', linestyle='--', label=f'real_t total Volume ({t_real_t_1} min, step={t_real_t_1//dt})')


# plt.xlabel("Time (minutes)")
# plt.ylabel("Volume (μm³) or Fraction")
# plt.title("Cell Volume Dynamics Over Time")
# plt.legend()
# plt.grid(True)
# plt.tight_layout()
# plt.show()

# #_-----------------------------------this is the equation I am using in my symplification to have just one volume type-------------------------------------------
import matplotlib.pyplot as plt

# Parameters
initial_volume = 2494.0
current_volume = initial_volume
target_volume = 2*initial_volume         
relaxation_rate = 0.054*(0.33/60)+0.75*(3.0/60)+0.196*(0.27/60)        
dt = 6                           
t_max = 500                      
steps = int(t_max / dt)

# Lists for plotting
time2 = []
volumes = []

# To track when volume real_ts
t_real_t_2= None

# Simulation loop
for step in range(steps):
    t = step * dt

    # Update the volume using the formula
    current_volume += (target_volume - current_volume) * relaxation_rate * dt
    volumes.append(current_volume)
    time2.append(t)

    # Check if volume has real_td
    if current_volume >= 2 * initial_volume and t_real_t_2 is None:
        t_real_t_2 = t
        print(f"Volume real_td at {t_real_t_2} minutes")

# # Plotting
# plt.figure(figsize=(10, 6))
# plt.plot(time, volumes, label='Volume')

# if t_real_t_2 is not None:
#     plt.axvline(x=t_real_t, color='red', linestyle='--', label=f'real_t Volume ({t_real_t} min, step={t_real_t//dt})')

# plt.xlabel("Time (minutes)")
# plt.ylabel("Volume")
# plt.title("Volume Evolution Over Time")
# plt.grid(True)
# plt.legend()
# plt.tight_layout()
# plt.show()

# === Plotting Side by Side ===
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), sharex=True, sharey=True)

# Plot first model
ax1.plot(time, total_volume, label='Total Volume')
ax1.plot(time, fluid_volume, label='Fluid Volume')
ax1.plot(time, solid_volume, label='Solid Volume')
ax1.plot(time, fluid_fraction, label='Fluid Fraction')
if t_real_t_1:
    ax1.axvline(t_real_t_1, color='red', linestyle='--', label=f'real_t at {t_real_t_1} min, step={t_real_t_1//dt}')
ax1.set_title("Model 1: Detailed Volume Dynamics")
ax1.set_xlabel("Time (minutes)")
ax1.set_ylabel("Volume (μm³) or Fraction")
ax1.grid(True)
ax1.legend()

# Plot second model
ax2.plot(time2, volumes, label='Volume', color='purple')
if t_real_t_2:
    ax2.axvline(t_real_t_2, color='red', linestyle='--', label=f'real_t at {t_real_t_2} min, step={t_real_t_2//dt}')
ax2.set_title("Model 2: Simplified Growth Model")
ax2.set_xlabel("Time (minutes)")
ax2.grid(True)
ax2.legend()

plt.tight_layout()
plt.show()