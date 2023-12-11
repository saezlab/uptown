from scipy.optimize import curve_fit
import pandas as pandas
import numpy as np
import matplotlib.pyplot as plt

# Read in the data
data = pandas.read_table('phosphorylated_prots_clean.tsv', sep='\t')
data.head()



# Parameter fitting
genes = data['Gene'].unique()

def logistic_function(x, A, B, C, D):
    return A + (C-A)/(1 + np.exp(-(x-B)/D))

for gene in genes:
    subset_data = data[data['Gene'] == gene]
    concentrations = subset_data['conc']
    ratios = subset_data['ratio']

    # Initial parameter estimation
    A = np.min(ratios)
    B = np.median(concentrations)
    C = np.max(ratios)
    D = 1
    initial_parameters = [A, B, C, D]

    optimal_params, cov = curve_fit(logistic_function, concentrations, ratios, p0=initial_parameters)



import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Definir la función logística de tres parámetros
def logistic_function(x, A, B, C, D):
    return A + (C - A) / (1 + np.exp(-(x - B) / D))

# Preparar una función que realice el ajuste para un TF dado y devuelva los parámetros y el ajuste
def fit_logistic_curve(data, gene):
    # Filtrar datos para el TF dado
    gene_data = data[data['Gene'] == gene]
    
    # Asegurarse de que haya suficientes puntos de datos para el ajuste
    if gene_data.shape[0] < 4:
        return None, None, None
    
    # Extraer las concentraciones y los ratios
    x = gene_data['conc'].values
    y = gene_data['ratio'].values

    # Ajustar la curva logística
    try:
        # Estimar valores iniciales para los parámetros
        p0_estimates = [min(y), np.median(x), max(y), 1]
        params_opt, params_cov = curve_fit(logistic_function, x, y, p0=p0_estimates, maxfev=10000)
    except RuntimeError:
        return None, None, None
    
    # Calcular los valores ajustados
    x_fit = np.linspace(min(x), max(x), 500)
    y_fit = logistic_function(x_fit, *params_opt)
    
    return x_fit, y_fit, params_opt

# Crear gráficos para cada TF
unique_genes = data['Gene'].unique()

# Configurar la figura y los ejes
fig, axes = plt.subplots(ncols=8, figsize=(30,30))

# Ajustar la curva para cada TF y plotear
for ax, gene in zip(axes, unique_genes):
    x_fit, y_fit, params_opt = fit_logistic_curve(data, gene)

    #add the parameters to the dataset
    
    # Si el ajuste fue exitoso, plotear los datos
    if x_fit is not None:
        data.loc[data['Gene'] == gene, 'A'] = params_opt[0]
        data.loc[data['Gene'] == gene, 'B'] = params_opt[1]
        data.loc[data['Gene'] == gene, 'C'] = params_opt[2]
        data.loc[data['Gene'] == gene, 'D'] = params_opt[3]
    else:
        data.loc[data['Gene'] == gene, 'A'] = np.nan
        data.loc[data['Gene'] == gene, 'B'] = np.nan
        data.loc[data['Gene'] == gene, 'C'] = np.nan
        data.loc[data['Gene'] == gene, 'D'] = np.nan

# Ajustar el layout
data.head()

# Save the data
data.to_csv('phosphorylated_prots_fitted.tsv', sep='\t', index=False)



from scipy.optimize import curve_fit

# Definir la función logística para el ajuste
def logistic_equation(x, A, B, C):
    return A + (C - A) / (1 + np.exp(-(x - B)))

# Generar algunos datos simulados para el ajuste
np.random.seed(0)  # Semilla para reproducibilidad
x_data = np.linspace(0, 100, 50)  # Valores de x
y_true = logistic_equation(x_data, A=1, B=50, C=100)  # Valores verdaderos de y
y_data = y_true + np.random.normal(scale=5, size=x_data.size)  # Datos con ruido

# Realizar el ajuste de la curva
popt, pcov = curve_fit(logistic_equation, x_data, y_data, p0=[1, 50, 100])

# Usar los parámetros optimizados para calcular la y ajustada
y_fit = logistic_equation(x_data, *popt)

# Crear la gráfica
plt.figure(figsize=(10, 6))
plt.scatter(x_data, y_data, label='Noise', color='blue')
plt.plot(x_data, y_fit, label='Adjusted curve', color='red', linewidth=2)
plt.xlabel('Concentration')
plt.ylabel('Ratio')
plt.legend()
plt.show()
