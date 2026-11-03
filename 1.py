import numpy as np
import time
import matplotlib.pyplot as plt
import checks


DX = 1e-6 # приближение для производной
MNEXP = 1e-10 # приближение для МН
N_MAX_ITER = 100
def f(teta, beta, M, k):
    return np.sin(teta)**2 - ((k+1)/2 * np.sin(teta)*np.sin(beta))/np.cos(teta-beta) - (1/M)**2
           
def df(teta, beta, M, k):
    return (f(teta + DX, beta, M, k) - f(teta - DX, beta, M, k)) / (2 * DX)

def d2f(teta, beta=np.radians(18), k=1.4):
    """Вторая производная f''(θ_c) — аналитическая формула"""
    term1 = 2 * np.cos(2 * teta)
    term2 = (k + 1) * np.sin(beta) * np.cos(beta) * np.sin(teta - beta) / (np.cos(teta - beta)**3)
    return term1 - term2

def newton_method(beta, M, k=1.4, tol=MNEXP, max_iter=N_MAX_ITER):
    teta = np.arcsin(1.0 / M)
    
    for n_iter in range(1, max_iter + 1):
        f_val = f(teta, beta, M, k)
        df_val = df(teta, beta, M, k)
        
        if abs(df_val) < 1e-15:
            print("Ошибка: производная близка к нулю!")
            break
        
        teta_new = teta - f_val / df_val
        
        if abs(teta_new - teta) < tol:
            teta = teta_new
            break
        
        teta = teta_new
    
    return teta, n_iter






k = 1.4 #показателем адиабаты 

min_interval = 0
max_interval = np.pi / 2

#Использовать метод и обозначенные индексом 1 параметры набегающего потока, приведенные для каждого варианта в таблице 1 (приложение А). 
P_1 = 63950
T_1 = 320 
M_1 = 2
beta_1 = 18                 #degrees
beta_1 = np.radians(beta_1) #radians
# Method МН – метод Ньютона.
ro_1 = 1.225    
    
    
teta0 = np.arcsin(1.0 / M_1)
checks.check_newton_conditions(f, df, d2f, min_interval, max_interval, teta0, beta_1, M_1, k)
    
    
    
start = time.perf_counter()
teta, niter = newton_method(beta_1, M_1, k)
total_time = time.perf_counter() - start
print(f"Решение найдено, количесво итераций: {niter}")
print(f"Время рассчёта: {total_time} сек")
   
# 2) Рассчитать параметры 

R = 287.05
a1 = np.sqrt(k * R * T_1)  # скорость звука ≈ 340.3 м/с
v1 = M_1 * a1 
v1n = v1 * np.sin(teta)   
u1 = 1 - ((np.sin(teta-beta_1)*np.cos(teta))/(np.sin(teta)*np.cos(teta-beta_1)))

P_2 = P_1 + ro_1*u1*v1n**2
ro_2 = ro_1 / (1 - u1)
T_2 = P_2 / (R * ro_2)

# 1. Тангенциальная скорость (не меняется)
v1_tau = v1 * np.cos(teta)
v2_tau = v1_tau
# 2. Нормальная скорость за скачком (из уравнения неразрывности)
v2n = v1n * (1 - u1)
# 3. Полная скорость за скачком
v2 = np.sqrt(v2n**2 + v2_tau**2)
# 4. Скорость звука за скачком (через найденные p_2 и ro_2)
a2 = np.sqrt(k * P_2 / ro_2)
# 5. Число Маха за скачком (по формуле из конспекта: M = v / a)
M_2 = v2 / a2

# Шаг 6: Определение режима течения
if M_2 < 1.0:
    regime = "Дозвуковой"
else:
    regime = "Сверхзвуковой"

# Проверка, что угол лежит в допустимом интервале
if teta < min_interval or teta > max_interval:
    print("Ошибка: угол скачка вне допустимого диапазона [0, π/2]!")
else:
    print(f"Угол скачка θ_c = {np.degrees(teta):.4f}°")
    
    
print("Рассчитать параметры (p_2,T_2,ρ_2,M_2) за КСУ.")
print(f"Параметры : P_2={P_2}, T_2={T_2}, ro_2={ro_2}, M_2={M_2}, {regime}")
print(f"Угол скачка θ_c = {teta}")

""" 

	определить угол наклона КСУ θ_с из уравнения (4) на интервале (0, π/2), результат привести в градусах. 
	Рассчитать параметры (p_2,T_2,ρ_2,M_2) за КСУ. 
	Определить режим течения после КСУ по полученному значению числа Маха M_2 и указать это в отчете. 
	Привести все полученные значения в одну таблицу с размерностями параметров.
 """
 
 
teta_vals = np.linspace(0.01, np.pi/2 - 0.01, 500)
f_vals = f(teta_vals, beta_1, M_1, k)
root_shock = 0.868914
teta0 = np.arcsin(1.0 / M_1)

plt.figure(figsize=(10, 6))
plt.plot(np.degrees(teta_vals), f_vals, 'b-', linewidth=2, label=r'$f(\theta_c)$')
plt.axhline(0, color='black', linewidth=1.5, linestyle='--') # Ось X (y=0)

# Отмечаем корень и начальное приближение
plt.plot(np.degrees(root_shock), 0, 'ro', markersize=8, label=f'Корень: {np.degrees(root_shock):.2f}°')
plt.plot(np.degrees(teta0), f(teta0, beta_1, M_1, k), 'go', markersize=8, label=f'Нач. приближение (угол Маха): {np.degrees(teta0):.2f}°')

plt.title('График функции для определения угла наклона КСУ', fontsize=14)
plt.xlabel('Угол $\\theta_c$, градусы', fontsize=12)
plt.ylabel('Значение функции $f(\\theta_c)$', fontsize=12)
plt.legend(loc='upper right')
plt.xlim(0, 90)
plt.tight_layout()
plt.savefig('graph_shock.png', dpi=300) # Сохраняет картинку в папку с кодом