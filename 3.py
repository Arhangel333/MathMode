import numpy as np
import time
import checks

# ============================================================
# 1. ИСХОДНЫЕ ДАННЫЕ
# ============================================================
a = 0.0            # левая граница отрезка
b = np.pi / 2            # правая граница отрезка
eps = 1e-6         # требуемая точность (можно изменить, если в задании другая)
max_iter = 1000    # максимальное число итераций (защита от бесконечного цикла)

beta, M, k = 18, 2, 1.4
beta = np.radians(beta)

def f(teta):
    return np.sin(teta)**2 - ((k+1)/2 * np.sin(teta)*np.sin(beta))/np.cos(teta-beta) - (1/M)**2

# ============================================================
# 2. МЕТОД ПОЛОВИННОГО ДЕЛЕНИЯ (МПД)
# ============================================================
def bisection_method(func, a, b, eps, max_iter):
    """
    Находит корень уравнения f(x) = 0 методом половинного деления.
    
    Параметры:
        func     - функция, корень которой ищем
        a, b     - границы отрезка [a, b]
        eps      - требуемая точность
        max_iter - максимальное число итераций
        
    Возвращает:
        root     - найденный корень
        n_iter   - количество итераций
        history  - список итераций для отчета (a, b, c, f(c))
    """
    # Проверка условия применимости: f(a) и f(b) должны иметь разные знаки
    if func(a) * func(b) > 0:
        raise ValueError("Ошибка: f(a) и f(b) имеют одинаковый знак. Корень на отрезке не гарантирован.")
    
    history = []
    
    for n_iter in range(1, max_iter + 1):
        c = (a + b) / 2.0          # середина отрезка
        fc = func(c)
        
        # Сохраняем историю для отчета (первые 5 итераций и последнюю)
        if n_iter <= 5 or (b - a) / 2 < eps:
            history.append((n_iter, a, b, c, fc))
        
        # Проверка достижения точности (полуширина отрезка меньше eps)
        if (b - a) / 2.0 < eps:
            return c, n_iter, history
        
        # Сужаем отрезок
        if func(a) * fc < 0:
            b = c    # корень в левой половине [a, c]
        else:
            a = c    # корень в правой половине [c, b]
            
    return c, max_iter, history


# --- Вспомогательная функция: поиск рабочего отрезка для МПД ---
def find_valid_interval_for_bisection(start, end, steps=1000):
    """Сканирует интервал и находит подынтервал, где f(x) меняет знак"""
    x_vals = np.linspace(start, end, steps)
    for i in range(len(x_vals) - 1):
        if f(x_vals[i]) * f(x_vals[i+1]) < 0:
            return x_vals[i], x_vals[i+1]
    return None, None


# ============================================================
# 3. ЗАПУСК И ВЫВОД РЕЗУЛЬТАТОВ
# ============================================================
if __name__ == "__main__":
    print("=" * 75)
    print("ЗАДАНИЕ №3: Решение нелинейного уравнения методом половинного деления")
    print("=" * 75)
    print(f"Уравнение: f(x) = sin(teta)**2 - ((k+1)/2 * sin(teta)*sin(beta))/cos(teta-beta) - (1/M)**2")
    print(f"Отрезок локализации корня: [{a}, {b}]")
    print(f"Требуемая точность: ε = {eps}")
    print("-" * 75)
    
    # Проверка условия перед запуском
    print(f"Проверка условия: f({a}) = {f(a):.4f}, f({b}) = {f(b):.4f}")
    print(f"f({a}) * f({b}) = {f(a)*f(b):.4f} < 0  \n")
    if f(a)*f(b) >= 0:
        a, b = find_valid_interval_for_bisection(a, b, steps=1000)
        print(f"Отрезок локализации корня: [{a}, {b}]")
        print(f"Проверка условия: f({a}) = {f(a):.4f}, f({b}) = {f(b):.4f}")
        print(f"f({a}) * f({b}) = {f(a)*f(b)} < 0  \n")
    if a is None:
        print("Ошибка: не удалось найти интервал со сменой знака!")
        exit()

    checks.check_root_existence(f, a, b, " (уравнение КСУ)")
    checks.check_uniqueness(f, a, b)
    
    # Вызов метода
    start = time.perf_counter()
    root, n_iter, history = bisection_method(f, a, b, eps, max_iter)
    total_time = time.perf_counter() - start

    print(f"Решение найдено, количесво итераций: {n_iter}")
    print(f"Время рассчёта: {total_time} сек")
    # Вывод таблицы первых итераций (для отчета)
    print("Ход решения (первые 5 итераций):")
    print(f"{'№':<3} | {'a':<10} | {'b':<10} | {'c (середина)':<15} | {'f(c)':<15}")
    print("-" * 75)
    for row in history[:5]:
        print(f"{row[0]:<3} | {row[1]:<10.6f} | {row[2]:<10.6f} | {row[3]:<15.6f} | {row[4]:<15.6e}")
    print("-" * 75)
    
    # Итоговый результат
    print("\nИТОГОВЫЙ РЕЗУЛЬТАТ:")
    print(f"Найденный корень:      x = {root:.8f}")
    print(f"Количество итераций:   n = {n_iter}")
    print(f"Проверка (невязка):    f(x) = {f(root):.2e}")
    print(f"Погрешность отрезка:   (b-a)/2 = {(history[-1][2]-history[-1][1])/2:.2e} < {eps}")
    print("=" * 75)
