import numpy as np




def check_root_existence(func, a, b, name=""):
    """
    Проверка условия Больцано-Коши (существование корня)
    """
    fa = func(a)
    fb = func(b)
    product = fa * fb
    
    print(f"\n{'='*60}")
    print(f"ПРОВЕРКА СУЩЕСТВОВАНИЯ КОРНЯ{name}")
    print(f"{'='*60}")
    print(f"Отрезок: [{a:.4f}, {b:.4f}]")
    print(f"f({a:.4f}) = {fa:.6f}")
    print(f"f({b:.4f}) = {fb:.6f}")
    print(f"f(a) · f(b) = {product:.6e}")
    
    if product < 0:
        print(" Условие Больцано-Коши ВЫПОЛНЕНО: корень существует")
        return True
    else:
        print(" Условие Больцано-Коши НЕ ВЫПОЛНЕНО: корень может отсутствовать")
        return False


def check_uniqueness(func, a, b, steps=100):
    """
    Проверка единственности корня (монотонность функции)
    Проверяем, что f'(x) не меняет знак на отрезке
    """
    print(f"\n{'='*60}")
    print("ПРОВЕРКА ЕДИНСТВЕННОСТИ КОРНЯ")
    print(f"{'='*60}")
    
    
    x_vals = np.linspace(a, b, steps)
    dx = 1e-6
    derivative_vals = []
    
    for x in x_vals:
        df = (func(x + dx) - func(x - dx)) / (2 * dx)
        derivative_vals.append(df)
    
    derivative_vals = np.array(derivative_vals)
    
    
    if np.all(derivative_vals > 0):
        print(" f'(x) > 0 на всём отрезке → функция возрастает")
        print(" Корень ЕДИНСТВЕННЫЙ")
        return True
    elif np.all(derivative_vals < 0):
        print(" f'(x) < 0 на всём отрезке → функция убывает")
        print(" Корень ЕДИНСТВЕННЫЙ")
        return True
    else:
        print(" f'(x) меняет знак → функция немонотонна")
        print(" Корней может быть несколько!")
        return False


def check_newton_conditions(func, df, d2f, a, b, x0, beta, M, k):
    """
    Проверка условий сходимости метода Ньютона:
    1. f'(x) ≠ 0 на отрезке
    2. f''(x) не меняет знак
    3. f(x0) · f''(x0) > 0 (правильный выбор начального приближения)
    """
    print(f"\n{'='*60}")
    print("ПРОВЕРКА УСЛОВИЙ СХОДИМОСТИ МЕТОДА НЬЮТОНА")
    print(f"{'='*60}")
    
    all_passed = True
    
    
    x_vals = np.linspace(a, b, 100)
    df_vals = [df(x, beta, M, k) for x in x_vals]
    
    if np.all(np.abs(df_vals) > 1e-10):
        print(" Условие 1: f'(x) ≠ 0 на отрезке [a, b]")
    else:
        print(" Условие 1: f'(x) обращается в ноль на отрезке!")
        all_passed = False
    
    
    d2f_vals = [d2f(x) for x in x_vals]
    
    if np.all(np.array(d2f_vals) > 0):
        print(" Условие 2: f''(x) > 0 (функция выпукла вниз)")
    elif np.all(np.array(d2f_vals) < 0):
        print(" Условие 2: f''(x) < 0 (функция выпукла вверх)")
    else:
        print(" Условие 2: f''(x) меняет знак на отрезке!")
        all_passed = False
    
    
    f_x0 = func(x0, beta, M, k)
    d2f_x0 = d2f(x0)
    product = f_x0 * d2f_x0
    
    print(f"\nПроверка начального приближения x0 = {x0:.4f}:")
    print(f"  f(x0) = {f_x0:.6f}")
    print(f"  f''(x0) = {d2f_x0:.6f}")
    print(f"  f(x0) · f''(x0) = {product:.6e}")
    
    if product > 0:
        print(" Условие 3: f(x0) · f''(x0) > 0 — начальное приближение ВЫБРАНО ВЕРНО")
    else:
        print(" Условие 3: f(x0) · f''(x0) < 0 — начальное приближение ВЫБРАНО НЕВЕРНО!")
        all_passed = False
    
    print(f"\n{'='*60}")
    if all_passed:
        print("ВСЕ УСЛОВИЯ СХОДИМОСТИ МЕТОДА НЬЮТОНА ВЫПОЛНЕНЫ ")
    else:
        print("НЕ ВСЕ УСЛОВИЯ ВЫПОЛНЕНЫ! Метод может не сойтись ")
    print(f"{'='*60}")
    
    return all_passed
