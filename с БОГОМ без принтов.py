import numpy as np
import copy
import time
import os

m, n = 103, 3
T1, T2 = 10, 40
Z, Zi = 5000, 300
Pk, Pm = 1.0, 1.0


def generator_heterogeneous_matrix(m_number, n_number, T1_number, T2_number):
    # Генерация массива случайных чисел
    mas = np.random.randint(T1_number, T2_number + 1, size=m_number * n_number)
    # Разделение массива на m подмассивов
    splits = [_.tolist() for _ in np.array_split(mas, m_number)]

    # Запись матрицы в файл
    with open(os.path.join(folder_path, f'main_matrix.txt'), 'w') as f:
        for row in splits:
            f.write(' '.join(map(str, row)) + '\n')


def find_scope_matrix():
    for _ in range(n):
        start = _ * 255 // n
        end = (_ + 1) * 255 // n
        if _ == 0:
            scope_matrix[_ * 2] = start
        else:
            scope_matrix[_ * 2] = start + 1
        scope_matrix[_ * 2 + 1] = end


def find_minimax_criterion():
    minimax_criterion_number_distribution = np.zeros((m, n), dtype=np.int32)  # распределение чисел по приборам
    minimax_criterion_one_dimensional_array_loads = np.zeros(n, dtype=np.int32)  # одномерный массив нагрузки

    print(f'\nодномерный массив нагрузки - {minimax_criterion_one_dimensional_array_loads}\n')

    for i_m in range(m):
        min_val = np.iinfo(np.int32).max  # Используем максимально возможное значение для int32
        min_index = -1

        for j_n in range(n):
            # переменная обхода
            traversal_variable = original_matrix[i_m, j_n] + minimax_criterion_one_dimensional_array_loads[j_n]

            # Нахождение минимального элемента в строке и его индекса
            if traversal_variable < min_val:
                min_val, min_index = traversal_variable, j_n

        # Обновляем сумму минимальных элементов по их индексам в matrix_min_element_method
        print(f'\n{minimax_criterion_one_dimensional_array_loads}\n')
        minimax_criterion_one_dimensional_array_loads[min_index] = np.int32(min_val)
        # записываем числа для просмотра распределения
        minimax_criterion_number_distribution[i_m, min_index] = np.int32(original_matrix[i_m, min_index])

        for k_z in range(Z):
            random_number = np.random.randint(scope_matrix[min_index * 2], scope_matrix[min_index * 2 + 1] + 1)
            print(f'{k_z + 1} - {random_number}')
            second_line_individual[k_z, i_m] = random_number
            print(second_line_individual)

    print(f'\nMAX из массива нагрузки --- {np.max(minimax_criterion_one_dimensional_array_loads)}')
    print(minimax_criterion_number_distribution)


def find_extra_minimax_criterion(barrier_min_numbers: int):
    index_for_extra, end_index_need_matrix = 0, 0

    extra_minimax_criterion_one_dimensional_array_loads = np.zeros(n, dtype=np.int32)  # одномерный массив нагрузки
    extra_minimax_criterion_additional_matrix = np.zeros((m, n), dtype=np.int32)  # дополнительная матрица
    extra_minimax_criterion_number_distribution = np.zeros((m, n), dtype=np.int32)  # распределение чисел по приборам
    extra_minimax_criterion_i_m_indices = []  # сохранение выпадающих индексов для построения генотипов

    for i_m in range(m):
        min_val = np.iinfo(np.int32).max  # Используем максимально возможное значение для int32
        min_index = -1

        for j_n in range(n):
            if original_matrix[i_m, j_n] < min_val:
                min_val = original_matrix[i_m, j_n]
                min_index = j_n

        traversal_variable_barrier = extra_minimax_criterion_one_dimensional_array_loads[min_index] + min_val

        if traversal_variable_barrier < barrier_min_numbers:
            extra_minimax_criterion_one_dimensional_array_loads[min_index] = traversal_variable_barrier
            extra_minimax_criterion_number_distribution[index_for_extra, min_index] = min_val

            for k_z in range(Z):
                random_number = np.random.randint(scope_matrix[min_index * 2], scope_matrix[min_index * 2 + 1] + 1)
                # закидываю по изначальному индексу i_m, чтобы распределение чисел в особи соответствовало
                # матрице, которая отсортирована в порядке убывания сумм строк
                second_line_individual[k_z, i_m] = random_number

            index_for_extra += 1
        else:
            extra_minimax_criterion_additional_matrix[end_index_need_matrix] = original_matrix[i_m]
            extra_minimax_criterion_i_m_indices.append(i_m)  # сохраняем пропущенный индекс
            end_index_need_matrix += 1

    for i_e in range(end_index_need_matrix):
        min_val = np.iinfo(np.int32).max  # Используем максимально возможное значение для int32
        min_index = -1

        for j_n in range(n):
            traversal_variable = extra_minimax_criterion_additional_matrix[i_e, j_n] + \
                                 extra_minimax_criterion_one_dimensional_array_loads[j_n]
            if traversal_variable < min_val:
                min_val = traversal_variable
                min_index = j_n

        extra_minimax_criterion_one_dimensional_array_loads[min_index] = min_val
        extra_minimax_criterion_number_distribution[index_for_extra, min_index] = \
            extra_minimax_criterion_additional_matrix[i_e, min_index]

        for k_z in range(Z):
            random_number = np.random.randint(scope_matrix[min_index * 2], scope_matrix[min_index * 2 + 1] + 1)
            second_line_individual[k_z, extra_minimax_criterion_i_m_indices[i_e]] = random_number

        index_for_extra += 1

    print(f'\nMAX из массива нагрузки --- {np.max(extra_minimax_criterion_one_dimensional_array_loads)}')
    print(extra_minimax_criterion_number_distribution)


def kernel_phenotype_detection(first_number: int, second_number: int, third_number: int):
    extra_minimax_criterion_number_distribution = np.zeros((first_number, second_number, third_number), dtype=np.int32)
    list_phenotype_detection = np.zeros((n,), dtype=np.int32)  # вычисление фенотипов особей
    for i_fn in range(first_number):          # Z
        for j_sn in range(second_number):     # m
            for k_tn in range(third_number):  # n
                if scope_matrix[k_tn * 2] <= second_line_individual[i_fn, j_sn] <= scope_matrix[k_tn * 2 + 1]:
                    list_phenotype_detection[k_tn] += original_matrix[j_sn, k_tn]
                    extra_minimax_criterion_number_distribution[i_fn, j_sn, k_tn] = original_matrix[j_sn, k_tn]
                    break

        max_value = np.max(list_phenotype_detection)
        individual_phenotypes[i_fn] = max_value
        list_phenotype_detection[:] = 0
    #print(extra_minimax_criterion_number_distribution[-1])


def randomize_individual_random_crossover(i_value: int, z_value: int, top_counter_value: int, slice_number_value: int):  # выбираем рандомную особь

    if top_counter_value % 10 < 5 and slice_number_value != -1:
        if i_value < slice_number_value:
            individual_random_crossover = np.random.randint(slice_number_value, z_value)
            while individual_random_crossover == i_value:
                individual_random_crossover = np.random.randint(slice_number_value, z_value)
            return individual_random_crossover
        else:
            individual_random_crossover = np.random.randint(0, slice_number_value)
            while individual_random_crossover == i_value:
                individual_random_crossover = np.random.randint(0, slice_number_value)
            return individual_random_crossover
    else:
        individual_random_crossover = np.random.randint(0, z_value)
        while individual_random_crossover == i_value:
            individual_random_crossover = np.random.randint(0, z_value)
        return individual_random_crossover


def randomize_number_random_crossover(m_value: int):  # выбираем рандомное число для кроссовера
    number_random_crossover = np.random.randint(0, m_value)
    while number_random_crossover == 0 or number_random_crossover == m_value:
        number_random_crossover = np.random.randint(0, m_value)
    return number_random_crossover


def kernel_crossover_detection(i_number: int, irc: int, nrc: int):  # кроссовер
    individual_random_crossover = irc
    number_random_crossover = nrc

    for j in range(m):
        if j < number_random_crossover:
            two_individual_second_line[0, j] = second_line_individual[i_number, j]
            two_individual_second_line[1, j] = second_line_individual[individual_random_crossover, j]
        else:
            two_individual_second_line[0, j] = second_line_individual[individual_random_crossover, j]
            two_individual_second_line[1, j] = second_line_individual[i_number, j]


def kernel_phenotype_detection_in_crossover_mutation(second_number: int, third_number: int):  # поиск фенотипов особей
    list_phenotype_detection_in_crossover_mutation = np.zeros(n, dtype=np.int32)  # вычисление фенотипов особей
    for i_fn in range(2):  # особи
        for j_sn in range(second_number):  # m
            for k_tn in range(third_number):  # n
                if scope_matrix[k_tn * 2] <= two_individual_second_line[i_fn, j_sn] <= scope_matrix[k_tn * 2 + 1]:
                    list_phenotype_detection_in_crossover_mutation[k_tn] += original_matrix[j_sn, k_tn]
                    break

        max_value = np.max(list_phenotype_detection_in_crossover_mutation)
        individual_phenotypes_in_crossover_mutation[i_fn] = max_value
        list_phenotype_detection_in_crossover_mutation[:] = 0


def kernel_mutation_detection():  # нахождение мутации
    list_phenotype_detection_in_crossover_mutation = np.zeros(n, dtype=np.int32)  # вычисление фенотипов особей
    for gene in range(2):

        random_gene = np.random.randint(m)  # выбираем рандомное число - ген
        random_gene_bit = np.random.randint(8)  # выбираем рандомное число - бит гена
        number_random_gene = two_individual_second_line[gene, random_gene]  # выбранное число по индексам

        binary_number = bin(number_random_gene)[2:].zfill(8)  # 8-битный формат
        # Заменяем выбранный бит на его противоположный
        binary_number_list = list(binary_number)
        binary_number_list[random_gene_bit] = '1' if binary_number_list[random_gene_bit] == '0' else '0'
        binary_number = ''.join(binary_number_list)

        # Переводим обратно в десятичную систему
        number_random_gene = int(binary_number, 2)

        new_number_random_gene, old_number_random_gene = 0, 0

        # проверка диапазона и пересчет фенотипа
        for f_number in range(n):
            if scope_matrix[f_number * 2] <= number_random_gene <= scope_matrix[f_number * 2 + 1]:
                new_number_random_gene = f_number
                break

        for f_number in range(n):
            if scope_matrix[f_number * 2] <= two_individual_second_line[gene, random_gene] <= scope_matrix[f_number * 2 + 1]:
                old_number_random_gene = f_number
                break

        #  перерасчет, если диапазоны не сходятся
        if new_number_random_gene != old_number_random_gene:

            two_individual_second_line[gene, random_gene] = number_random_gene  # замена числа в строке генов

            for pj_number in range(m):
                for pk_number in range(n):
                    if scope_matrix[pk_number * 2] <= two_individual_second_line[gene, pj_number] <= scope_matrix[pk_number * 2 + 1]:
                        list_phenotype_detection_in_crossover_mutation[pk_number] += original_matrix[pj_number, pk_number]
                        break

            individual_phenotypes_in_crossover_mutation[gene] = np.max(list_phenotype_detection_in_crossover_mutation)
            list_phenotype_detection_in_crossover_mutation[:] = 0
        else:
            two_individual_second_line[gene, random_gene] = number_random_gene  # замена числа в строке генов


def kernel_parent_child_comparison(i_number: int):  # сравнение ребенка и родителя

    min_index = np.argmin(individual_phenotypes_in_crossover_mutation)
    min_value = int(individual_phenotypes_in_crossover_mutation[min_index])

    if individual_phenotypes[i_number] < min_value:
        parent_child_comparison[i_number] = int(individual_phenotypes[i_number])
        second_line_individual_for_rounds[i_number] = copy.deepcopy(second_line_individual[i_number].astype(np.int32))
    else:
        parent_child_comparison[i_number] = min_value
        second_line_individual_for_rounds[i_number] = copy.deepcopy(two_individual_second_line[min_index].astype(np.int32))


desktop_path = os.path.join(os.path.join(os.environ['USERPROFILE']), 'Desktop')  # Путь к папке на рабочем столе
folder_name = "неоднородные матрицы"  # Имя папки, в которой будут храниться матрицы
folder_path = os.path.join(desktop_path, folder_name)  # Полный путь к папке
# Проверяем существование папки и создаем ее, если не существует
if not os.path.exists(folder_path):
    os.makedirs(folder_path)

original_matrix = np.zeros((m, n), dtype=np.int32)  # неоднородная матрица
scope_matrix = np.zeros(2 * n, dtype=np.int32)  # диапазон чисел

# начальные переменные особей, их не изменять по ходу программы
second_line_individual = np.zeros((Z, m), dtype=np.int32)  # вторая линия особей поколения
second_line_individual_for_rounds = np.zeros((Z, m), dtype=np.int32)  # для обработки внутри кода
two_individual_second_line = np.zeros((2, m), dtype=np.int32)  # две новых особи для кроссовера и мутации

slice_number = -1

choice_matrix = int(input('(1) -> Генерация новой матрицы\n(2) -> Матрица задана\n-> '))

if choice_matrix == 1:
    generator_heterogeneous_matrix(m, n, T1, T2)

matrix_of_all_iterations = np.genfromtxt(os.path.join(folder_path, f'main_matrix.txt'), delimiter=' ', dtype=np.int32)
print(f'\n--- Сгенерированная матрица ---\n{matrix_of_all_iterations}\n')
original_matrix = copy.deepcopy(np.array(sorted(matrix_of_all_iterations, key=sum, reverse=True)))
print(f'\n--- Отсортированная ---\n')
for _ in original_matrix:
    print(_, sum(_))
find_scope_matrix()

print(f'диапазон - {scope_matrix}')

choice_method = int(input('(1) -> Рандомные особи\n(2) -> Генерация особей с помощью Плотникова-Зверева\n-> '))

if choice_method == 1:
    second_line_individual = np.random.randint(0, 256, (Z, m)).astype(np.int32)

if choice_method == 2:

    choice_algorithm = int(input('Какой начальный алгоритм Вы хотите выбрать?\n(1) -> П-З по минимаксному критерию\n(2) -> Экстра П-З по минимаксному критерию\n-> '))

    if choice_algorithm == 1:
        find_minimax_criterion()  # подсчет минимаксного критерия

    if choice_algorithm == 2:
        min_elements = np.min(original_matrix, axis=1)  # находим минимальные элементы в каждой строке
        sum_of_min_elements = (np.sum(min_elements) // n) + 1  # находим барьер
        find_extra_minimax_criterion(sum_of_min_elements)

    choice_individual = int(input('Подмешивать особи?\n(1) -> да\n(2) -> нет\n-> '))

    if choice_individual == 1:
        choice_slice_number = int(input('Какое % соотношение?\n-> '))
        slice_number = int(Z * choice_slice_number / 100)
        second_line_individual[:slice_number] = np.random.randint(0, 256, (slice_number, m)).astype(np.int32)

top_counter = np.int32(0)  # счетчик лучших особей для остановки при достаточном повторе
counter_individual = np.int32(0)  # подсчет кол-ва пройденных поколений
best_phenotype_for_top_counter = np.int32(0)  # лучший фенотип для повторов особи

start_time = time.perf_counter()

while top_counter < Zi:

    individual_phenotypes = np.zeros((Z,), dtype=np.int32)  # max числа из подсчета фенотипов особей одного поколения
    individual_phenotypes_in_crossover_mutation = np.zeros(2, dtype=np.int32)  # числа фенотипов детей
    parent_child_comparison = np.zeros(Z, dtype=np.int32)  # итоговые фенотипы в конце нового поколения

    if top_counter == 0:  # чтобы изначально top_counter был единицей, т.к. в первом поколении уже есть лучший фенотип, если во втором будет такой же, то это уже будет второй повтор
        top_counter = np.int32(1)
        counter_individual = np.int32(1)

    kernel_phenotype_detection(Z, m, n)

    best_phenotype_for_top_counter = np.min(individual_phenotypes)

    if counter_individual == 1:
        print(f'\nПоколение №{counter_individual}  --- лучший фенотип - {best_phenotype_for_top_counter} --- кол-во повторов - {top_counter}')

    for i in range(Z):

        random_number_pk = np.float16(np.random.uniform(0.0, 1.0))
        random_number_pm = np.float16(np.random.uniform(0.0, 1.0))

        if random_number_pk >= Pk and random_number_pm >= Pm:

            second_line_individual_for_rounds[i] = copy.deepcopy(second_line_individual[i])
            parent_child_comparison[i] = int(individual_phenotypes[i])

        else:

            if random_number_pk < Pk:

                first_irc = randomize_individual_random_crossover(i, Z, top_counter, slice_number)
                second_nrc = randomize_number_random_crossover(m)

                kernel_crossover_detection(i, first_irc, second_nrc)
                kernel_phenotype_detection_in_crossover_mutation(m, n)

            if random_number_pk >= Pk:

                first_irc = randomize_individual_random_crossover(i, Z, top_counter, slice_number)

                two_individual_second_line[0] = copy.deepcopy(second_line_individual[i])
                two_individual_second_line[1] = copy.deepcopy(second_line_individual[first_irc])

                kernel_phenotype_detection_in_crossover_mutation(m, n)

            if random_number_pk < Pm:

                kernel_mutation_detection()

            kernel_parent_child_comparison(i)

    second_line_individual = copy.deepcopy(second_line_individual_for_rounds)
    individual_phenotypes = copy.deepcopy(parent_child_comparison)

    min_value_of_phenotype = np.min(individual_phenotypes)
    if best_phenotype_for_top_counter == min_value_of_phenotype:
        top_counter += 1
    else:
        best_phenotype_for_top_counter = min_value_of_phenotype
        top_counter = 1

    counter_individual += 1

    print(f'Поколение №{counter_individual} --- лучший фенотип - {best_phenotype_for_top_counter} --- кол-во повторов - {top_counter}')

    second_line_individual_for_rounds[:], two_individual_second_line[:] = 0, 0

end_time = time.perf_counter()
elapsed_time = end_time - start_time
hours, rem = divmod(elapsed_time, 3600)
minutes, seconds = divmod(rem, 60)

print(f"\n{int(hours):02}:{int(minutes):02}:{seconds:05.2f}")
