import streamlit as st
import numpy as np
import copy
import time


def find_scope_matrix(n, scope_matrix):
    for i in range(n):
        start = i * 255 // n
        end = (i + 1) * 255 // n
        if i == 0:
            scope_matrix[i * 2] = start
        else:
            scope_matrix[i * 2] = start + 1
        scope_matrix[i * 2 + 1] = end


def kernel_phenotype_detection(Z, m, n, scope_matrix, second_line_individual, original_matrix, individual_phenotypes):
    extra_minimax_criterion_number_distribution = np.zeros((Z, m, n), dtype=np.int32)
    list_phenotype_detection = np.zeros(n, dtype=np.int32)
    for i_fn in range(Z):
        for j_sn in range(m):
            for k_tn in range(n):
                if scope_matrix[k_tn * 2] <= second_line_individual[i_fn, j_sn] <= scope_matrix[k_tn * 2 + 1]:
                    list_phenotype_detection[k_tn] += original_matrix[j_sn, k_tn]
                    extra_minimax_criterion_number_distribution[i_fn, j_sn, k_tn] = original_matrix[j_sn, k_tn]
                    break
        max_value = np.max(list_phenotype_detection)
        individual_phenotypes[i_fn] = max_value
        list_phenotype_detection[:] = 0


def randomize_individual_random_crossover(i_value, z_value):
    individual_random_crossover = np.random.randint(0, z_value)
    while individual_random_crossover == i_value:
        individual_random_crossover = np.random.randint(0, z_value)
    return individual_random_crossover


def randomize_number_random_crossover(m_value):
    number_random_crossover = np.random.randint(0, m_value)
    while number_random_crossover == 0 or number_random_crossover == m_value:
        number_random_crossover = np.random.randint(0, m_value)
    return number_random_crossover


def kernel_crossover_detection(i_number, irc, nrc, m, second_line_individual, two_individual_second_line):
    individual_random_crossover = irc
    number_random_crossover = nrc
    for j in range(m):
        if j < number_random_crossover:
            two_individual_second_line[0, j] = second_line_individual[i_number, j]
            two_individual_second_line[1, j] = second_line_individual[individual_random_crossover, j]
        else:
            two_individual_second_line[0, j] = second_line_individual[individual_random_crossover, j]
            two_individual_second_line[1, j] = second_line_individual[i_number, j]


def kernel_phenotype_detection_in_crossover_mutation(m, n, scope_matrix, original_matrix, two_individual_second_line,
                                                     individual_phenotypes_in_crossover_mutation):
    list_phenotype_detection_in_crossover_mutation = np.zeros(n, dtype=np.int32)
    for i_fn in range(2):
        for j_sn in range(m):
            for k_tn in range(n):
                if scope_matrix[k_tn * 2] <= two_individual_second_line[i_fn, j_sn] <= scope_matrix[k_tn * 2 + 1]:
                    list_phenotype_detection_in_crossover_mutation[k_tn] += original_matrix[j_sn, k_tn]
                    break
        max_value = np.max(list_phenotype_detection_in_crossover_mutation)
        individual_phenotypes_in_crossover_mutation[i_fn] = max_value
        list_phenotype_detection_in_crossover_mutation[:] = 0


def kernel_mutation_detection(m, n, scope_matrix, original_matrix, two_individual_second_line,
                              individual_phenotypes_in_crossover_mutation):
    list_phenotype_detection_in_crossover_mutation = np.zeros(n, dtype=np.int32)
    for gene in range(2):
        random_gene = np.random.randint(m)
        random_gene_bit = np.random.randint(8)
        number_random_gene = two_individual_second_line[gene, random_gene]
        binary_number = bin(number_random_gene)[2:].zfill(8)
        binary_number_list = list(binary_number)
        binary_number_list[random_gene_bit] = '1' if binary_number_list[random_gene_bit] == '0' else '0'
        binary_number = ''.join(binary_number_list)
        number_random_gene = int(binary_number, 2)
        new_number_random_gene, old_number_random_gene = 0, 0
        for f_number in range(n):
            if scope_matrix[f_number * 2] <= number_random_gene <= scope_matrix[f_number * 2 + 1]:
                new_number_random_gene = f_number
                break
        for f_number in range(n):
            if scope_matrix[f_number * 2] <= two_individual_second_line[gene, random_gene] <= scope_matrix[
                f_number * 2 + 1]:
                old_number_random_gene = f_number
                break
        if new_number_random_gene != old_number_random_gene:
            two_individual_second_line[gene, random_gene] = number_random_gene
            for pj_number in range(m):
                for pk_number in range(n):
                    if scope_matrix[pk_number * 2] <= two_individual_second_line[gene, pj_number] <= scope_matrix[
                        pk_number * 2 + 1]:
                        list_phenotype_detection_in_crossover_mutation[pk_number] += original_matrix[
                            pj_number, pk_number]
                        break
            individual_phenotypes_in_crossover_mutation[gene] = np.max(list_phenotype_detection_in_crossover_mutation)
            list_phenotype_detection_in_crossover_mutation[:] = 0
        else:
            two_individual_second_line[gene, random_gene] = number_random_gene


def kernel_parent_child_comparison(i_number, individual_phenotypes, individual_phenotypes_in_crossover_mutation,
                                   second_line_individual, two_individual_second_line,
                                   second_line_individual_for_rounds, parent_child_comparison):
    min_index = np.argmin(individual_phenotypes_in_crossover_mutation)
    min_value = int(individual_phenotypes_in_crossover_mutation[min_index])
    if individual_phenotypes[i_number] < min_value:
        parent_child_comparison[i_number] = int(individual_phenotypes[i_number])
        second_line_individual_for_rounds[i_number] = copy.deepcopy(second_line_individual[i_number].astype(np.int32))
    else:
        parent_child_comparison[i_number] = min_value
        second_line_individual_for_rounds[i_number] = copy.deepcopy(
            two_individual_second_line[min_index].astype(np.int32))


def run_algorithm(Z, m, n, Zi, Pk, Pm, original_matrix, scope_matrix, progress_placeholder):
    second_line_individual = np.random.randint(0, 256, (Z, m)).astype(np.int32)
    two_individual_second_line = np.zeros((2, m), dtype=np.int32)
    second_line_individual_for_rounds = np.zeros((Z, m), dtype=np.int32)
    top_counter = np.int32(0)
    counter_individual = np.int32(0)
    best_phenotype_for_top_counter = np.int32(0)
    start_time = time.perf_counter()
    progress_log = st.session_state.progress_log  # Используем сессию для хранения логов
    while top_counter < Zi:
        individual_phenotypes = np.zeros((Z,), dtype=np.int32)
        individual_phenotypes_in_crossover_mutation = np.zeros(2, dtype=np.int32)
        parent_child_comparison = np.zeros(Z, dtype=np.int32)
        if top_counter == 0:
            top_counter = np.int32(1)
            counter_individual = np.int32(1)
        kernel_phenotype_detection(Z, m, n, scope_matrix, second_line_individual, original_matrix,
                                   individual_phenotypes)
        best_phenotype_for_top_counter = np.min(individual_phenotypes)
        if counter_individual == 1:
            progress_msg = f'Поколение №{counter_individual}  --- лучший фенотип - {best_phenotype_for_top_counter} --- кол-во повторов - {top_counter}'
            progress_placeholder.text(progress_msg)
            progress_log.append(progress_msg)
        for i in range(Z):
            random_number_pk = np.float16(np.random.uniform(0.0, 1.0))
            random_number_pm = np.float16(np.random.uniform(0.0, 1.0))
            if random_number_pk >= Pk and random_number_pm >= Pm:
                second_line_individual_for_rounds[i] = copy.deepcopy(second_line_individual[i])
                parent_child_comparison[i] = int(individual_phenotypes[i])
            else:
                if random_number_pk < Pk:
                    first_irc = randomize_individual_random_crossover(i, Z)
                    second_nrc = randomize_number_random_crossover(m)
                    kernel_crossover_detection(i, first_irc, second_nrc, m, second_line_individual,
                                               two_individual_second_line)
                    kernel_phenotype_detection_in_crossover_mutation(m, n, scope_matrix, original_matrix,
                                                                     two_individual_second_line,
                                                                     individual_phenotypes_in_crossover_mutation)
                if random_number_pk >= Pk:
                    first_irc = randomize_individual_random_crossover(i, Z)
                    two_individual_second_line[0] = copy.deepcopy(second_line_individual[i])
                    two_individual_second_line[1] = copy.deepcopy(second_line_individual[first_irc])
                    kernel_phenotype_detection_in_crossover_mutation(m, n, scope_matrix, original_matrix,
                                                                     two_individual_second_line,
                                                                     individual_phenotypes_in_crossover_mutation)
                if random_number_pk < Pm:
                    kernel_mutation_detection(m, n, scope_matrix, original_matrix, two_individual_second_line,
                                              individual_phenotypes_in_crossover_mutation)
                kernel_parent_child_comparison(i, individual_phenotypes, individual_phenotypes_in_crossover_mutation,
                                               second_line_individual, two_individual_second_line,
                                               second_line_individual_for_rounds, parent_child_comparison)
        second_line_individual = copy.deepcopy(second_line_individual_for_rounds)
        individual_phenotypes = copy.deepcopy(parent_child_comparison)
        min_value_of_phenotype = np.min(individual_phenotypes)
        if best_phenotype_for_top_counter == min_value_of_phenotype:
            top_counter += 1
        else:
            best_phenotype_for_top_counter = min_value_of_phenotype
            top_counter = 1
        counter_individual += 1
        progress_msg = f'Поколение №{counter_individual} --- лучший фенотип - {best_phenotype_for_top_counter} --- кол-во повторов - {top_counter}'
        progress_placeholder.text(progress_msg)
        progress_log.append(progress_msg)
        second_line_individual_for_rounds[:], two_individual_second_line[:] = 0, 0
    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    hours, rem = divmod(elapsed_time, 3600)
    minutes, seconds = divmod(rem, 60)
    final_msg = f"\n{int(hours):02}:{int(minutes):02}:{seconds:05.2f}"
    st.write(final_msg)
    progress_log.append(final_msg)

    # Сохраняем лог в файл
    log_file = "\n".join(progress_log)

    # Добавляем кнопку для скачивания файла
    st.download_button(
        label="Скачать лог",
        data=log_file,
        file_name="progress_log.txt",
        mime="text/plain"
    )

    st.session_state.progress_log = []


def app():
    st.title('Алгоритм Голдберга')

    if 'progress_log' not in st.session_state:
        st.session_state.progress_log = []

    matrix_option = st.selectbox(
        "Выберите источник матрицы:",
        ("", "Генерация новой матрицы", "Матрица задана"),
        key="matrix_option"
    )

    if matrix_option == "Генерация новой матрицы":

        m = st.number_input("Количество строк матрицы (m)", min_value=1, value=20, step=1)
        n = st.number_input("Количество столбцов матрицы (n)", min_value=1, value=3, step=1)
        T1 = st.number_input("Начало диапазона значений (T1)", min_value=1, value=10, step=1)
        T2 = st.number_input("Конец диапазона значений (T2)", min_value=T1 + 1, value=20, step=1)
        Z = st.number_input("Количество особей (Z)", min_value=100, value=1000, step=1)
        Zi = st.number_input("Количество повторов (Zi)", min_value=1, value=100, step=1)
        Pk = st.number_input("Вероятность кроссовера (Pk)", min_value=0.0, max_value=1.0, value=0.5, step=0.01)
        Pm = st.number_input("Вероятность мутации (Pm)", min_value=0.0, max_value=1.0, value=0.5, step=0.01)

        if st.button("Начать процесс"):
            original_matrix = np.random.randint(T1, T2 + 1, size=(m, n)).astype(np.int32)
            sorted_matrix = np.array(sorted(original_matrix, key=lambda x: sum(x), reverse=True))
            col1, col2 = st.columns(2)
            with col1:
                st.write("Сгенерированная матрица:")
                st.write(original_matrix)
            with col2:
                st.write("Отсортированная матрица:")
                st.write(sorted_matrix)
            scope_matrix = np.zeros(2 * n, dtype=np.int32)
            find_scope_matrix(n, scope_matrix)
            st.write(f'Диапазон - {scope_matrix}')
            st.session_state.progress_log.append(f'Диапазон - {scope_matrix}')
            progress_placeholder = st.empty()
            run_algorithm(Z, m, n, Zi, Pk, Pm, sorted_matrix, scope_matrix, progress_placeholder)

    elif matrix_option == "Матрица задана":
        uploaded_file = st.file_uploader("Загрузите файл матрицы", type="txt", key="file_uploader")

        if uploaded_file is not None:
            file_content = uploaded_file.getvalue().decode("utf-8").strip()
            try:
                matrix_of_all_iterations = np.loadtxt(file_content.splitlines(), dtype=int)
                m, n = matrix_of_all_iterations.shape  # Определяем размеры матрицы
                sorted_matrix = np.array(sorted(matrix_of_all_iterations, key=lambda x: sum(x), reverse=True))
                col1, col2 = st.columns(2)
                with col1:
                    st.write("Загруженная матрица:")
                    st.write(matrix_of_all_iterations)
                with col2:
                    st.write("Отсортированная матрица:")
                    st.write(sorted_matrix)

                Z = st.number_input("Количество особей (Z)", min_value=100, value=1000, step=1)
                Zi = st.number_input("Количество повторов (Zi)", min_value=1, value=100, step=1)
                Pk = st.number_input("Вероятность кроссовера (Pk)", min_value=0.0, max_value=1.0, value=0.5, step=0.01)
                Pm = st.number_input("Вероятность мутации (Pm)", min_value=0.0, max_value=1.0, value=0.5, step=0.01)

                if st.button("Начать процесс"):
                    scope_matrix = np.zeros(2 * n, dtype=np.int32)
                    find_scope_matrix(n, scope_matrix)
                    st.write(f'Диапазон - {scope_matrix}')
                    st.session_state.progress_log.append(f'Диапазон - {scope_matrix}')
                    progress_placeholder = st.empty()
                    run_algorithm(Z, m, n, Zi, Pk, Pm, sorted_matrix, scope_matrix, progress_placeholder)

            except ValueError:
                st.warning("Файл содержит недопустимые символы. Пожалуйста, загрузите файл, содержащий только числа.")


if __name__ == "__main__":
    app()
