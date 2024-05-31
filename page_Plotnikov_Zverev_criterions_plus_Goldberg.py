import streamlit as st
import numpy as np
import copy
import time
import zipfile
import io


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


def kernel_phenotype_detection_and_save(Z, m, n, scope_matrix, second_line_individual, sorted_matrix,
                                        individual_phenotypes):
    extra_minimax_criterion_number_distribution = np.zeros((Z, m, n), dtype=np.int32)
    list_phenotype_detection = np.zeros(n, dtype=np.int32)
    phenotype_distribution_list = np.zeros((Z, n), dtype=np.int32)

    for i_fn in range(Z):
        for j_sn in range(m):
            for k_tn in range(n):
                if scope_matrix[k_tn * 2] <= second_line_individual[i_fn, j_sn] <= scope_matrix[k_tn * 2 + 1]:
                    list_phenotype_detection[k_tn] += sorted_matrix[j_sn, k_tn]
                    extra_minimax_criterion_number_distribution[i_fn, j_sn, k_tn] = sorted_matrix[j_sn, k_tn]
                    break
        max_value = np.max(list_phenotype_detection)
        individual_phenotypes[i_fn] = max_value
        phenotype_distribution_list[i_fn] = list_phenotype_detection
        list_phenotype_detection[:] = 0

    max_index = np.argmin(individual_phenotypes)
    distribution_list = extra_minimax_criterion_number_distribution[max_index]
    phenotype_detection_list = phenotype_distribution_list[max_index]

    distribution_file_content = f"- Распределение заданий по приборам -\n\n{distribution_list}\n\n{phenotype_detection_list}"
    return distribution_file_content


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


def initialize_algorithm(Z, m):
    two_individual_second_line = np.zeros((2, m), dtype=np.int32)
    second_line_individual_for_rounds = np.zeros((Z, m), dtype=np.int32)
    return two_individual_second_line, second_line_individual_for_rounds


def initialize_tracking():
    top_counter = np.int32(0)
    counter_individual = np.int32(0)
    best_phenotype_for_top_counter = np.int32(0)
    progress_log_Goldberg = []
    progress_placeholder = st.empty()
    return top_counter, counter_individual, best_phenotype_for_top_counter, progress_log_Goldberg, progress_placeholder


def track_progress(counter_individual, best_phenotype_for_top_counter, top_counter, progress_placeholder,
                   progress_log_Goldberg):
    progress_msg = f'Поколение №{counter_individual} --- лучший фенотип - {best_phenotype_for_top_counter} --- кол-во повторов - {top_counter}'
    progress_placeholder.text(progress_msg)
    progress_log_Goldberg.append(progress_msg)


def finalize_tracking(progress_log_Goldberg, start_time, end_time, second_line_individual, m, n, scope_matrix,
                      sorted_matrix, individual_phenotypes):
    if 'finalize_called' not in st.session_state:
        st.session_state['finalize_called'] = True
        elapsed_time = end_time - start_time
        hours, rem = divmod(elapsed_time, 3600)
        minutes, seconds = divmod(rem, 60)
        final_msg = f"\n{int(hours):02}:{int(minutes):02}:{seconds:05.2f}"
        progress_log_Goldberg.append(final_msg)

        progress_distribution_file = kernel_phenotype_detection_and_save(len(second_line_individual), m, n,
                                                                         scope_matrix,
                                                                         second_line_individual, sorted_matrix,
                                                                         individual_phenotypes)

        st.session_state['progress_log'] = progress_log_Goldberg
        st.session_state['distribution_file'] = progress_distribution_file
        st.session_state['final_msg'] = final_msg
        st.session_state['algorithm_completed'] = True

        zip_buffer = download_files_as_zip(progress_log_Goldberg, progress_distribution_file)
        st.session_state['zip_buffer'] = zip_buffer

        # Вывод последнего сообщения лога и итогового сообщения с временем
        st.write(progress_log_Goldberg[-2])  # Вывод предпоследнего сообщения, которое является последним результатом
        st.write(final_msg)

        st.download_button(
            label="Скачать все файлы",
            data=zip_buffer,
            file_name="results.zip",
            mime="application/zip",
            key="final_download_button"
        )


def run_algorithm(Z, m, n, Zi, Pk, Pm, sorted_matrix, scope_matrix, number_of_alg, second_line_individual):
    two_individual_second_line, second_line_individual_for_rounds = initialize_algorithm(Z, m)
    top_counter, counter_individual, best_phenotype_for_top_counter, progress_log_Goldberg, progress_placeholder = initialize_tracking()

    start_time = time.perf_counter()
    st_scope_matrix = f'Диапазон --- | {" | ".join(f"{scope_matrix[i]} - {scope_matrix[i + 1]}" for i in range(0, len(scope_matrix), 2))} |\n'
    progress_log_Goldberg.append(st_scope_matrix)

    individual_phenotypes = np.zeros((Z,), dtype=np.int32)

    while top_counter < Zi:
        individual_phenotypes_in_crossover_mutation = np.zeros(2, dtype=np.int32)
        parent_child_comparison = np.zeros(Z, dtype=np.int32)

        if top_counter == 0:
            top_counter = np.int32(1)
            counter_individual = np.int32(1)

        kernel_phenotype_detection(Z, m, n, scope_matrix, second_line_individual, sorted_matrix, individual_phenotypes)
        best_phenotype_for_top_counter = np.min(individual_phenotypes)

        if counter_individual == 1:
            track_progress(counter_individual, best_phenotype_for_top_counter, top_counter, progress_placeholder,
                           progress_log_Goldberg)

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
                    kernel_phenotype_detection_in_crossover_mutation(m, n, scope_matrix, sorted_matrix,
                                                                     two_individual_second_line,
                                                                     individual_phenotypes_in_crossover_mutation)
                if random_number_pk >= Pk:
                    first_irc = randomize_individual_random_crossover(i, Z)
                    two_individual_second_line[0] = copy.deepcopy(second_line_individual[i])
                    two_individual_second_line[1] = copy.deepcopy(second_line_individual[first_irc])
                    kernel_phenotype_detection_in_crossover_mutation(m, n, scope_matrix, sorted_matrix,
                                                                     two_individual_second_line,
                                                                     individual_phenotypes_in_crossover_mutation)
                if random_number_pk < Pm:
                    kernel_mutation_detection(m, n, scope_matrix, sorted_matrix, two_individual_second_line,
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
        track_progress(counter_individual, best_phenotype_for_top_counter, top_counter, progress_placeholder,
                       progress_log_Goldberg)

        second_line_individual_for_rounds[:], two_individual_second_line[:] = 0, 0

    end_time = time.perf_counter()
    finalize_tracking(progress_log_Goldberg, start_time, end_time, second_line_individual, m, n, scope_matrix,
                      sorted_matrix, individual_phenotypes)


def generate_matrix(m, n, T1, T2):
    original_matrix = np.random.randint(T1, T2 + 1, size=(m, n)).astype(np.int32)
    sorted_matrix = np.array(sorted(original_matrix, key=lambda x: sum(x), reverse=True))
    return original_matrix, sorted_matrix


def process_uploaded_file(uploaded_file):
    if uploaded_file is not None:
        file_content = uploaded_file.getvalue().decode("utf-8").strip()
        try:
            matrix_of_all_iterations = np.loadtxt(file_content.splitlines(), dtype=int)
            return matrix_of_all_iterations
        except ValueError:
            st.warning("Файл содержит недопустимые символы. Пожалуйста, загрузите файл, содержащий только числа.")
    return None


def download_files_as_zip(log_file_content, distribution_file_content):
    zip_buffer = io.BytesIO()
    with zipfile.ZipFile(zip_buffer, "a", zipfile.ZIP_DEFLATED) as zip_file:
        zip_file.writestr("progress_log.txt", "\n".join(log_file_content))  # Преобразуем список в строку
        zip_file.writestr("distribution_list.txt", distribution_file_content)  # Это уже строка
    zip_buffer.seek(0)
    return zip_buffer


def initialize_state():
    keys = ["generation_method_selected", "original_matrix", "sorted_matrix", "scope_matrix", "m", "n", "Z", "Zi", "Pk",
            "Pm", "generation_method", "algorithm_option", "second_line_individual"]
    for key in keys:
        if key not in st.session_state:
            st.session_state[key] = None
    if "initialized" not in st.session_state:
        st.session_state["initialized"] = True
    if "generation_method_selected" not in st.session_state:
        st.session_state["generation_method_selected"] = False
    if "progress_log_Goldberg" not in st.session_state:
        st.session_state["progress_log_Goldberg"] = []


def reset_state():
    keys = ["generation_method_selected", "original_matrix", "sorted_matrix", "scope_matrix", "m", "n", "Z", "Zi", "Pk",
            "Pm", "generation_method", "algorithm_option", "second_line_individual", "progress_log_Goldberg"]
    for key in keys:
        if key in st.session_state:
            del st.session_state[key]
    initialize_state()


def find_minimax_criterion(second_line_individual, scope_matrix, sorted_matrix, m, n, Z):
    minimax_criterion_number_distribution = np.zeros((m, n), dtype=np.int32)
    minimax_criterion_one_dimensional_array_loads = np.zeros(n, dtype=np.int32)

    for i_m in range(m):
        min_val = np.iinfo(np.int32).max
        min_index = -1

        for j_n in range(n):
            traversal_variable = sorted_matrix[i_m, j_n] + minimax_criterion_one_dimensional_array_loads[j_n]

            if traversal_variable < min_val:
                min_val, min_index = traversal_variable, j_n

        minimax_criterion_one_dimensional_array_loads[min_index] = np.int32(min_val)
        minimax_criterion_number_distribution[i_m, min_index] = np.int32(sorted_matrix[i_m, min_index])

        for k_z in range(Z):
            random_number = np.random.randint(scope_matrix[min_index * 2], scope_matrix[min_index * 2 + 1] + 1)
            second_line_individual[k_z, i_m] = random_number

    max_load = np.max(minimax_criterion_one_dimensional_array_loads)
    return minimax_criterion_number_distribution, max_load


def find_extra_minimax_criterion(barrier_min_numbers: int, second_line_individual, scope_matrix):
    index_for_extra, end_index_need_matrix = 0, 0

    extra_minimax_criterion_one_dimensional_array_loads = np.zeros(st.session_state["n"], dtype=np.int32)
    extra_minimax_criterion_additional_matrix = np.zeros((st.session_state["m"], st.session_state["n"]), dtype=np.int32)
    extra_minimax_criterion_number_distribution = np.zeros((st.session_state["m"], st.session_state["n"]),
                                                           dtype=np.int32)
    extra_minimax_criterion_i_m_indices = []

    for i_m in range(st.session_state["m"]):
        min_val = np.iinfo(np.int32).max
        min_index = -1

        for j_n in range(st.session_state["n"]):
            if st.session_state["sorted_matrix"][i_m, j_n] < min_val:
                min_val = st.session_state["sorted_matrix"][i_m, j_n]
                min_index = j_n

        traversal_variable_barrier = extra_minimax_criterion_one_dimensional_array_loads[min_index] + min_val

        if traversal_variable_barrier < barrier_min_numbers:
            extra_minimax_criterion_one_dimensional_array_loads[min_index] = traversal_variable_barrier
            extra_minimax_criterion_number_distribution[index_for_extra, min_index] = min_val

            for k_z in range(st.session_state["Z"]):
                random_number = np.random.randint(scope_matrix[min_index * 2], scope_matrix[min_index * 2 + 1] + 1)
                second_line_individual[k_z, i_m] = random_number

            index_for_extra += 1
        else:
            extra_minimax_criterion_additional_matrix[end_index_need_matrix] = st.session_state["sorted_matrix"][i_m]
            extra_minimax_criterion_i_m_indices.append(i_m)
            end_index_need_matrix += 1

    for i_e in range(end_index_need_matrix):
        min_val = np.iinfo(np.int32).max
        min_index = -1

        for j_n in range(st.session_state["n"]):
            traversal_variable = extra_minimax_criterion_additional_matrix[i_e, j_n] + \
                                 extra_minimax_criterion_one_dimensional_array_loads[j_n]
            if traversal_variable < min_val:
                min_val = traversal_variable
                min_index = j_n

        extra_minimax_criterion_one_dimensional_array_loads[min_index] = min_val
        extra_minimax_criterion_number_distribution[index_for_extra, min_index] = \
            extra_minimax_criterion_additional_matrix[i_e, min_index]

        for k_z in range(st.session_state["Z"]):
            random_number = np.random.randint(scope_matrix[min_index * 2], scope_matrix[min_index * 2 + 1] + 1)
            second_line_individual[k_z, extra_minimax_criterion_i_m_indices[i_e]] = random_number

        index_for_extra += 1

    max_load = np.max(extra_minimax_criterion_one_dimensional_array_loads)
    return extra_minimax_criterion_number_distribution, max_load


def app():
    if "initialized" not in st.session_state:
        reset_state()

    def reset_on_change():
        keys_to_reset = [
            "original_matrix", "sorted_matrix", "scope_matrix", "m", "n", "Z", "Zi", "Pk", "Pm",
            "generation_method_selected", "generation_method", "algorithm_option", "progress_log", "distribution_file",
            "final_msg", "algorithm_completed", "zip_buffer", "finalize_called"
        ]
        for key in keys_to_reset:
            if key in st.session_state:
                del st.session_state[key]

    st.title('Алгоритм П-З + Голдберга')

    matrix_option = st.selectbox("Выберите источник матрицы:", ("", "Генерация новой матрицы", "Матрица задана"))

    if matrix_option == "Генерация новой матрицы":
        m = st.number_input("Количество строк матрицы (m)", min_value=1, value=20, step=1, on_change=reset_on_change)
        n = st.number_input("Количество столбцов матрицы (n)", min_value=1, value=3, step=1, on_change=reset_on_change)
        T1 = st.number_input("Начало диапазона значений (T1)", min_value=1, value=10, step=1, on_change=reset_on_change)
        T2 = st.number_input("Конец диапазона значений (T2)", min_value=T1 + 1, value=20, step=1,
                             on_change=reset_on_change)
        Z = st.number_input("Количество особей (Z)", min_value=100, value=1000, step=1, on_change=reset_on_change)
        Zi = st.number_input("Количество повторов (Zi)", min_value=1, value=100, step=1, on_change=reset_on_change)
        Pk = st.number_input("Вероятность кроссовера (Pk)", min_value=0.0, max_value=1.0, value=0.5, step=0.01,
                             on_change=reset_on_change)
        Pm = st.number_input("Вероятность мутации (Pm)", min_value=0.0, max_value=1.0, value=0.5, step=0.01,
                             on_change=reset_on_change)

        if st.button("Сгенерировать матрицу"):
            original_matrix, sorted_matrix = generate_matrix(m, n, T1, T2)
            scope_matrix = np.zeros(2 * n, dtype=np.int32)
            find_scope_matrix(n, scope_matrix)

            st.session_state["original_matrix"] = original_matrix
            st.session_state["sorted_matrix"] = sorted_matrix
            st.session_state["scope_matrix"] = scope_matrix
            st.session_state["m"] = m
            st.session_state["n"] = n
            st.session_state["Z"] = Z
            st.session_state["Zi"] = Zi
            st.session_state["Pk"] = Pk
            st.session_state["Pm"] = Pm
            st.session_state["generation_method_selected"] = True

    elif matrix_option == "Матрица задана":
        uploaded_file = st.file_uploader("Загрузите файл матрицы", type="txt")
        if uploaded_file is not None:
            matrix_of_all_iterations = process_uploaded_file(uploaded_file)
            if matrix_of_all_iterations is not None:
                m, n = matrix_of_all_iterations.shape
                sorted_matrix = np.array(sorted(matrix_of_all_iterations, key=lambda x: sum(x), reverse=True))
                scope_matrix = np.zeros(2 * n, dtype=np.int32)
                find_scope_matrix(n, scope_matrix)

                Z = st.number_input("Количество особей (Z)", min_value=100, value=1000, step=1,
                                    on_change=reset_on_change)
                Zi = st.number_input("Количество повторов (Zi)", min_value=1, value=100, step=1,
                                     on_change=reset_on_change)
                Pk = st.number_input("Вероятность кроссовера (Pk)", min_value=0.0, max_value=1.0, value=0.5, step=0.01,
                                     on_change=reset_on_change)
                Pm = st.number_input("Вероятность мутации (Pm)", min_value=0.0, max_value=1.0, value=0.5, step=0.01,
                                     on_change=reset_on_change)

                if st.button("Продолжить"):
                    st.session_state["original_matrix"] = matrix_of_all_iterations
                    st.session_state["sorted_matrix"] = sorted_matrix
                    st.session_state["scope_matrix"] = scope_matrix
                    st.session_state["m"] = m
                    st.session_state["n"] = n
                    st.session_state["Z"] = Z
                    st.session_state["Zi"] = Zi
                    st.session_state["Pk"] = Pk
                    st.session_state["Pm"] = Pm
                    st.session_state["generation_method_selected"] = True
        else:
            if "generation_method_selected" in st.session_state:
                del st.session_state["generation_method_selected"]
            if "generation_method" in st.session_state:
                del st.session_state["generation_method"]
            if "algorithm_option" in st.session_state:
                del st.session_state["algorithm_option"]

    if st.session_state.get("generation_method_selected"):
        generation_method = st.selectbox("Выберите метод генерации особей:",
                                         ("", "Рандомные особи", "Генерация особей с помощью Плотникова-Зверева"))

        if generation_method:
            st.session_state["generation_method"] = generation_method

        if st.session_state.get("generation_method") == "Рандомные особи":
            if st.button("Начать работу"):
                second_line_individual = np.random.randint(0, 256,
                                                           (st.session_state["Z"], st.session_state["m"])).astype(
                    np.int32)
                run_algorithm(st.session_state["Z"], st.session_state["m"], st.session_state["n"],
                              st.session_state["Zi"], st.session_state["Pk"], st.session_state["Pm"],
                              st.session_state["sorted_matrix"], st.session_state["scope_matrix"], 0,
                              second_line_individual)
                st.session_state["progress_log_Goldberg"].append("Алгоритм завершен.")
                st.experimental_rerun()

        elif st.session_state.get("generation_method") == "Генерация особей с помощью Плотникова-Зверева":
            algorithm_option = st.selectbox("Выберите алгоритм Плотникова-Зверева:",
                                            ("", "Алгоритм Плотникова-Зверева по минимаксному критерию",
                                             "Алгоритм Плотникова-Зверева по минимаксному критерию с барьером"))

            if algorithm_option:
                st.session_state["algorithm_option"] = algorithm_option

            if st.session_state.get("algorithm_option"):
                if st.session_state["algorithm_option"] == "Алгоритм Плотникова-Зверева по минимаксному критерию":
                    st.write("Выполнение алгоритма Плотникова-Зверева по минимаксному критерию...")

                    second_line_individual = np.zeros((st.session_state["Z"], st.session_state["m"]), dtype=np.int32)
                    minimax_criterion_number_distribution, max_load = find_minimax_criterion(
                        second_line_individual, st.session_state["scope_matrix"], st.session_state["sorted_matrix"],
                        st.session_state["m"], st.session_state["n"], st.session_state["Z"]
                    )

                    col1, col2, col3 = st.columns(3)

                    with col1:
                        st.write("Изначальная матрица")
                        st.write(st.session_state["original_matrix"])

                    with col2:
                        st.write("Отсортированная матрица")
                        st.write(st.session_state["sorted_matrix"])

                    with col3:
                        st.write("Распределение чисел")
                        st.write(minimax_criterion_number_distribution)

                    st.markdown(
                        f"<p style='font-size:24px; font-weight:bold; text-align: center;'>MAX из массива нагрузки: {max_load}</p>",
                        unsafe_allow_html=True
                    )

                    choice_slice_number = st.number_input("Введите % соотношение подмешивания особей", min_value=1,
                                                          max_value=99, value=10, key="choice_slice_number_minimax")

                    if st.button("Подтвердить и начать алгоритм"):
                        slice_number = int(st.session_state["Z"] * choice_slice_number / 100)
                        second_line_individual[:slice_number] = np.random.randint(0, 256, (
                            slice_number, st.session_state["m"])).astype(np.int32)

                        run_algorithm(st.session_state["Z"], st.session_state["m"], st.session_state["n"],
                                      st.session_state["Zi"], st.session_state["Pk"], st.session_state["Pm"],
                                      st.session_state["sorted_matrix"], st.session_state["scope_matrix"], 1,
                                      second_line_individual)

                elif st.session_state[
                    "algorithm_option"] == "Алгоритм Плотникова-Зверева по минимаксному критерию с барьером":
                    st.write("Выполнение алгоритма Плотникова-Зверева по минимаксному критерию с барьером...")

                    second_line_individual = np.zeros((st.session_state["Z"], st.session_state["m"]), dtype=np.int32)
                    min_elements = np.min(st.session_state["sorted_matrix"], axis=1)
                    sum_of_min_elements = (np.sum(min_elements) // st.session_state["n"]) + 1
                    extra_minimax_criterion_number_distribution, max_load = find_extra_minimax_criterion(
                        sum_of_min_elements, second_line_individual, st.session_state["scope_matrix"]
                    )

                    col1, col2, col3 = st.columns(3)

                    with col1:
                        st.write("Изначальная матрица")
                        st.write(st.session_state["original_matrix"])

                    with col2:
                        st.write("Отсортированная матрица")
                        st.write(st.session_state["sorted_matrix"])

                    with col3:
                        st.write("Распределение чисел")
                        st.write(extra_minimax_criterion_number_distribution)

                    st.markdown(
                        f"<p style='font-size:24px; font-weight:bold; text-align: center;'>MAX из массива нагрузки: {max_load}</p>",
                        unsafe_allow_html=True
                    )

                    choice_slice_number = st.number_input("Введите % соотношение подмешивания особей", min_value=1,
                                                          max_value=99, value=10,
                                                          key="choice_slice_number_extra_minimax")

                    if st.button("Подтвердить и начать алгоритм"):
                        slice_number = int(st.session_state["Z"] * choice_slice_number / 100)
                        second_line_individual[:slice_number] = np.random.randint(0, 256, (
                            slice_number, st.session_state["m"])).astype(np.int32)

                        run_algorithm(st.session_state["Z"], st.session_state["m"], st.session_state["n"],
                                      st.session_state["Zi"], st.session_state["Pk"], st.session_state["Pm"],
                                      st.session_state["sorted_matrix"], st.session_state["scope_matrix"], 1,
                                      second_line_individual)

    if st.session_state.get('algorithm_completed'):
        st.write(st.session_state['final_msg'])
        st.download_button(
            label="Скачать все файлы",
            data=st.session_state['zip_buffer'],
            file_name="results.zip",
            mime="application/zip"
        )


if __name__ == "__main__":
    app()
