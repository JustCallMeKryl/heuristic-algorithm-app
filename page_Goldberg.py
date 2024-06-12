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


def initialize_algorithm(Z, m, n):
    second_line_individual = np.random.randint(0, 256, (Z, m)).astype(np.int32)
    two_individual_second_line = np.zeros((2, m), dtype=np.int32)
    second_line_individual_for_rounds = np.zeros((Z, m), dtype=np.int32)
    return second_line_individual, two_individual_second_line, second_line_individual_for_rounds


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


def finalize_tracking(progress_log_Goldberg, start_time, end_time, second_line_individual, m, n, scope_matrix, sorted_matrix, individual_phenotypes):
    elapsed_time = end_time - start_time
    hours, rem = divmod(elapsed_time, 3600)
    minutes, seconds = divmod(rem, 60)
    final_msg = f"\n{int(hours):02}:{int(minutes):02}:{seconds:05.2f}"
    progress_log_Goldberg.append(final_msg)

    progress_distribution_file = kernel_phenotype_detection_and_save(len(second_line_individual), m, n, scope_matrix, second_line_individual, sorted_matrix, individual_phenotypes)

    st.write(final_msg)

    progress_Goldberg_file = "\n".join(progress_log_Goldberg)

    zip_buffer = download_files_as_zip(progress_Goldberg_file, progress_distribution_file)

    st.download_button(
        label="Скачать все файлы",
        data=zip_buffer,
        file_name="results.zip",
        mime="application/zip"
    )


def run_algorithm(Z, m, n, Zi, Pk, Pm, sorted_matrix, scope_matrix):
    second_line_individual, two_individual_second_line, second_line_individual_for_rounds = initialize_algorithm(Z, m,
                                                                                                                 n)
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

        kernel_phenotype_detection(Z, m, n, scope_matrix, second_line_individual, sorted_matrix,
                                   individual_phenotypes)
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


def display_matrix(matrix, title):
    st.write(title)
    st.write(matrix)


def download_files_as_zip(log_file_content, distribution_file_content):
    zip_buffer = io.BytesIO()
    with zipfile.ZipFile(zip_buffer, "a", zipfile.ZIP_DEFLATED) as zip_file:
        zip_file.writestr("progress_log.txt", log_file_content)
        zip_file.writestr("distribution_list.txt", distribution_file_content)
    zip_buffer.seek(0)
    return zip_buffer


def app():
    st.title('Алгоритм Голдберга')

    matrix_option = st.selectbox("Выберите источник матрицы:", ("", "Генерация новой матрицы", "Матрица задана"),
                                 key="matrix_option")

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
            original_matrix, sorted_matrix = generate_matrix(m, n, T1, T2)

            m, n = original_matrix.shape
            if m <= 30 and n <= 5:
                col1, col2 = st.columns(2)
                with col1:
                    display_matrix(original_matrix, "Сгенерированная матрица:")
                with col2:
                    display_matrix(sorted_matrix, "Отсортированная матрица:")

            scope_matrix = np.zeros(2 * n, dtype=np.int32)
            find_scope_matrix(n, scope_matrix)

            run_algorithm(Z, m, n, Zi, Pk, Pm, sorted_matrix, scope_matrix)

    elif matrix_option == "Матрица задана":
        uploaded_file = st.file_uploader("Загрузите файл матрицы", type="txt", key="file_uploader")
        matrix_of_all_iterations = process_uploaded_file(uploaded_file)

        if matrix_of_all_iterations is not None:
            m, n = matrix_of_all_iterations.shape
            sorted_matrix = np.array(sorted(matrix_of_all_iterations, key=lambda x: sum(x), reverse=True))

            if m <= 30 and n <= 5:

                col1, col2 = st.columns(2)
                with col1:
                    display_matrix(matrix_of_all_iterations, "Загруженная матрица:")
                with col2:
                    display_matrix(sorted_matrix, "Отсортированная матрица:")

            Z = st.number_input("Количество особей (Z)", min_value=100, value=1000, step=1)
            Zi = st.number_input("Количество повторов (Zi)", min_value=1, value=100, step=1)
            Pk = st.number_input("Вероятность кроссовера (Pk)", min_value=0.0, max_value=1.0, value=0.5, step=0.01)
            Pm = st.number_input("Вероятность мутации (Pm)", min_value=0.0, max_value=1.0, value=0.5, step=0.01)

            if st.button("Начать процесс"):
                scope_matrix = np.zeros(2 * n, dtype=np.int32)
                find_scope_matrix(n, scope_matrix)

                run_algorithm(Z, m, n, Zi, Pk, Pm, sorted_matrix, scope_matrix)


if __name__ == "__main__":
    app()
