import streamlit as st
import numpy as np


def find_minimax_criterion(m_number, n_number, original_matrix_number):
    minimax_criterion_one_dimensional_array_loads = np.zeros(n_number, dtype=np.int32)  # одномерный массив нагрузки
    minimax_criterion_number_distribution = np.zeros((m_number, n_number),
                                                     dtype=np.int32)  # распределение чисел по приборам

    for i_m in range(m_number):
        min_val = np.iinfo(np.int32).max  # Используем максимально возможное значение для int32
        min_index = -1

        for j_n in range(n_number):
            # переменная обхода
            traversal_variable = original_matrix_number[i_m, j_n] + minimax_criterion_one_dimensional_array_loads[j_n]

            # Нахождение минимального элемента в строке и его индекса
            if traversal_variable < min_val:
                min_val, min_index = traversal_variable, j_n

        # Обновляем сумму минимальных элементов по их индексам в matrix_min_element_method
        minimax_criterion_one_dimensional_array_loads[min_index] = np.int32(min_val)
        # записываем числа для просмотра распределения
        minimax_criterion_number_distribution[i_m, min_index] = np.int32(original_matrix_number[i_m, min_index])

    max_load = np.max(minimax_criterion_one_dimensional_array_loads)
    return minimax_criterion_number_distribution, max_load


def find_extra_minimax_criterion(m, n, original_matrix, barrier_min_numbers):
    index_for_extra, end_index_need_matrix = 0, 0

    extra_minimax_criterion_one_dimensional_array_loads = np.zeros(n, dtype=np.int32)  # одномерный массив нагрузки
    extra_minimax_criterion_additional_matrix = np.zeros((m, n), dtype=np.int32)  # дополнительная матрица
    extra_minimax_criterion_number_distribution = np.zeros((m, n), dtype=np.int32)  # распределение чисел по приборам
    extra_minimax_criterion_i_m_indices = []  # сохраняем пропущенные индексы

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

        index_for_extra += 1

    max_load = np.max(extra_minimax_criterion_one_dimensional_array_loads)
    return extra_minimax_criterion_number_distribution, max_load


def reset_state():
    st.session_state.pop("original_matrix", None)
    st.session_state.pop("sorted_matrix", None)
    st.session_state.pop("minimax_criterion_number_distribution", None)
    st.session_state.pop("max_load", None)
    st.session_state.pop("uploaded_file_content", None)
    st.session_state.pop("uploaded_file_name", None)


def display_matrices():
    if "original_matrix" in st.session_state and "sorted_matrix" in st.session_state and "minimax_criterion_number_distribution" in st.session_state and "max_load" in st.session_state:
        original_matrix = st.session_state["original_matrix"]
        sorted_matrix = st.session_state["sorted_matrix"]
        minimax_criterion_number_distribution = st.session_state["minimax_criterion_number_distribution"]
        max_load = st.session_state["max_load"]

        col1, col2, col3 = st.columns(3)

        with col1:
            st.write("Сгенерированная матрица")
            st.write(original_matrix)

            matrix_str = '\n'.join(' '.join(map(str, row)) for row in original_matrix)
            st.download_button(
                label="Скачать",
                data=matrix_str,
                file_name='matrix.txt',
                mime='text/plain'
            )

        with col2:
            st.write("Отсортированная матрица")
            st.write(sorted_matrix)

            matrix_str = '\n'.join(' '.join(map(str, row)) for row in sorted_matrix)
            st.download_button(
                label="Скачать",
                data=matrix_str,
                file_name='sorted_matrix.txt',
                mime='text/plain'
            )

        with col3:
            st.write("Распределение по приборам")
            st.write(minimax_criterion_number_distribution)

            matrix_str = '\n'.join(' '.join(map(str, row)) for row in minimax_criterion_number_distribution)
            st.download_button(
                label="Скачать",
                data=matrix_str,
                file_name='distribution.txt',
                mime='text/plain'
            )

        st.markdown(
            f"<p style='font-size:24px; font-weight:bold; text-align: center;'>MAX из массива нагрузки: {max_load}</p>",
            unsafe_allow_html=True)


def process_uploaded_file(uploaded_file):
    if uploaded_file is not None:
        file_content = uploaded_file.getvalue().decode("utf-8").strip()
        if file_content:
            st.session_state["uploaded_file_content"] = file_content
            st.session_state["uploaded_file_name"] = uploaded_file.name
        else:
            st.warning("Файл пустой. Пожалуйста, выберите другой файл.")
            reset_state()
    else:
        if "uploaded_file_name" in st.session_state:
            reset_state()


def process_file_content(algorithm):
    if "uploaded_file_content" in st.session_state and "uploaded_file_name" in st.session_state:
        file_content = st.session_state["uploaded_file_content"]
        try:
            original_matrix = np.loadtxt(file_content.splitlines(), dtype=int)
            m, n = original_matrix.shape

            sorted_matrix = np.array(sorted(original_matrix, key=lambda x: sum(x), reverse=True))

            minimax_criterion_number_distribution = None
            max_load = None

            if algorithm == "minimax":
                minimax_criterion_number_distribution, max_load = find_minimax_criterion(m, n, sorted_matrix)
            elif algorithm == "extra_minimax":
                min_elements = np.min(sorted_matrix, axis=1)  # находим минимальные элементы в каждой строке
                sum_of_min_elements = (np.sum(min_elements) // n) + 1  # находим барьер
                minimax_criterion_number_distribution, max_load = find_extra_minimax_criterion(m, n, sorted_matrix,
                                                                                               sum_of_min_elements)

            if minimax_criterion_number_distribution is not None and max_load is not None:
                st.session_state["original_matrix"] = original_matrix
                st.session_state["sorted_matrix"] = sorted_matrix
                st.session_state["minimax_criterion_number_distribution"] = minimax_criterion_number_distribution
                st.session_state["max_load"] = max_load

        except ValueError:
            st.warning("Файл содержит недопустимые символы. Пожалуйста, загрузите файл, содержащий только числа.")
            reset_state()


def app():
    st.title("Алгоритм Плотникова-Зверева")

    algorithm_option = st.selectbox(
        "Выберите алгоритм:",
        ("", "Алгоритм Плотникова-Зверева по минимаксному критерию",
         "Алгоритм Плотникова-Зверева по минимаксному критерию с барьером"),
        key="algorithm_option",
        on_change=reset_state
    )

    if algorithm_option == "Алгоритм Плотникова-Зверева по минимаксному критерию":
        matrix_option = st.selectbox(
            "Выберите источник матрицы:",
            ("", "Генерация новой матрицы", "Матрица задана"),
            key="matrix_option_minimax",  # уникальный ключ для каждого алгоритма
            on_change=reset_state
        )

        if matrix_option == "Генерация новой матрицы":
            m = st.number_input("Количество строк матрицы", min_value=1, value=21, step=1)
            n = st.number_input("Количество столбцов матрицы", min_value=1, value=3, step=1)
            T1 = st.number_input("Начало диапазона значений", value=10)
            T2 = st.number_input("Конец диапазона значений", value=20)

            if st.button("Сгенерировать матрицу"):
                original_matrix = np.random.randint(T1, T2 + 1, size=(m, n)).astype(np.int32)
                sorted_matrix = np.array(sorted(original_matrix, key=lambda x: sum(x), reverse=True))
                minimax_criterion_number_distribution, max_load = find_minimax_criterion(m, n, sorted_matrix)

                st.session_state["original_matrix"] = original_matrix
                st.session_state["sorted_matrix"] = sorted_matrix
                st.session_state["minimax_criterion_number_distribution"] = minimax_criterion_number_distribution
                st.session_state["max_load"] = max_load

            display_matrices()

        elif matrix_option == "Матрица задана":
            uploaded_file = st.file_uploader("Загрузите файл матрицы", type="txt", key="file_uploader_minimax")

            process_uploaded_file(uploaded_file)
            process_file_content("minimax")
            display_matrices()

    elif algorithm_option == "Алгоритм Плотникова-Зверева по минимаксному критерию с барьером":
        matrix_option = st.selectbox(
            "Выберите источник матрицы:",
            ("", "Генерация новой матрицы", "Матрица задана"),
            key="matrix_option_extra_minimax",  # уникальный ключ для каждого алгоритма
            on_change=reset_state
        )

        if matrix_option == "Генерация новой матрицы":
            m = st.number_input("Количество строк матрицы", min_value=1, value=21, step=1)
            n = st.number_input("Количество столбцов матрицы", min_value=1, value=3, step=1)
            T1 = st.number_input("Начало диапазона значений", value=10)
            T2 = st.number_input("Конец диапазона значений", value=20)

            if st.button("Сгенерировать матрицу"):
                original_matrix = np.random.randint(T1, T2 + 1, size=(m, n)).astype(np.int32)
                sorted_matrix = np.array(sorted(original_matrix, key=lambda x: sum(x), reverse=True))
                min_elements = np.min(sorted_matrix, axis=1)  # находим минимальные элементы в каждой строке
                sum_of_min_elements = (np.sum(min_elements) // n) + 1  # находим барьер
                minimax_criterion_number_distribution, max_load = find_extra_minimax_criterion(m, n, sorted_matrix,
                                                                                               sum_of_min_elements)

                st.session_state["original_matrix"] = original_matrix
                st.session_state["sorted_matrix"] = sorted_matrix
                st.session_state["minimax_criterion_number_distribution"] = minimax_criterion_number_distribution
                st.session_state["max_load"] = max_load

            display_matrices()

        elif matrix_option == "Матрица задана":
            uploaded_file = st.file_uploader("Загрузите файл матрицы", type="txt", key="file_uploader_extra_minimax")

            process_uploaded_file(uploaded_file)
            process_file_content("extra_minimax")
            display_matrices()


if __name__ == "__main__":
    app()
