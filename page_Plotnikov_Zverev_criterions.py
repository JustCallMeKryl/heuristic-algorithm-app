# pages/page_Plotnikov_Zverev_criterions.py
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


def reset_state():
    st.session_state.pop("original_matrix", None)
    st.session_state.pop("minimax_criterion_number_distribution", None)
    st.session_state.pop("max_load", None)


def app():
    st.title("Алгоритм Плотникова-Зверева")

    algorithm_option = st.selectbox(
        "Выберите алгоритм:",
        ("", "Алгоритм Плотникова-Зверева по минимаксному критерию",
         "Алгоритм Плотникова-Зверева по минимаксному критерию с барьером"),
        on_change=reset_state
    )

    if algorithm_option == "Алгоритм Плотникова-Зверева по минимаксному критерию":

        matrix_option = st.selectbox(
            "Выберите источник матрицы:",
            ("", "Генерация новой матрицы", "Матрица задана"),
            on_change=reset_state
        )

        if matrix_option == "Генерация новой матрицы":
            m = st.number_input("Количество строк матрицы", min_value=1, value=21, step=1)
            n = st.number_input("Количество столбцов матрицы", min_value=1, value=3, step=1)
            T1 = st.number_input("Начало диапазона значений", value=10)
            T2 = st.number_input("Конец диапазона значений", value=20)

            if st.button("Сгенерировать матрицу"):
                original_matrix = np.random.randint(T1, T2 + 1, size=(m, n)).astype(np.int32)
                minimax_criterion_number_distribution, max_load = find_minimax_criterion(m, n, original_matrix)

                st.session_state["original_matrix"] = original_matrix
                st.session_state["minimax_criterion_number_distribution"] = minimax_criterion_number_distribution
                st.session_state["max_load"] = max_load

        if "original_matrix" in st.session_state and "minimax_criterion_number_distribution" in st.session_state and "max_load" in st.session_state:
            original_matrix = st.session_state["original_matrix"]
            minimax_criterion_number_distribution = st.session_state["minimax_criterion_number_distribution"]
            max_load = st.session_state["max_load"]

            # Вывод двух матриц и max_load в одну строку
            col1, col2 = st.columns(2)

            with col1:
                st.write("Сгенерированная матрица:")
                st.write(original_matrix)

                # Конвертируем матрицу в строку для сохранения в файл
                matrix_str = '\n'.join(' '.join(map(str, row)) for row in original_matrix)

                # Создаем ссылку для скачивания файла
                st.download_button(
                    label="Сгенерированная матрица",
                    data=matrix_str,
                    file_name='generated_matrix.txt',
                    mime='text/plain'
                )

            with col2:
                st.write("Распределение чисел по приборам:")
                st.write(minimax_criterion_number_distribution)

                # Конвертируем матрицу в строку для сохранения в файл
                matrix_str = '\n'.join(' '.join(map(str, row)) for row in minimax_criterion_number_distribution)

                # Создаем ссылку для скачивания файла
                st.download_button(
                    label="Распределенная матрица",
                    data=matrix_str,
                    file_name='distributed_matrix.txt',
                    mime='text/plain'
                )

            st.write(f"MAX из массива нагрузки: {max_load}")

        elif matrix_option == "Матрица задана":
            uploaded_file = st.file_uploader("Загрузите файл матрицы", type="txt")

            if uploaded_file is not None:
                file_content = uploaded_file.getvalue().decode("utf-8").strip()

                if not file_content:
                    st.warning("Файл пустой. Пожалуйста, выберите другой файл.")
                else:
                    try:
                        original_matrix = np.loadtxt(uploaded_file, dtype=int)
                        minimax_criterion_number_distribution, max_load = find_minimax_criterion(
                            original_matrix.shape[0], original_matrix.shape[1], original_matrix)

                        # Сохраняем матрицы и max_load в session_state
                        st.session_state["original_matrix"] = original_matrix
                        st.session_state["minimax_criterion_number_distribution"] = minimax_criterion_number_distribution
                        st.session_state["max_load"] = max_load

                    except ValueError:
                        st.warning(
                            "Файл содержит недопустимые символы. Пожалуйста, загрузите файл, содержащий только числа.")

            if "original_matrix" in st.session_state and "minimax_criterion_number_distribution" in st.session_state and "max_load" in st.session_state:
                original_matrix = st.session_state["original_matrix"]
                minimax_criterion_number_distribution = st.session_state["minimax_criterion_number_distribution"]
                max_load = st.session_state["max_load"]

                # Вывод двух матриц и max_load в одну строку
                col1, col2 = st.columns(2)

                with col1:
                    st.write("Загруженная матрица:")
                    st.write(original_matrix)

                    # Конвертируем матрицу в строку для сохранения в файл
                    matrix_str = '\n'.join(' '.join(map(str, row)) for row in original_matrix)

                    # Создаем ссылку для скачивания файла
                    st.download_button(
                        label="Сгенерированная матрица",
                        data=matrix_str,
                        file_name='generated_matrix.txt',
                        mime='text/plain'
                    )

                with col2:
                    st.write("Распределение чисел по приборам:")
                    st.write(minimax_criterion_number_distribution)

                    # Конвертируем матрицу в строку для сохранения в файл
                    matrix_str = '\n'.join(' '.join(map(str, row)) for row in minimax_criterion_number_distribution)

                    # Создаем ссылку для скачивания файла
                    st.download_button(
                        label="Распределенная матрица",
                        data=matrix_str,
                        file_name='distributed_matrix.txt',
                        mime='text/plain'
                    )

                st.write(f"MAX из массива нагрузки: {max_load}")


if __name__ == "__main__":
    app()
