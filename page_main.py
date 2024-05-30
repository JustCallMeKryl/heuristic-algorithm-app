import streamlit as st
from streamlit_option_menu import option_menu

import page_Plotnikov_Zverev_criterions
import page_Goldberg
import page_Plotnikov_Zverev_criterions_plus_Goldberg


# Функция для сброса состояния
def reset_state():
    for key in st.session_state.keys():
        del st.session_state[key]


def main():
    st.sidebar.title("Навигация")

    # Получаем текущий выбор, чтобы сохранить его после сброса состояния
    current_selection = st.session_state.get("current_selection", "Главная")

    with st.sidebar:
        selected = option_menu(
            menu_title="",
            options=["Главная", "Алгоритм Плотникова-Зверева", "Алгоритм Голдберга",
                     "Алгоритм Плотникова-Зверева + Голдберг"],
            icons=["house", "graph-up", "book", "bar-chart"],
            menu_icon="cast",
            default_index=0,
            key="menu",  # Указываем ключ для сохранения выбора
        )

    # Если выбор изменился, сбрасываем состояние
    if selected != current_selection:
        reset_state()
        st.session_state["current_selection"] = selected  # Обновляем выбор после сброса состояния

    if selected == "Главная":
        st.title('Главная страница')
    elif selected == "Алгоритм Плотникова-Зверева":
        page_Plotnikov_Zverev_criterions.app()
    elif selected == "Алгоритм Голдберга":
        page_Goldberg.app()
    elif selected == "Алгоритм Плотникова-Зверева + Голдберг":
        page_Plotnikov_Zverev_criterions_plus_Goldberg.app()


if __name__ == "__main__":
    main()
