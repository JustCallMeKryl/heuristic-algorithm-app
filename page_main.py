import streamlit as st
from streamlit_option_menu import option_menu

import page_Plotnikov_Zverev_criterions
import page_Goldberg
import page_Plotnikov_Zverev_criterions_plus_Goldberg
import page_analysis


def main():
    st.sidebar.title("Навигация")

    with st.sidebar:
        selected = option_menu(
            menu_title="",
            options=["Главная", "Алгоритм Плотникова-Зверева", "Алгоритм Голдберга",
                     "Алгоритм Плотникова-Зверева + Голдберг", "Анализ"],
            icons=["house", "graph-up", "book", "bar-chart", "search"],
            menu_icon="cast",
            default_index=0,
        )

    if selected == "Главная":
        st.title('Главная страница')
    elif selected == "Алгоритм Плотникова-Зверева":
        page_Plotnikov_Zverev_criterions.app()
    elif selected == "Алгоритм Голдберга":
        page_Goldberg.app()
    elif selected == "Алгоритм Плотникова-Зверева + Голдберг":
        page_Plotnikov_Zverev_criterions_plus_Goldberg.app()
    elif selected == "Анализ":
        page_analysis.app()


if __name__ == "__main__":
    main()
